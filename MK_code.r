library(Seurat)
library(stringr)
library(data.table)
library(tidyverse)
library(tibble)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(monocle3)
library(patchwork)
library(msigdbr)
library(fgsea)
library(ComplexHeatmap)

setwd("C:/Users/LIBRA/Desktop/R_workshop/temp/ITP_raw_data")
sample = dir("./")
mt.genes = c("MT-ND6","MT-CO2","MT-CYB","MT-ND2","MT-ND5","MT-CO1","MT-ND3","MT-ND4","MT-ND1","MT-ATP6","MT-CO3","MT-ND4L","MT-ATP8")
file_path = getwd()

options(future.globals.maxSize = 8000 * 1024^2)



###读取、构造、合并Seurat对象
for(s in sample){
    data=dir(paste0("./",s,"/"))
    for(d in data){
        if(str_detect(d,'tpm')){
            print(d)
            #data.table=F表明是data.frame格式，=T表明是data.table格式
            tmp= data.table::fread(paste0("./",s,"/",d), data.table = F)
            #去重第一列，空白第一列名会被赋值V1（注意是大写V）
            tmp = distinct(tmp,V1,.keep_all = T)
            #把V1即第一列对应的基因名作为数据的行名，否则后面seurat对象值全都是NA
            tmp = column_to_rownames(tmp,"V1") #或者用rownames(tmp) = tmp[,1]

            tmp= CreateSeuratObject(counts = tmp,project=sample,min.cells = 3, min.features = 200)
            #？为什么下面这行代码会报错：
            #Error in `[<-.data.frame`(`*tmp*`, , i, value = new("Seurat", assays = list( : 替换数据里有21081行，但数据有191
            #tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-", col.name = "percent.mt")
            tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, features = intersect(mt.genes, rownames(tmp)))
            #s这里不要写成sample啊啊啊
            tmp$sample = s
            tmp$group = strsplit(s,'_')[[1]][2]
            p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
            
            #自动创建目录
            folder_path = paste0(dirname(file_path),"/","01_QC/")
            if(!dir.exists(folder_path)){
                dir.create(folder_path)
            }
            png(paste0(folder_path, s, '_qc.png'), height = 500, width = 1500, res=100)
            print(p)
            dev.off()
            tmp=subset(tmp, subset = nFeature_RNA > 2000 & percent.mt < 1 & nCount_RNA > 10000)
        }
        
    }
    assign(s, tmp)
}

#？？？这里也许可以先合成int.list，再用lapply进行（normalize）、find variable features、scale、以及后面整合需要的pca
for (i in sample){
  tmp <- get(i)
  #tpm数据需要跳过NormalizeData，以免操作重复 
  tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2500)
  #注意这里只有feature genes会被scaled，要scale所有基因，参看官方教程
  #tmp <- ScaleData(tmp) 这里scale也可以，在后面整合的时候用lapply也可以
  assign(i, tmp)
}

int.list = list()
for (tmp in sample){
  int.list = append(int.list, get(tmp))
}

names(int.list) = sample

folder_path2 = paste0(dirname(file_path),"/","00_RDS/")
            if(!dir.exists(folder_path2)){
                dir.create(folder_path2)
            }
saveRDS(int.list,paste0(folder_path2, "int.list.Scaled.rds"))



###数据整合，消除批次效应
int.features = SelectIntegrationFeatures(object.list = int.list, nfeatures = 2500)
#int.list = PrepSCTIntegration(object.list = int.list, anchor.features = int.features, verbose = FALSE)
int.list <- lapply(X = int.list, FUN = function(x) {   
    x <- ScaleData(x, features = int.features, verbose = FALSE)
    x <- RunPCA(x, features = int.features, verbose = FALSE)
})
#删掉method=SCT
int.anchors = FindIntegrationAnchors(object.list = int.list, anchor.features = int.features, verbose = FALSE,reduction = "rpca")
int.obj = IntegrateData(anchorset = int.anchors, verbose = FALSE)

saveRDS(int.obj,paste0(folder_path2, "int.obj.IntegrateData.rds"))

DefaultAssay(int.obj) <- "integrated"



###降维
#注意这里要对整合后的feature genes再scale一遍，前面那个lapply里的scale只是对每一个的scale
#int.obj=readRDS(paste0(folder_path2, "int.obj.IntegrateData.rds"))
int.obj = ScaleData(int.obj, verbose = FALSE)
int.obj = RunPCA(int.obj, verbose = FALSE)
p = ElbowPlot(int.obj, ndims=50)

folder_path3 = paste0(dirname(file_path),"/","02_PCA/")
            if(!dir.exists(folder_path3)){
                dir.create(folder_path3)
            }
png(paste0(folder_path3,"pca_elbow.png"), height=500, width=600)
print(p)
dev.off()

int.obj = RunUMAP(int.obj, reduction = "pca", dims = 1:20)
int.obj = FindNeighbors(int.obj, reduction = "pca", dims = 1:20)
int.obj = FindClusters(int.obj, resolution = 0.05)

folder_path4 = paste0(dirname(file_path),"/","03_UMAP/")
            if(!dir.exists(folder_path4)){
                dir.create(folder_path4)
            }

p1 <- DimPlot(int.obj, reduction = "umap", label = TRUE, repel = TRUE)
png(paste0(folder_path4,"umap.png"), height=500, width=600)
print(p1)
dev.off()

p2 <- DimPlot(int.obj, reduction = "umap", group.by="group", label = TRUE, repel = TRUE)
png(paste0(folder_path4,"umap2.png"), height=500, width=600)
print(p2)
dev.off()

p3 <- DimPlot(int.obj, reduction = "umap", group.by="sample", label = TRUE, repel = TRUE)
png(paste0(folder_path4,"umap3.png"), height=500, width=600)
print(p3)
dev.off()

saveRDS(int.obj,paste0(folder_path2, "int.obj.Original.rds"))



###markers
DefaultAssay(int.obj) <- "RNA"
int.obj.markers <- FindAllMarkers(int.obj, only.pos = TRUE)
int.obj.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

int.obj.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

int.obj.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% -> ma

saveRDS(int.obj.markers,paste0(folder_path2, "int.obj.Markers.rds"))

folder_path5 = paste0(dirname(file_path),"/","04_Markers/")
            if(!dir.exists(folder_path5)){
                dir.create(folder_path5)
            }
write.csv(top10, paste0(folder_path5,"top10_markers.csv"))

#p6 <- DoHeatmap(int.obj, features = top10$gene) + NoLegend()
#DefaultAssay(int.obj) <- "RNA"
#如果前面没有先对int.obj scale，而是在整合的步骤scale，这里需要先对assay=RNA的int.obj scale一下
#int.obj <- ScaleData(int.obj)
p6 <- DoHeatmap(int.obj, features = ma$gene) + theme(axis.text.y = element_blank())#隐藏y轴基因名
pdf(paste0(folder_path5,"DoHeatmap.pdf"))
print(p6)
dev.off()



###cell type anno
new.cluster.ids <- c("Cycling_MK", "Immune_MK", "PLT_MK", "Niche_MK")
names(new.cluster.ids) <- levels(int.obj)
anno <- RenameIdents(int.obj, new.cluster.ids)
a1.new <- DimPlot(anno, reduction = "umap", label = TRUE, repel = TRUE)

folder_path6 = paste0(dirname(file_path),"/","05_Anno/")
            if(!dir.exists(folder_path6)){
                dir.create(folder_path6)
            }
png(paste0(folder_path6,"anno1.png"), height=500, width=600)
print(a1.new)
dev.off()

saveRDS(anno, paste0(folder_path2, "anno.rds"))

features = c("POLD2", "POLE2", "LIG1", "CCNE2", "CDK4", "CD226", "NFE2", 
"P2RY1", "P2RY12", "F2R", "HDGF", "BDNF", "GJA4", "TPM4", "ZYX", "TLR4", 
"ANXA1", "HLA-DRB1", "HLA-DMA", "CD74")
p7 <- DotPlot(anno, features = features) + RotatedAxis()
png(paste0(folder_path5,"Dotplot_anno.png"), height=800, width=800)
print(p7)
dev.off()

p8 <- VlnPlot(anno, features = c("PF4", "ITGA2B", "GP9", "GP1BA"), group.by = "group", pt.size = 0)
png(paste0(folder_path5,"MK_marker_vinplot.png"), height=500, width=600)
print(p8)
dev.off()






###题外 如果想，例如对上面p8小提琴图，加入统计学检验
#基因列表
genes <- c("PF4", "ITGA2B", "GP9", "GP1BA")

#创建一个空列表来存储图形
plot_list <- list()

#确定要比较的组对
comparisons <- list(c("Cycling_MK", "Immune_MK"), c("Cycling_MK", "PLT_MK"), c("Cycling_MK", "Niche_MK"), c("Immune_MK", "PLT_MK"), c("Immune_MK", "Niche_MK"), c("PLT_MK", "Niche_MK"))

#绘制每个基因的小提琴图并添加显著性统计
for (gene in genes) {
    vln_data <- as.data.frame(FetchData(anno, vars = gene))
    vln_data$cluster <- Idents(anno)
    
    p <- ggplot(vln_data, aes_string(x = "cluster", y = gene, fill = "cluster")) + 
        geom_violin(trim = FALSE) + 
        geom_boxplot(width = 0.1, fill = "white") + 
        theme_minimal() + 
        labs(y = gene, x = "Cluster", title = paste("Expression of", gene)) +
        theme(legend.position = "none")
    
    #添加显著性统计
    try({
        p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", hide.ns = FALSE)
    }, silent = TRUE)
    
    #使用assign创建变量并存储图形
    assign(paste0("p_", gene), p)
    
    #将图形添加到列表中
    plot_list[[gene]] <- get(paste0("p_", gene))
}

#将所有小提琴图组合在一张图上
combined_plot <- plot_grid(plotlist = plot_list, ncol = 1)

#绘制组合图形
print(combined_plot)






###monocle3拟时序分析
data<-GetAssayData(anno,assay ='RNA',slot ='counts')#slot='counts'（原始的）/'data'（归一化的）/'scale.data'（标准化的）
cell_metadata <-anno@meta.data
gene_annotation <-data.frame(gene_short_name =rownames(data))
rownames(gene_annotation)<-rownames(data)
cds <-new_cell_data_set(data,cell_metadata =cell_metadata,gene_metadata =gene_annotation)

cds@int_colData$reducedDims$PCA <- anno@reductions$pca@cell.embeddings
cds@int_colData$reducedDims$UMAP <- Embeddings(anno, reduction = "umap")
cds@clusters$UMAP$clusters <- Idents(anno)[rownames(colData(cds))]

#seurat的counts需要进行size factor计算（？data/scale.data）
cds <- estimate_size_factors(cds)

cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(cds)

#way1 根据umap进行分区
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves=FALSE, 
            label_branch_points=FALSE, graph_label_size=1.5)->pp1

#way2 合并分区
cds@clusters$UMAP$partitions <- factor(rep(1, ncol(cds)))
names(cds@clusters$UMAP$partitions) <- colnames(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves=FALSE, 
            label_branch_points=FALSE, graph_label_size=1.5)->pp2

folder_path7 = paste0(dirname(file_path),"/","07_Pseudotime/")
            if(!dir.exists(folder_path7)){
                dir.create(folder_path7)
            }
png(paste0(folder_path7,"pseudotime.png"), height=500, width=1200)
print(pp1+pp2)
dev.off()

saveRDS(cds, paste0(folder_path2, "pseudotime.cds.rds"))

genes <- c("CD34", "ITGA2B", "GP9", "GP1BA", "ITGB3", "GP5")
plot_genes_in_pseudotime(cds[genes,],
                         color_cells_by="celltype",
                         min_expr=0.5) -> pp3
png(paste0(folder_path7,"pseudo_marker.png"), height=1000, width=600)
print(pp3)
dev.off()






###添加倍性到meta.data
#自己写的版本：
j=c()
a=Cells(anno)#"D0ITP34_LY28_2N_1_1" ……
for(i in a){
    temp=strsplit(i,"_")[[1]][3]#2N……
    j=append(j,temp)
}
names(j)=Cells(anno)

#检查细胞名和倍性是不是匹配的
for(temp2 in a){
    if(strsplit(temp2,"_")[[1]][3]!=j[temp2]){
        print("error")
    }
}
#添加到meta.data
anno$ploid=j

#chatgpt优化版本：
#获取所有细胞的名称
cell_names <- Cells(anno)

#提取倍性信息
ploid <- sapply(cell_names, function(x) strsplit(x, "_")[[1]][3])

#检查提取的倍性信息是否正确（可选）
if (!all(sapply(cell_names, function(x) strsplit(x, "_")[[1]][3]) == ploid)) {
    stop("Mismatch in ploid information")
}

#将倍性信息添加到 meta.data 中
anno$ploid <- ploid
#或者：
anno <- AddMetaData(anno, metadata = ploid, col.name = "ploid")






###细胞比例统计堆叠图

##按照ploid统计细胞比例
cell.prop=prop.table(table(anno$celltype,anno$ploid),margin = 2)#margin=2表示按照列分组统计各细胞类型
Cellratio <- as.data.frame(cell.prop)#var1是celltype，var2是ploid
colnames(Cellratio) <- c("celltype", "ploid", "Freq")#对var命名，这样下面ggplot直接填对应名字就行了

#设置ploid和celltype的因子水平，以控制顺序
Cellratio$ploid <- factor(Cellratio$ploid, levels = c("2N", "4N", "8N", "16N", "32N"))
Cellratio$celltype <- factor(Cellratio$celltype, levels = c("Cycling_MK", "Niche_MK", "PLT_MK", "Immune_MK"))

#设置特定颜色
celltype_colors <- c("Cycling_MK" = "yellow", "Niche_MK" = "purple", "PLT_MK" = "red", "Immune_MK" = "green")

colourCount = length(unique(Cellratio$celltype))
library(ggplot2)
p <- ggplot(Cellratio) + 
    geom_bar(aes(x = ploid, y= Freq, fill = celltype), stat = "identity", position = "stack", width = 0.25, size = 0.5, colour = "#222222")+
    theme_classic() + #经典主题
    labs(x='Ploid', y = 'Ratio')+ #添加坐标轴名称
    scale_fill_manual(values = celltype_colors) +  # 设置特定颜色
    #coord_flip()+ #翻转坐标轴
    #theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid")) #添加外框 +
    scale_y_continuous(expand = c(0, 0)) #确保条形图紧贴 y 轴 0 值

print(p)

folder_path8 = paste0(dirname(file_path),"/","08_CellRatio/")
if(!dir.exists(folder_path8)){
    dir.create(folder_path8)
}
png(paste0(folder_path8,"CellRatio.png"), height=600, width=600)
print(p)
dev.off()



##按照样本统计细胞比例
#cell.prop2 <- prop.table(table(anno$celltype, anno$ploid, anno$group), margin = 2) X 没有按照不同分组统计
# 确保 meta_data 是一个数据框
meta_data <- as.data.frame(anno@meta.data)

# 将 group 和 ploid 列转换为因子类型
meta_data$group <- as.factor(meta_data$group)
meta_data$ploid <- as.factor(meta_data$ploid)

# 分开计算不同 group 下的 ploid 的 cell type 比例
cell.prop2 <- meta_data %>%
    group_by(group, ploid,celltype) %>%
    summarise(n = n()) %>%
    mutate(Freq = n / sum(n)) %>%
    ungroup()

print(cell.prop2)

Cellratio2 <- as.data.frame(cell.prop2)
colnames(Cellratio2) <- c("group", "ploid", "celltype", "ncount", "Freq")
#Cellratio2 %>% group_by(group,ploid) %>% summarise(sum(ncount))查看每个分组所有类型细胞总数
#Cellratio2 %>% group_by(group,ploid) %>% summarise(n=sum(ncount)) %>% mutate(freq=n/sum(n))查看每个分组不同倍性细胞占比
#由上发现这里3个分组（ITP和healthy）单细胞数据的倍性分布与流式统计不符，如“ITP Ⅱ的2N细胞极端增多”在这里是不符的
#这应该是由于在单细胞转录组建库过程中为独立分选不同倍性巨核细胞以及手工挑取单细胞建库
#且这里的目的是为了关联倍性和巨核细胞异质群，通过前面绘制的倍性UMAP将不同巨核异质群与倍性关联起来
#接着这里细胞比例统计，目的是看不同ITP/健康样本的不同倍性对应的巨核细胞异质群，同样是关联
#最后的下游分析，则结合流式结果（如ITP Ⅱ 2N细胞极端增加）和关联结果（2N细胞主要是Immune_MK和Cycling_MK），进行针对性的分析

Cellratio2$ploid <- factor(Cellratio2$ploid, levels = c("2N", "4N", "8N", "16N", "32N"))
Cellratio2$celltype <- factor(Cellratio2$celltype, levels = c("Cycling_MK", "Niche_MK", "PLT_MK", "Immune_MK"))
celltype_colors <- c("Cycling_MK" = "yellow", "Niche_MK" = "purple", "PLT_MK" = "red", "Immune_MK" = "green")

pp <- ggplot(Cellratio2, aes(x = ploid, y = Freq, fill = celltype)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.25, size = 0.5, colour = "#222222") + 
  theme_classic() +
  labs(x = 'Ploidy', y = 'Ratio') +
  scale_fill_manual(values = celltype_colors) +
  facet_wrap(~ group, ncol = 3) +  # 按样本分组展示 ncol表示分几列排布图表
  #coord_flip() +
  #theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid")) +
  scale_y_continuous(expand = c(0, 0))

print(pp)

file_path=getwd()
folder_path8 = paste0(dirname(file_path),"/","08_CellRatio/")
if(!dir.exists(folder_path8)){
    dir.create(folder_path8)
}
png(paste0(folder_path8,"CellRatio_groupby.png"), height=600, width=600)
print(pp)
dev.off()




####################下游分析###################
#################### ⅠITP PLT_MK ###################
setwd("../09_downstream/ⅠITP/PLT_MK/")

###GSEA分析
plt_mk_cells <- subset(anno, subset = celltype == "PLT_MK")

#获取MSigDB中的基因集
msigdbr_data <- msigdbr(species = "Homo sapiens", category = "C5")
gene_sets <- msigdbr_data %>%
    filter(gs_name %in% c("GOBP_HEMOSTASIS", "GOBP_PLATELET_ACTIVATION", "GOBP_PLATELET_AGGREGATION")) %>%
    split(x = .$gene_symbol, f = .$gs_name)

#重命名通路
names(gene_sets) <- c("Hemostasis", "Platelet_activation", "Platelet_aggregation")

#检查新的通路名称
print(names(gene_sets))


##***不区分样本的GSEA分析：

#获取差异基因
#min.pct和logfc.threshold:该基因至少在一个群体中有25%的细胞表达(默认10%)，logfc.threshold = 0 表示所有基因都将被考虑，而不管它们的表达变化有多小（默认值0.25表示筛选 绝对值>=0.25）
markers_ITP_I_vs_healthy <- FindMarkers(plt_mk_cells, ident.1 = "Ⅰtype", group.by = "group", ident.2 = "healthy", min.pct = 0.25, logfc.threshold = 0)
markers_ITP_II_vs_healthy <- FindMarkers(plt_mk_cells, ident.1 = "Ⅱtype", group.by = "group", ident.2 = "healthy", min.pct = 0.25, logfc.threshold = 0)

#准备排名数据
#way1
#Ⅰ
markers_ITP_I_vs_healthy <- markers_ITP_I_vs_healthy %>%
    rownames_to_column(var = "gene") #把行名转换成列

ranks_ITP_I_vs_healthy <- markers_ITP_I_vs_healthy %>% 
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    deframe() #转变成list才能进行fgsea
#之前的方法用names(ranks_ITP_I_vs_healthy)=rownames(markers_ITP_I_vs_healthy)出错，是因为后者没有降序，导致基因名不匹配

#Ⅱ
markers_ITP_II_vs_healthy <- markers_ITP_II_vs_healthy %>%
    rownames_to_column(var = "gene")

ranks_ITP_II_vs_healthy <- markers_ITP_II_vs_healthy %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    deframe()

#way2
#直接合并 以Ⅰ为例（如果不加dplyr::，dplyr::，tibble:: ，会报错）
ranks_ITP_I_vs_healthy <- markers_ITP_I_vs_healthy %>%
  rownames_to_column(var = "gene") %>%
  dplyr::arrange(desc(avg_log2FC)) %>%  # 确保使用dplyr中的arrange
  dplyr::select(gene, avg_log2FC) %>%  # 确保使用dplyr中的select
  tibble::deframe()

#进行GSEA分析
fgseaRes_ITP_I <- fgsea(pathways = gene_sets, stats = ranks_ITP_I_vs_healthy)
fgseaRes_ITP_II <- fgsea(pathways = gene_sets, stats = ranks_ITP_II_vs_healthy)

#绘制GSEA折线图
gsea_plot_1 <- plotEnrichment(gene_sets$Hemostasis, ranks_ITP_I_vs_healthy) + 
  labs(title = "GSEA Hemostasis ITP I vs Healthy") -> g1

gsea_plot_2 <- plotEnrichment(gene_sets$Platelet_activation, ranks_ITP_I_vs_healthy) + 
  labs(title = "GSEA Platelet Activation ITP I vs Healthy") -> g2

gsea_plot_3 <- plotEnrichment(gene_sets$Platelet_aggregation, ranks_ITP_I_vs_healthy) + 
  labs(title = "GSEA Platelet Aggregation ITP I vs Healthy") -> g3

gsea_plot_4 <- plotEnrichment(gene_sets$Hemostasis, ranks_ITP_II_vs_healthy) + 
  labs(title = "GSEA Hemostasis ITP II vs Healthy") -> g4

gsea_plot_5 <- plotEnrichment(gene_sets$Platelet_activation, ranks_ITP_II_vs_healthy) + 
  labs(title = "GSEA Platelet Activation ITP II vs Healthy") -> g5

gsea_plot_6 <- plotEnrichment(gene_sets$Platelet_aggregation, ranks_ITP_II_vs_healthy) + 
  labs(title = "GSEA Platelet Aggregation ITP II vs Healthy") -> g6

#组合折线图
combined_gsea_plots <- plot_grid(gsea_plot_1, gsea_plot_2, gsea_plot_3, gsea_plot_4, gsea_plot_5, gsea_plot_6, ncol = 2)
print(combined_gsea_plots)

file_path=getwd()
folder_path9 = paste0(file_path,"/","GSEA/")
if(!dir.exists(folder_path9)){
    dir.create(folder_path9)
}
png(paste0(folder_path9,"combined_gsea_plots.png"), height=2000, width=2000)
print(combined_gsea_plots)
dev.off()


##***区分样本的GSEA分析：

#构造一个nes提取函数
extract_nes <- function(fgsea_result) {
  fgsea_result %>%
    filter(pathway %in% c("Hemostasis", "Platelet_activation", "Platelet_aggregation")) %>%
    select(pathway, sample, NES) %>%
    pivot_wider(names_from = pathway, values_from = NES)
}

#提取plt_mk_cells中的sample 列
healthy_cell <- subset(plt_mk_cells, group == "healthy")
samples_h <- unique(healthy_cell$sample)

ITPcell_1 <- subset(plt_mk_cells, group == "Ⅰtype")
samples_1 <- unique(ITPcell_1$sample)

ITPcell_2 <- subset(plt_mk_cells, group == "Ⅱtype")
samples_2 <- unique(ITPcell_2$sample)


#创建一个空的数据框用于存储所有样本的 NES 值
nes_values <- data.frame()
# 创建一个空的列表用于存储每个样本的 NES 值
nes_list <- list()

#对ITP Ⅰ每个样本进行fgse分析
for (s in samples_1) {
  #******哭死了啊啊啊啊啊，不要写sample in samples_1 ……sample==sample，这样取出来的细胞根本没变化，换成s就好了！！！
  sample_cells <- subset(ITPcell_1, sample == s)#提取特定样本的数据
  combined_cells <- merge(sample_cells, y = healthy_cell)#***注意这里没有做批次效应消除

  markers_ITP_vs_healthy <- FindMarkers(combined_cells, ident.1 = "Ⅰtype", group.by = "group", ident.2 = "healthy", min.pct = 0.25, logfc.threshold = 0)
  ranks_ITP_vs_healthy <- markers_ITP_vs_healthy %>%
    rownames_to_column(var = "gene") %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    tibble::deframe()
  fgseaRes_ITP <- fgsea(pathways = gene_sets, stats = ranks_ITP_vs_healthy)
  
  #提取NES值并添加样本信息
  #先mutate添加样本列，否则构造函数取不到sample会报错
  #extract_nes()，而非extract_nes(fgseaRes_ITP)，才能正确传递管道参数
  nes_ITP <- fgseaRes_ITP %>% mutate(sample = s) %>% extract_nes() 

  #合并到nes_values数据框中
  #way1：循环内直接一个个bind合并进去
  #nes_values <- bind_rows(nes_values, nes_ITP) 
  #way2：循环内列表，循环外bind
  #way3：更加严谨，用assign防止变量被覆盖（实际上应该不会覆盖，都是前面sample惹的祸）

  #使用 assign 函数创建全局独立变量
  assign(paste0("nes_", s), nes_ITP, envir = .GlobalEnv)

  #将结果添加到列表中
  nes_list[[s]] <- get(paste0("nes_", s))
}

#对ITP Ⅱ每个样本进行fgse分析
for (s in samples_2) {
  #提取特定样本的数据
  sample_cells <- subset(ITPcell_2, sample == s)
  combined_cells <- merge(sample_cells, y = healthy_cell)
  markers_ITP_vs_healthy <- FindMarkers(combined_cells, ident.1 = "Ⅱtype", group.by = "group", ident.2 = "healthy", min.pct = 0.25, logfc.threshold = 0)
  ranks_ITP_vs_healthy <- markers_ITP_vs_healthy %>%
    rownames_to_column(var = "gene") %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    tibble::deframe()
  fgseaRes_ITP <- fgsea(pathways = gene_sets, stats = ranks_ITP_vs_healthy)
  nes_ITP <- fgseaRes_ITP %>% mutate(sample = s) %>% extract_nes() 
  #nes_values <- bind_rows(nes_values, nes_ITP)
  assign(paste0("nes_", s), nes_ITP, envir = .GlobalEnv)
  nes_list[[s]] <- get(paste0("nes_", s))
}
# 将列表中的所有元素合并为一个数据框
nes_values <- do.call(bind_rows, nes_list)

#转换数据框格式适应热图绘制
nes_matrix <- nes_values %>%
  pivot_longer(cols = -sample, names_to = "pathway", values_to = "NES") %>%
  pivot_wider(names_from = sample, values_from = NES) %>%
  column_to_rownames(var = "pathway") %>%
  as.matrix()
#***示例理解上述代码：
#第一次：长格式->宽格式
#pivot_longer 函数用于将数据框从宽格式转换为长格式。
#cols = -sample 指定除了 sample 列以外的所有列都进行转换。
#names_to = "pathway" 将列名转换为新列 pathway 的值。
#values_to = "NES" 将原始列的值转换为新列 NES 的值。

#第一次转换前：
sample  Hemostasis  Platelet_activation  Platelet_aggregation
sample1 1.5         2.0                  1.8
sample2 1.3         1.9                  1.7

#第一次转换后&第二次转换前：
sample  pathway              NES
sample1 Hemostasis           1.5
sample1 Platelet_activation  2.0
sample1 Platelet_aggregation 1.8
sample2 Hemostasis           1.3
sample2 Platelet_activation  1.9
sample2 Platelet_aggregation 1.7

#第二次：宽格式->长格式
#pivot_wider 函数用于将数据框从长格式转换为宽格式。
#names_from = sample 指定 sample 列的值作为新的列名。
#values_from = NES 指定 NES 列的值作为新的列值。

#第二次转换后：
pathway              sample1 sample2
Hemostasis           1.5     1.3
Platelet_activation  2.0     1.9
Platelet_aggregation 1.8     1.7

#column_to_rownames 函数将指定的列 pathway 设置为行名。
#这样做可以使矩阵的行名表示路径名称，方便热图绘制。
                    sample1 sample2
Hemostasis              1.5     1.3
Platelet_activation     2.0     1.9
Platelet_aggregation    1.8     1.7

#as.matrix 函数将数据框转换为矩阵。
#矩阵格式是绘制热图时所需的数据格式。

#创建NES热图

#设置column_split，以实现分组展示热图
samples=c(samples_1, samples_2)

#way1
#way1
has_Ⅰ_type <- grepl("_Ⅰtype",samples)
has_Ⅱ_type <- grepl("_Ⅱtype", samples)
sample_types1 <- factor(ifelse(has_Ⅰ_type, "ITP I", ifelse(has_Ⅱ_type, "ITP II", "Unknown")))

#way2
sample_groups <- sapply(samples, function(sample) {
  if (grepl("_Ⅰtype", sample)) {
    return("ITP I")
  } else if (grepl("_Ⅱtype", sample)) {
    return("ITP II")
  } else {
    return("Unknown")
  }
})
sample_types2 <- factor(sample_groups, levels = c("ITP I", "ITP II"))

heatmap <- Heatmap(
  nes_matrix,
  name = "NES",
  row_names_side = "right",#left
  column_names_side = "bottom",#top
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = sample_types1,
  row_split = factor(rep("Pathway", 3)),
  column_title = "",
  column_title_side = "top",
  row_title = "",
  row_title_side = "left",
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = NA))
  }
)

#绘制热图
draw(heatmap)

png("../PLT_MK/GSEA/NES_heatmap",width=600,height=600)
print(heatmap)
dev.off()




###差异分析
aging=subset(sce, seurat_clusters==9)
om_deg=FindMarkers(aging, ident.1 = "OM", group.by = 'group', ident2 = "YM")
head(om_deg)

degdf=om_deg
degdf$symbol <- rownames(om_deg)
logFC_t=0

P.Value_t = 0.05
degdf$change = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < 0,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 0,"up","stable"))
ggplot(degdf, aes(avg_log2FC,  -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()



Singlecellratio_plotstat(anno, group_by = "group",
                         meta.include = c("group","orig.ident"),
                         color_by = "cell.type")
