library(Seurat)
library(harmony)
library(ggplot2)
library(stringr)
library(future)
library(Matrix)
library(data.table)
library(tidyverse)
library(ggpubr)
#library(DoubletFinder)
library(ROGUE)

#******LINUX******

#并行处理
plan("multicore", workers = 5)
options(future.globals.maxSize = 50000 * 1024^2)

#路径下文件夹以组织命名
tissue=dir("./")
#可能产生链接文件需要删除
samples=dir(list.files()) samples=samples[-19] samples=samples[-18]

#定义线粒体基因（这里是猕猴的）
mt.genes = c("ND1", "ND2",  "COX1", "COX2", "ATP8", "ATP6", "COX3","ND3",  "ND4L", "ND4",  "ND5",  "ND6",  "CYTB")

#遍历组织样本 分别创建seurat对象
for(t in tissue){
print(t)
s=dir(paste0("./",t,"/"))
print(s)
for(sample in s){
#根据样本名修改 跳过链接文件
if(str_detect(sample,'Frontal_lobe')){break}
if(str_detect(sample,'Spinal_cord')){break}
if(str_detect(sample,'Thalamus')){break}

#读取10X数据、创建Seurat对象
tmp= Read10X(paste0("./",t,"/",sample,"/"))
tmp= CreateSeuratObject(counts = tmp,project=sample,min.cells = 3, min.features = 200)

#加入线粒体百分比、分组、样本等信息
tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, features = intersect(mt.genes, rownames(tmp)))
tmp$tissue = t
tmp$sample = sample

#根据样本分组修改
tmp$group = ifelse(str_detect(sample, 'YM'), 'YM', 
                     ifelse(str_detect(sample, 'OM'), 'OM'))

p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#质控结果输出文件路径
png(paste0('../01_QC/', tissue, '_', sample, '_qc.png'), height = 500, width = 1500, res=100)
print(p)
dev.off()

#根据质控结果图修改参数
tmp=subset(tmp, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & nCount_RNA<10000)

#样本重命名+赋值给sample保存的“”里的变量,相当于创建了样本名字的变量
tmp <- RenameCells(tmp, add.cell.id = paste0(tissue, '_', sample))
assign(sample, tmp)
}
}


samples=dir(list.files()) 
#测试数据有两个链接名需要删除
samples=samples[-19] 
samples=samples[-18]

#SCTransform包括了 NormalizeData  FindVariableFeatures ScaleData
for (i in samples){
  #get可以取出assign赋给变量的值
  tmp <- get(i)
  tmp <- SCTransform(tmp, verbose = FALSE, vars.to.regress = 'percent.mt')
  assign(i, tmp)
}

#构造列表合并所有样本变量（seurat对象）
int.list = list()
for (tmp in samples){
  #append向列表添加对象，注意这里添加的是get(tmp)，也就是seurat对象值
  int.list = append(int.list, get(tmp))
}

#names之后，print的int.list会带两行，上行是names，下行是值
names(int.list) = samples

#存储SCTransform结果
saveRDS(int.list,"../00_RDS/int.list.SCTransform.rds")


#开始整合所有样本数据
#int.list=readRDS("../test2_RDS/int.list.SCTransform.rds")

#nfeatures设置保留top前多少个特征基因
int.features = SelectIntegrationFeatures(object.list = int.list, nfeatures = 3000)

int.list = PrepSCTIntegration(object.list = int.list, anchor.features = int.features, verbose = FALSE)

int.list <- lapply(X = int.list, FUN = function(x) {   
    x <- RunPCA(x, features = int.features, verbose = FALSE)
})
#method是SCT
int.anchors = FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT",anchor.features = int.features, verbose = FALSE,reduction = "rpca")

int.obj = IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

saveRDS(int.obj,"../00_RDS/int.obj.IntegrateData.rds")

#findmarkers时要改成rna
DefaultAssay(int.obj) <- "integrated"


#降维
int.obj = RunPCA(int.obj, verbose = FALSE)
p = ElbowPlot(int.obj, ndims=50)
png('../02_PCA/pca_elbow.png', height=500, width=600)
print(p)
dev.off()
#维数设置参考pca_elbow图
int.obj = RunUMAP(int.obj, reduction = "pca", dims = 1:30)
int.obj = FindNeighbors(int.obj, reduction = "pca", dims = 1:30)
int.obj = FindClusters(int.obj, resolution = 0.5)

# Visualization 各种可视化模式自选
p1 <- DimPlot(int.obj, reduction = "umap", split.by="group",group.by = "tissue")
png('../03_UMAP/umap1.png', height=500, width=1200)
print(p1)
dev.off()

p2 <- DimPlot(int.obj, reduction = "umap", group.by = "tissue")
png('../03_UMAP/umap2.png', height=500, width=600)
print(p2)
dev.off()

p3 <- DimPlot(int.obj, reduction = "umap", label = TRUE, repel = TRUE)
png('../03_UMAP/umap3.png', height=500, width=600)
print(p3)
dev.off()

p4 <- DimPlot(int.obj, reduction = "umap", split.by="group", label = TRUE,repel = TRUE)
png('../03_UMAP/umap4.png', height=500, width=1200)
print(p4)
dev.off()

#测试是否消除了批次效应
p5 <- DimPlot(int.obj, reduction = "umap", split.by="group", group.by="sample",label = TRUE,repel = TRUE)
png('../03_UMAP/umap5.png', height=500, width=1200)
print(p5)
dev.off()

saveRDS(int.obj,"../00_RDS/int.obj.Original.rds")


#markers
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

#取特定cluster子集的top10 markers，注意是“==”，而不是“=”
top10_clusterxx = subset(top10,cluster==9)

saveRDS(int.obj.markers,"../00_RDS/int.obj.Markers.rds")
write.csv(top10, "../04_markers/top10_markers.csv")


#cell type anno
a1.new <- DimPlot(anno, reduction = "umap", label = TRUE, repel = TRUE)
png('../test2_anno/anno_subset1.png'), height=500, width=600)
print(a1.new)
dev.off()

a2.new <- DimPlot(anno, reduction = "umap", split.by="group", label = TRUE,repel = TRUE)
png('../test2_anno/anno_subset2.png'), height=500, width=1200)
print(a2.new)
dev.off()


###差异分析
aging=subset(sce, seurat_clusters==9)
#这里做的是对aging这个单细胞亚群，做的年老和年轻的差异分析
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