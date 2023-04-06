x.paralle <- T
if (x.paralle)  {
  library(future)
  options(future.globals.maxSize = 120 * 1024^3)
  plan("multiprocess", workers = 4)
  plan()
}
library(xlsx)
library(Seurat)
library(stringr)
library(RColorBrewer)
library(plyr)
library(ggplot2)
library(Cairo)
library(dplyr)
library(data.table)
library(doMC)

rds.dir = file.path("/seurat_object_dir")
rds.file = file.path(rds.dir, "Seurat.celltype.RDS")

rds.data = readRDS(rds.file)
##### Set output directory
out.dir = file.path("/seurat_output_dir")
dir.create(out.dir)
##### Add patient information
rds.data@meta.data$patient <- ifelse(grepl("OCC", rds.data@meta.data$orig.ident), 
                                     str_extract(string = rds.data@meta.data$orig.ident, pattern = "OCC-[0-9]+"), 
                                     rds.data@meta.data$orig.ident)
##### Add survival information
short.term.surv = paste0("OCC-", c("05", "07", "09"))
long.term.surv = paste0("OCC-", c("01", "03", "04", "06", "08", "10", "11"))
rds.data@meta.data$Survival = "Unknown"
rds.data@meta.data[which(rds.data@meta.data$patient %in% short.term.surv),"Survival"] = "Short-term"
rds.data@meta.data[which(rds.data@meta.data$patient %in% long.term.surv),"Survival"] = "Long-term"

##### Filtering of contaminated cells containing high levels of mitochondrial transcripts (>=10% of total UMI counts)
rds.data <- subset(x = rds.data, subset = percent.mt < 10)

DefaultAssay(rds.data) = "RNA"

##### Data normalization
rds.data <- NormalizeData(object = rds.data, normalization.method = "LogNormalize", scale.factor = 10000)


##### Detection of variable genes across the single cells
rds.data <- FindVariableFeatures(object = rds.data, selection.method = "vst")
png(filename = paste0(out.dir, "/FindVariableGenes.png"), width = 900, height = 900)
print(VariableFeaturePlot(object = rds.data))
dev.off()


##### Regressing out cell cycle scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
rds.data <- CellCycleScoring(rds.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)


#####   Scaling the data
rds.data <- ScaleData(object = rds.data, features = rownames(rds.data), vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

x.meta = rds.data@meta.data

##### Saving the result
saveRDS(object = rds.data, file = paste0(out.dir, "/Seurat.normalization.QC.RDS"))

##### Performing linear dimensional reduction
x.dims <- 50
rds.data <- RunPCA(object = rds.data, features = VariableFeatures(object = rds.data), npcs = x.dims)


x.resolution <- 1.0

x.genes <- c("PTPRC", "EPCAM", "FAP", "PECAM1" , "CD3D", "CD4", "CD8A", "CD14", "CD79A")
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))

##### Dimension check using parallel computing
registerDoMC(50)

x.outdir <- paste0(out.dir, "/checked_all_PCA")
if (!dir.exists(x.outdir)) { dir.create(x.outdir) }

registerDoMC(50)

foreach (i = 4:x.dims) %dopar% {
  temp.obj.filter <- FindNeighbors(object = rds.data, dims = 1:i, force.recalc = T)
  temp.obj.filter <- FindClusters(object = temp.obj.filter, resolution = x.resolution)
  
  temp.obj.filter <- RunUMAP(object = temp.obj.filter, dims = 1:i)
  
  temp.cluster.color <- x.getPalette(length(unique(temp.obj.filter@meta.data$seurat_clusters)))
  temp.sample.color <- x.getPalette(length(unique(temp.obj.filter@meta.data$orig.ident)))
  
  png(filename = paste(x.outdir,"/UMAP.cell.", i, ".png", sep=""),  width = 600, height = 600)
  print(DimPlot(object = temp.obj.filter, reduction = "umap", label = T, pt.size = 0.5, label.size = 6, cols = temp.cluster.color))
  dev.off()
  
  png(filename = paste(x.outdir,"/UMAP.group.", i, ".png", sep=""), width = 600, height = 600)
  print(DimPlot(object = temp.obj.filter, reduction = "umap", group.by = "orig.ident", pt.size = 0.5, label.size = 6, cols = temp.sample.color))
  dev.off()
  
  png(filename = paste(x.outdir,"/UMAP.genes.", i, ".png", sep=""), width = 1800, height = 1800)
  print(FeaturePlot(object = temp.obj.filter, features = x.genes, pt.size = 0.5, order = T))
  dev.off()
}

#####   Select a dimension
x.dims <- 34
x.resolution <- 1.0

x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))


##### FindClusters
temp.obj.filter <- FindNeighbors(object = rds.data, dims = 1:x.dims, force.recalc = T)
temp.obj.filter <- FindClusters(object = temp.obj.filter, resolution = x.resolution)

temp.obj.filter <- RunUMAP(object = temp.obj.filter, dims = 1:x.dims)
temp.obj.filter <- RunTSNE(object = temp.obj.filter, dims = 1:x.dims)

##### Plot UMAP
x.cluster.color <- x.getPalette(length(unique(temp.obj.filter@meta.data$seurat_clusters)))
x.sample.color <- x.getPalette(length(unique(temp.obj.filter@meta.data$orig.ident)))

png(filename = paste0(out.dir, "/UMAP.cell.png"),  width = 1300, height = 1200)
DimPlot(object = temp.obj.filter, reduction = "umap", label = T, pt.size = 0.5, label.size = 10, cols = x.cluster.color)
dev.off()

png(filename = paste0(out.dir, "/UMAP.group.png"), width = 1400, height = 1200)
DimPlot(object = temp.obj.filter, reduction = "umap", group.by = "orig.ident", pt.size = 0.5, label.size = 10, cols = x.sample.color)
dev.off()

png(filename = paste0(out.dir, "/TSNE.cell.png"),  width = 1300, height = 1200)
DimPlot(object = temp.obj.filter, reduction = "tsne", label = T, pt.size = 0.5, label.size = 10, cols = x.cluster.color)
dev.off()

png(filename = paste0(out.dir, "/TSNE.group.png"), width = 1400, height = 1200)
DimPlot(object = temp.obj.filter, reduction = "tsne", group.by = "orig.ident", pt.size = 0.5, label.size = 10, cols = x.sample.color)
dev.off()


##### QC metric
temp.obj.filter@meta.data$log10.nCount_RNA <- log10(temp.obj.filter@meta.data$nCount_RNA)
temp.obj.filter@meta.data$log10.nFeature_RNA <- log10(temp.obj.filter@meta.data$nFeature_RNA)

png(filename = paste0(out.dir, "/UMAP.qc.png"), width = 1200, height = 1200)
FeaturePlot(object = temp.obj.filter, features = c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "S.Score"), order = T)
dev.off()

png(filename = paste0(out.dir, "/TSNE.qc.png"), width = 1200, height = 1200)
FeaturePlot(object = temp.obj.filter, features = c("log10.nCount_RNA", "log10.nFeature_RNA", "percent.mt", "S.Score"), order = T, reduction = "tsne")
dev.off()



##### Saving the Seurat object
saveRDS(object = temp.obj.filter, file = paste0(out.dir, "/Seurat.RDS"))
