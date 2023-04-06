library(Seurat)
library(monocle)
library(dplyr)
library(RColorBrewer)
library(doMC)
library(data.table)
library(ggplot2)
library(enrichR)
library(Cairo)
library(ggridges)
options(bitmapType='cairo')
##### Set monocle2 output directory
monocle.outdir <- file.path("/monocle2_output_dir")
##### Read Seurat object
seurat.outdir <- file.path("/seurat_output_dir")
x.file <- paste0(seurat.outdir, "/Seurat.celltype.RDS")
x.seurat.obj <- readRDS(x.file)
##### Find Markers of all clusters
x.marker.df <- FindAllMarkers(object = x.seurat.obj, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", assay="RNA", slot="data")
##### Save the marker file
write.table(x = x.marker.df, file = paste0(seurat.outdir,"/Seurat.cluster.markers.txt"), quote = F, sep = "\t", row.names = F)

##### Importing data with Seurat
### a numeric matrix of expression values, where rows are genes, and columns are cells
x.data <- as(as.matrix(x.seurat.obj@assays$RNA@data), 'sparseMatrix')
x.data <- x.data[rowSums(x.data) > 0,]

##### Make a phenoData
x.pd <- new('AnnotatedDataFrame', data = x.seurat.obj@meta.data)
##### Make a featureData
x.fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(x.data), row.names = row.names(x.data)))

x.cds <- newCellDataSet(x.data, phenoData = x.pd, featureData = x.fd, expressionFamily = negbinomial.size())


##### Estimate size factors and dispersions
x.cds <- estimateSizeFactors(x.cds)
x.cds <- estimateDispersions(x.cds)


##### Filter off low-quality cells 
set.seed(123)

x.cds <- detectGenes(x.cds, min_expr = 0.1)

##### Select the top 50 most differentially expressed genes for each cluster
x.marker.df <- read.table(x.deg.file, header = T, sep = "\t")
##### select cluster information (Ex: C0: Epithelial --> 0)
x.marker.df$cluster = stringr::str_extract(x.marker.df$cluster, pattern = "C[0-9]+")
x.marker.df$cluster = stringr::str_remove(x.marker.df$cluster, pattern = "^C")
x.marker.df$cluster = factor(x.marker.df$cluster, levels = unique(x.marker.df$cluster))
x.marker.df <- x.marker.df[x.marker.df$cluster %in% as.character(unique(x.seurat.obj@meta.data$seurat_clusters)),]
x.marker.df <- subset(x.marker.df, p_val_adj < 0.05)
x.marker.top.df <- x.marker.df %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
x.order.genes <- unique(as.character(x.marker.top.df$gene))
x.cds <- setOrderingFilter(x.cds, x.order.genes)

#####	Run Monocle2 and make a plot
##### Set color
x.getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
x.cluster.color <- x.getPalette(length(unique(x.seurat.obj@meta.data$seurat_clusters)))
x.sample.color <- x.getPalette(length(unique(x.seurat.obj@meta.data$orig.ident)))
x.type.color = x.getPalette2(length(levels(x.seurat.obj@meta.data$ClusterLabel)))
names(x.type.color) = levels(x.seurat.obj@meta.data$ClusterLabel)

##### Run monocle2 for multiple components
foreach (i = 2:10) %dopar%  {
  sprintf("Component: %d", i)
  temp.outdir <- paste0(monocle.outdir, "/max_components_", i)
  
  if (!dir.exists(temp.outdir)) { dir.create(temp.outdir) }
  
  temp.cds <- reduceDimension(cds = x.cds, max_components = i, method = 'DDRTree')
  temp.cds <- orderCells(temp.cds)
  
  
  ##### Plot trajectory
  png(filename = paste0(temp.outdir, "/monocle.samples.png"), width = 1000, height = 900)
  print(plot_cell_trajectory(temp.cds, color_by = "orig.ident", cell_size = 3) + scale_color_manual(values = x.sample.color) + theme_classic(base_size = 20))
  dev.off()
  
  png(filename = paste0(temp.outdir, "/monocle.clusters.png"), width = 1000, height = 900)
  print(plot_cell_trajectory(temp.cds, color_by = "seurat_clusters", cell_size = 3) + scale_color_manual(values = x.cluster.color) + theme_classic(base_size = 20))
  dev.off()
  
  png(filename = paste0(temp.outdir, "/monocle.clusters_by_samples.png"), width = 2100, height = 900)
  print(plot_cell_trajectory(temp.cds, color_by = "seurat_clusters", cell_size = 0.1, show_branch_points = FALSE) + facet_wrap(~orig.ident, nrow = 2) + scale_color_manual(values = x.cluster.color) + theme_classic(base_size = 20))
  dev.off()
  
  png(filename = paste0(temp.outdir, "/monocle.pseudotime.png"), width = 900, height = 900)
  print(plot_cell_trajectory(temp.cds, color_by = "Pseudotime", cell_size = 3, show_branch_points = FALSE) + theme_classic(base_size = 20))
  dev.off()
  
  png(filename = paste0(temp.outdir, "/monocle.celltypes.png"), width = 900, height = 900)
  print(plot_cell_trajectory(temp.cds, color_by = "ClusterLabel", cell_size = 3, show_branch_points = FALSE) + scale_color_manual(values = x.type.color) + theme_classic(base_size = 20))
  dev.off()
  
  png(filename = paste0(temp.outdir, "/monocle.celltypes_by_samples.png"), width = 2100, height = 900)
  print(plot_cell_trajectory(temp.cds, color_by = "ClusterLabel", cell_size = 3, show_branch_points = FALSE) + 
          facet_wrap(~orig.ident, nrow = 2) + scale_color_manual(values = x.type.color) + theme_classic(base_size = 20))
  dev.off()
  
  saveRDS(object = temp.cds, file = paste0(temp.outdir, "/monocle.RDS"))
}
