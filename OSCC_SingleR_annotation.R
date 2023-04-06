library(Seurat)
library(SingleR)
library(scater)

##### Read Seurat object
out.dir <- file.path("/seurat_output_dir")
x.object <- readRDS(paste0(out.dir, "/Seurat.RDS"))

x.count <- x.object@assays$RNA@data
hpca.se <- HumanPrimaryCellAtlasData()

common <- intersect(rownames(x.count), rownames(hpca.se))

hpca.se <- hpca.se[common,]

x.count <- x.count[common,]

pred.hpca <- SingleR(test = x.count, ref = hpca.se, labels = hpca.se$label.main)

x.object@meta.data$labels <- pred.hpca$labels

cell.type = data.frame()

for(cluster in unique(x.object@meta.data$seurat_clusters)){
  singler.ex <- subset(x.object@meta.data, seurat_clusters == cluster)
  singler.data <- data.frame(table(singler.ex$labels))
  singler.data$Freq <- 100*singler.data$Freq/sum(singler.data$Freq)
  singler.data <- singler.data[order(singler.data$Freq, decreasing = T),]
  cell.type <- rbind(cell.type, data.frame(clusters = cluster, cell.types = singler.data$Var1[1:3], ratio = singler.data$Freq[1:3]))
}

cell.type$clusters <- as.numeric(as.character(cell.type$clusters))

write.table(x = cell.type, file = paste0(out.dir, "/singleR.txt"), quote = F, sep = "\t", row.names = F)
