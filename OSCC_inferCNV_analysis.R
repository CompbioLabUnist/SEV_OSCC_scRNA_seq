# Run InferCNV for each patient
##### Make an input file
{
  library(Seurat)
  library(infercnv)
  library(plyr)
  library(data.table)
  library(doMC)
  #####   Set parameter
  x.dir <- file.path("/Individual_patient_seurat_obj_dir")
  x.outdir <- file.path("/Individual_patient_InferCNV_output_dir")
  x.patients <- paste0("P", c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10"))
  registerDoMC(11)
  foreach (i = 1:10) %dopar%  {
    x.patient = x.patients[i]
    print(x.patient)
    patient.outdir <- file.path(x.outdir, x.patient)
    if (!dir.exists(patient.outdir)) { dir.create(patient.outdir) }
    ##### Read seurat object
    x.file <- paste0(x.dir, "/", x.patient, "/", x.patient, ".Seurat.celltype.RDS")
    x.obj.filter <- readRDS(x.file)
    #####   Make input files
    x.case.celltyeps <- c("Epithelial Cell")
    x.ctl.celltypes <- c("Endothelial Cell")
    temp.outdir <- patient.outdir
    ##### Extract case and control cells
    temp.subset <- subset(x.obj.filter, subset = ClusterLabel %in% c(x.case.celltyeps, x.ctl.celltypes))
    temp.subset@meta.data$celltype <- paste0("C", temp.subset@meta.data$seurat_clusters, "-", temp.subset@meta.data$ClusterLabel)
    temp.count <- temp.subset@assays$RNA@counts #Raw UMI count
    temp.celltype.df <- data.frame(id = rownames(temp.subset@meta.data), celltype = temp.subset@meta.data$celltype)
    temp.count.df <- data.frame(table(temp.celltype.df$celltype))
    temp.count.df <- temp.count.df[ temp.count.df$Freq >= 10,]
    temp.celltype.df <- temp.celltype.df[temp.celltype.df$celltype %in% temp.count.df$Var1,]
    temp.count <- temp.count[, colnames(temp.count) %in% temp.celltype.df$id]
    ##### Convert to a dataframe
    temp.count.df <- as.data.frame(temp.count)
    ##### Save as a text file
    fwrite(x = temp.count.df, file = paste0(temp.outdir, "/", x.patient, ".Seurat.case_epi_control_endo.matrix"), quote = FALSE, sep = "\t", row.names = T)
    fwrite(x = temp.celltype.df, file = paste0(temp.outdir, "/", x.patient, ".Seurat.case_epi_control_endo.annotation"), quote = F, sep = "\t", row.names = F, col.names = F)
  }
}

##### Run inferCNV
{
  library(infercnv)
  library(doMC)
  registerDoMC(11)
  dirs = Sys.glob(file.path("/Individual_patient_InferCNV_output_dir/P*"))
  thread = 15
  foreach (i = 1:10) %dopar%  {
    # for(i in 1:length(dirs)){
    x.dir = dirs[i]
    x.id <- basename(x.dir)
    x.gene.file <- "reference_dir/genes.gtf.pos.rm_duplication"
    ##### Read expression matrix and annotation file
    x.matrix.file <- paste0(x.dir, "/", x.id, ".Seurat.case_epi_control_endo.matrix")
    x.ann.file <- paste0(x.dir, "/", x.id, ".Seurat.case_epi_control_endo.annotation")
    x.outdir <- paste0(x.dir, "/cluster_by_CNV_using_case_epi_control_endo_wo_figure")
    if (!dir.exists(x.outdir)) { dir.create(x.outdir) }
    setwd(x.outdir)
    ##### Run inferCNV
    x.ann.df <- read.table(file = x.ann.file, header = F, sep = "\t")
    x.ref.group <- unique(x.ann.df$V2[grep(pattern= "Endothelial", x= x.ann.df$V2)])
    ##### Create InferCNV object
    x.obj <- CreateInfercnvObject(raw_counts_matrix = x.matrix.file, gene_order_file = x.gene.file, annotations_file = x.ann.file,
                                  ref_group_names = x.ref.group,
                                  delim = "\t", chr_exclude = c("X", "Y"))
    print(x.outdir)
    x.infercnv.obj <- infercnv::run(infercnv_obj = x.obj, cutoff = 0.1, cluster_by_groups = T, denoise = T, HMM = F, num_threads = thread, out_dir = "./", no_plot=T, analysis_mode='subclusters')
  }
}
