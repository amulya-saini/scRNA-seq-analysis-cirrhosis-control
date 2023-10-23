# set the working directory
setwd("data")

# reading the working directory into a variable
wd <- getwd()

wd
# loading required packages
library(Seurat)
library(Matrix)
library(patchwork)
library(SingleR)
library(SingleCellExperiment)
library(harmony)
library(ggplot2)

#getting the folder names in the directory into a variable
loc <- as.list(list.files(path = "data", 
                          recursive = FALSE, full.names = FALSE))

# Creating an empty list to store Seurat objects
seurat_obj_list <- list()

# writing a for loop to iterate over the folders
for (i in loc) {
  
  # Getting the path to the current folder
  folder_path <- paste(wd, i, sep = "/")
  
  # reading the matrix file as counts
  counts <- readMM(paste(folder_path, "matrix.mtx", sep = "/"))
  
  # reading the barcodes.tsv file as barcodes
  barcodes <- read.table(paste(folder_path, "barcodes.tsv", sep = "/"), header = FALSE, 
                         col.names = "barcodes")
  
  #attaching the barcodes to counts as column names
  colnames(counts) <- barcodes$barcodes
  
  # Reading the genes.tsv file as gene
  genes <- read.delim(paste(folder_path, "genes.tsv", sep = "/"), header = FALSE)
  
  # attaching the genes files as rownames to the counts matrix
  rownames(counts) <- genes$V2
  
  # Creating a Seurat object and filtering out cells with less 200 genes and 
  # genes that are expressed in less than 3 cells
  seurat_obj <- CreateSeuratObject(counts = counts,
                                   assay = "RNA",
                                   project = i,
                                   min.cells = 3,
                                   min.features = 200
  )
  
  # Add the Seurat object to the list
  seurat_obj_list[[i]] <- seurat_obj
}

# creating individual variables to merge the datasets
cirr_1 <- seurat_obj_list$cirr_1
cirr_a <- seurat_obj_list$cirr_a
cirr_b <- seurat_obj_list$cirr_b

ctrl_1 <- seurat_obj_list$ctrl_1
ctrl_a <- seurat_obj_list$ctrl_a
ctrl_b <- seurat_obj_list$ctrl_b

# adding metadata "sample"
cirr_1$sample <- "p"
cirr_a$sample <- "a"
cirr_b$sample <- "b"

ctrl_1$sample <- "p"
ctrl_a$sample <- "a"
ctrl_b$sample <- "b"

#merging the data
# cirrhotic
cirr_data <- merge(cirr_1,c(cirr_a, cirr_b), 
                 add.cell.ids = c("cirr_1","cirr_a","cirr_b"))

# control
ctrl_data <- merge(ctrl_1,c(ctrl_a, ctrl_b), 
                   add.cell.ids = c("ctrl_1","ctrl_a","ctrl_b"))

# quality control cirrhotic
cirr_data <- PercentageFeatureSet(cirr_data, "^MT-", col.name = "percent.mt")
cirr_data <- PercentageFeatureSet(cirr_data, "^RP[SL]",col.name = "percent.ribo")

# quality control control
ctrl_data <- PercentageFeatureSet(ctrl_data, "^MT-", col.name = "percent.mt")
ctrl_data <- PercentageFeatureSet(ctrl_data, "^RP[SL]",col.name = "percent.ribo")

# Visualizing QC metrics as a violin plot
feat <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
VlnPlot(cirr_data, group.by = "orig.ident", features = feat, pt.size = 0.1, ncol = 3) +NoLegend()
VlnPlot(ctrl_data, group.by = "orig.ident", features = feat, pt.size = 0.1, ncol = 3) +NoLegend()

# visualization as scatter plot
FeatureScatter(cirr_data, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
FeatureScatter(ctrl_data, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

# filtering out cells with high mitochondrial genes
cirr_mito <- WhichCells(cirr_data, expression = percent.mt < 5)
ctrl_mito <- WhichCells(ctrl_data, expression = percent.mt < 5)

# subsetting the data after quality control
cirr.filt <- subset(cirr_data, cells = cirr_mito)
ctrl.filt <- subset(ctrl_data, cells = ctrl_mito)

#visualization after removing mitochondrial genes
VlnPlot(cirr.filt, group.by = "orig.ident", features = feat, pt.size = 0.1, ncol = 3) +NoLegend()
VlnPlot(ctrl.filt, group.by = "orig.ident", features = feat, pt.size = 0.1, ncol = 3) +NoLegend()

# normalizing the data
cirr.filt <- NormalizeData(cirr.filt)
ctrl.filt <- NormalizeData(ctrl.filt)

# cell cycle scoring
cirr.filt <- CellCycleScoring(object = cirr.filt, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)
ctrl.filt <- CellCycleScoring(object = ctrl.filt, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)
# visualization
VlnPlot(cirr.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
         ncol = 4, pt.size = 0.1)
VlnPlot(ctrl.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
         ncol = 4, pt.size = 0.1)

# Feature selection
cirr.filt = FindVariableFeatures(cirr.filt, verbose = F)
ctrl.filt = FindVariableFeatures(ctrl.filt, verbose = F)

# visualization the top 10 features
# cirrhotic
top10 <- head(VariableFeatures(cirr.filt), 10)
plot1 <- VariableFeaturePlot(cirr.filt)
x <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


# #control
top10_a <- head(VariableFeatures(ctrl.filt), 10)
plot2 <- VariableFeaturePlot(ctrl.filt)
y <-LabelPoints(plot = plot2, points = top10_a, repel = TRUE)



# sclaing the data
cirr.filt = ScaleData(cirr.filt, vars.to.regress = c("nFeature_RNA", "percent.mt"),
                      verbose = F)
ctrl.filt = ScaleData(ctrl.filt, vars.to.regress = c("nFeature_RNA", "percent.mt"),
                      verbose = F)

# Principal Component Analysis
cirr.filt = RunPCA(cirr.filt, verbose = F)
ctrl.filt = RunPCA(ctrl.filt, verbose = F)

# # Elbow plot 
a <- ElbowPlot(cirr.filt)
b <- ElbowPlot(ctrl.filt)
a+b

# Finding Neighbors and clusters
cirr.filt <- FindNeighbors(cirr.filt, dims = 1:15)
cirr.filt <- FindClusters(cirr.filt, resolution = c(0.5))

ctrl.filt <- FindNeighbors(ctrl.filt, dims = 1:15)
ctrl.filt <- FindClusters(ctrl.filt, resolution = c(0.5))

# Umap cluster visualization (dimensionality reduction)
cirr.filt = RunUMAP(cirr.filt, dims = 1:15, verbose = F, umap.method = "uwot")
DimPlot(cirr.filt, reduction = "umap")

# chceking for batch effects in cirrhotic samples
DimPlot(cirr.filt, reduction = "umap", group.by = "sample")

ctrl.filt = RunUMAP(ctrl.filt, dims = 1:15, verbose = F, umap.method = "uwot")
DimPlot(ctrl.filt, reduction = "umap")

# checking for batch effects in control samples
DimPlot(ctrl.filt, reduction = "umap", group.by = "sample")

# integrating the data
# using Harmony
cirr.int <- cirr.filt %>% RunHarmony(group.by.vars = "sample", plot_convergence = FALSE)

# calculating the neighbors, clusters and UMAP again
cirr.int <- cirr.int %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualization
DimPlot(cirr.int, reduction = "umap", group.by = "sample")

#control
ctrl.int <- ctrl.filt %>% RunHarmony(group.by.vars = "sample", plot_convergence = FALSE)

# calculating the neighbors, clusters and UMAP again
ctrl.int <- ctrl.int %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualization
#DimPlot(ctrl.int, reduction = "umap", group.by = "sample")

#cluster visualization
DimPlot(cirr.int, reduction = "umap", label = TRUE)
DimPlot(ctrl.int, reduction = "umap", label = TRUE)

#finding markers
#cirrhotic
for (i in unique(Idents(cirr.int))) {
  
  # Run the Find Markers function for the current cluster
  markers <- FindMarkers(cirr.int, ident.1 = i, grouping.var = "type", 
                         min.pct = 0.25, logfc.threshold = log(1))
  
  # Add a new column with the gene names and reset the rownames
  markers <- cbind(gene = rownames(markers), markers)
  rownames(markers) <- NULL
  
  # Assign the markers dataframe to a new variable with a dynamic name based on the cluster identity
  assign(paste0("cirr.marker", i), markers)
  
}
cirr.int <- RenameIdents(cirr.int, "0" = "T cells")
cirr.int <- RenameIdents(cirr.int, "1" = "Natural Killer cells")
cirr.int <- RenameIdents(cirr.int, "2" = "E cells")
cirr.int <- RenameIdents(cirr.int, "3" = "GD T cells")
cirr.int <- RenameIdents(cirr.int, "4" = "Dendritic cells")
cirr.int <- RenameIdents(cirr.int, "5" = "Hepatocytes")
cirr.int <- RenameIdents(cirr.int, "6" = "E cells")
cirr.int <- RenameIdents(cirr.int, "7" = "Dendritic cells")
cirr.int <- RenameIdents(cirr.int, "8" = "Plasma cells")
cirr.int <- RenameIdents(cirr.int, "9" = "Fibroblasts")
cirr.int <- RenameIdents(cirr.int, "10" = "B cells")
cirr.int <- RenameIdents(cirr.int, "11" = "LyECs")
cirr.int <- RenameIdents(cirr.int, "12" = "Mast cells")
cirr.int <- RenameIdents(cirr.int, "13" = "E cells")

# plotting the umap with cell annotations
e <- DimPlot(cirr.int, reduction = "umap", label = TRUE)

#Control
for (i in unique(Idents(ctrl.int))) {
  
  # Run the Find Markers function for the current cluster
  markers <- FindMarkers(ctrl.int, ident.1 = i, grouping.var = "type", 
                         min.pct = 0.25, logfc.threshold = log(1))
  
  # Add a new column with the gene names and reset the rownames
  markers <- cbind(gene = rownames(markers), markers)
  rownames(markers) <- NULL
  
  # Assign the markers dataframe to a new variable with a dynamic name based on the cluster identity
  assign(paste0("ctrl.marker", i), markers)
  
}

# renaming the clusters
ctrl.int <- RenameIdents(ctrl.int, "0" = "T cells")
ctrl.int <- RenameIdents(ctrl.int, "1" = "Dendritic cells")
ctrl.int <- RenameIdents(ctrl.int, "2" = "E cells")
ctrl.int <- RenameIdents(ctrl.int, "3" = "NK cells")
ctrl.int <- RenameIdents(ctrl.int, "4" = "E cells")
ctrl.int <- RenameIdents(ctrl.int, "5" = "E cells")
ctrl.int <- RenameIdents(ctrl.int, "6" = "B cells")
ctrl.int <- RenameIdents(ctrl.int, "7" = "GD T cells")
ctrl.int <- RenameIdents(ctrl.int, "8" = "Dendritic cells")
ctrl.int <- RenameIdents(ctrl.int, "9" = "unknown cells(RP^)")
ctrl.int <- RenameIdents(ctrl.int, "10" = "Dendritic cells")
ctrl.int <- RenameIdents(ctrl.int, "11" = "Fibroblasts")
ctrl.int <- RenameIdents(ctrl.int, "12" = "E cells")
ctrl.int <- RenameIdents(ctrl.int, "13" = "Plasmacytoid dendritic cells")

# plotting the umap with cell annotations
f <- DimPlot(ctrl.int, reduction = "umap", label = TRUE)
e + f

# cell type annotation using SingleR
# Cirrhosis
ref <- celldex::HumanPrimaryCellAtlasData()
results_main <- SingleR(test = as.SingleCellExperiment(cirr.int), ref = ref, labels = ref$label.main)
cirr.int$main.labels <- results_main$labels
g <- DimPlot(cirr.int, reduction = "umap", group.by = "main.labels", label = TRUE)

# Control
results_main <- SingleR(test = as.SingleCellExperiment(ctrl.int), ref = ref, labels = ref$label.main)
ctrl.int$main.labels <- results_main$labels
h <- DimPlot(ctrl.int, reduction = "umap", group.by = "main.labels", label = TRUE)

# searching marker
search_value <- "PDPN"
search_value <- "COL1A1"
search_value <- "FSP1"
search_value <- "GFAP"
search_value <- "MMP2"

# getting a list of all marker dataframes
df_nums <- 0:13
df_names <- paste0("cirr.marker", df_nums)
df_names
df_list <- mget(df_names)

# cirrhotic
# loop over each dataframe in the list to check if a marker is present in the cluster

for (i in df_list) {
  # check if the search value is present in the "gene" column of the current dataframe
  if (search_value %in% i$gene) {
    # print a message if the value is found
    cat("Found", search_value, "in", deparse(substitute(i)), "\n")
  } else {
    # print a message if the value is not found
    cat("Did not find", search_value, "in", deparse(substitute(i)), "\n")
  }
}

# Control
# loop over each dataframe in the list to check if a marker is present in the cluster
df1_names <- paste0("ctrl.marker", df_nums)
df1_names
df1_list <- mget(df1_names)

for (df in df1_list) {
  # check if the search value is present in the "gene" column of the current dataframe
  if (search_value %in% df$gene) {
    # print a message if the value is found
    cat("Found", search_value, "in", deparse(substitute(df)), "\n")
  } else {
    # print a message if the value is not found
    cat("Did not find", search_value, "in", deparse(substitute(df)), "\n")
  }
}
which(cirr.marker11[, 1] == "PDPN")[1]
which(cirr.marker9[, 1] == "COL1A1")[1]
which(cirr.marker6[, 1] == "MMP2")[1]
which(ctrl.marker2[, 1] == "CD34")[1]
x <-VlnPlot(cirr.int, features = cirr.marker6[235,1]) + theme(legend.position = 'none')
y <-VlnPlot(ctrl.int, features = ctrl.marker2[29,1]) + theme(legend.position = 'none')
x+y
x
y










