library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# find markers for every cluster compared to all remaining cells, report only the positive ones
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.75)
cluster0_1.markers <- FindMarkers(pbmc, ident.1 = 0,ident.2 = 1, min.pct = 0.75)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
    "CD8A"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

########## run monocle on 
library(monocle)

sample.data <- as.matrix(pbmc@assays$RNA@counts[,c(pbmc@meta.data$seurat_clusters==2|pbmc@meta.data$seurat_clusters==5)])

sample_sheet <- data.frame(colnames(sample.data))
rownames(sample_sheet) <- colnames(sample.data)
gene_annotation <- data.frame(rownames(sample.data))
colnames(gene_annotation) <- c("gene_short_name")
rownames(gene_annotation) <- rownames(sample.data)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
sample.data <- newCellDataSet(as.matrix(sample.data),
                       phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())

sample.data <- estimateSizeFactors(sample.data)
sample.data <- estimateDispersions(sample.data)

sample.data <- detectGenes(sample.data, min_expr = 0.1)

gene1_id <- row.names(subset(fData(sample.data), gene_short_name == "CD14"))
gene2_id <- row.names(subset(fData(sample.data),
                             gene_short_name == "FCGR3A"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "CD14+", classify_func =
                     function(x) { x[gene1_id,] >= 1 & x[gene2_id,] <= 1 })
cth <- addCellType(cth, "FCGR3A+", classify_func = function(x)
{ x[gene1_id,] < 1 & x[gene2_id,] > 1 })
sample.data <- classifyCells(sample.data, cth, 0.1)

sample.data <- reduceDimension(sample.data, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
sample.data <- clusterCells(sample.data, num_clusters = 2)
plot_cell_clusters(sample.data, 1, 2, color = "CellType",
                   markers = c("CD14", "FCGR3A"))

#Ordering for clustering
diff_test_res <- differentialGeneTest(sample.data,
                                      fullModelFormulaStr = "~CellType")

ordering_genes <- row.names (subset(diff_test_res, qval <0.01))
sample.data <- setOrderingFilter(sample.data, ordering_genes)

sample.data <- reduceDimension(sample.data, max_components = 2,
                            method = 'DDRTree')
sample.data <- orderCells(sample.data)
plot_cell_trajectory(sample.data, color_by = "State")
plot_cell_trajectory(sample.data, color_by = "Pseudotime")

blast_genes <- row.names(subset(fData(sample.data),
                                gene_short_name %in% c("CD14", "FCGR3A")))
plot_genes_jitter(sample.data[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)
                  
my_genes <- row.names(subset(fData(sample.data),
          gene_short_name %in% c("CD14", "FCGR3A","LYZ")))
cds_subset <- sample.data[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "CellType")
