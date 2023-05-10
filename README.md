## mouse_microvascular_skin_2022-2023
## soup correction in r script 
library(Seurat)
library(sctransform)
library(harmony)
library(SoupX)
library(tidyverse)
library(ggplot2)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_17H/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_17H/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_17H/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_17H.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_39Q/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_39Q/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_39Q/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_39Q.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_41H/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_41H/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_41H/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_41H.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_B3/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_B3/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_B3/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_B3.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_C2/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_C2/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_C2/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_C2.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_D4/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_D4/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_D4/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_D4.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_E1/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_E1/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_E1/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_E1.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_F3/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_F3/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_F3/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_F3.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_I4/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_I4/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_I4/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_I4.data = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_Y3/outs/filtered_feature_bc_matrix'))
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_Y3/outs/raw_feature_bc_matrix'))
sc = SoupChannel(tod, toc)
cluster = read.table('/lustre03/project/6003727/sequencing_datastore/analyses/scRNAseq_skin-senescense_20220928/sample_Y3/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T)
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc)
soup_Y3.data = adjustCounts(sc)

## create seurat object and metadata
soup_skin_samples <- CreateSeuratObject(counts = cbind(soup_17H.data , soup_39Q.data, soup_41H.data, soup_B3.data, soup_C2.data, soup_D4.data, soup_E1.data, soup_F3.data, soup_I4.data, soup_Y3.data), project = "inkmouseskin", min.cells = 5)

soup_skin_samples@meta.data$sample <- c(rep("s39Q", ncol(soup_39Q.data)), rep("sD4", ncol(soup_D4.data)), rep("sC2", ncol(soup_C2.data)), rep("sY3", ncol(soup_Y3.data)), rep("sB3", ncol(soup_B3.data)), rep("sE1", ncol(soup_E1.data)), rep("sI4", ncol(soup_I4.data)), rep("sF3", ncol(soup_F3.data)), rep("s17H", ncol(soup_17H.data)), rep("s41H", ncol(soup_41H.data)))

soup_skin_samples@meta.data$condition <- c(rep("ntr", ncol(soup_39Q.data)), rep("ntr", ncol(soup_D4.data)), rep("ntr", ncol(soup_C2.data)), rep("ntr", ncol(soup_Y3.data)), rep("veh", ncol(soup_B3.data)), rep("veh", ncol(soup_E1.data)), rep("veh", ncol(soup_I4.data)), rep("tr", ncol(soup_F3.data)), rep("tr", ncol(soup_17H.data)), rep("tr", ncol(soup_41H.data)))

soup_skin_samples@meta.data$age <- c(rep("old", ncol(soup_17H.data)), rep("young", ncol(soup_39Q.data)), rep("old", ncol(soup_41H.data)), rep("old", ncol(soup_B3.data)), rep("young", ncol(soup_C2.data)), rep("young", ncol(soup_D4.data)), rep("old", ncol(soup_E1.data)), rep("old", ncol(soup_F3.data)), rep("old", ncol(soup_I4.data)), rep("young", ncol(soup_Y3.data)))

soup_skin_samples@meta.data$batch <- c(rep("batch3", ncol(soup_17H.data)), rep("batch1", ncol(soup_39Q.data)), rep("batch3", ncol(soup_41H.data)), rep("batch2", ncol(soup_B3.data)), rep("batch2", ncol(soup_C2.data)), rep("batch2", ncol(soup_D4.data)), rep("batch3", ncol(soup_E1.data)), rep("batch3", ncol(soup_F3.data)), rep("batch3", ncol(soup_I4.data)), rep("batch2", ncol(soup_Y3.data)))

soup_skin_samples[["percent.mt"]] <- PercentageFeatureSet(soup_skin_samples, pattern = "^mt-")
soup_skin_samples <- PercentageFeatureSet(soup_skin_samples, pattern = "^mt-", col.name = "percent.mt")

saveRDS(soup_skin_samples, file = "soup_skin_samples.rds")

## Quality control & doublet correction & filter & normalization

soup_skin_samples <- readRDS("soup_skin_samples.rds")

VlnPlot(soup_skin_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(soup_skin_samples, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(soup_skin_samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


##Quality control metrics
# Add number of genes per UMI for each cell to metadata

soup_skin_samples$log10GenesPerUMI <- log10(soup_skin_samples$nFeature_RNA) / log10(soup_skin_samples$nCount_RNA)

# Compute percent mito ratio

soup_skin_samples$mitoRatio <- PercentageFeatureSet(object = soup_skin_samples, pattern = "^mt-")
soup_skin_samples$mitoRatio <- soup_skin_samples@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- soup_skin_samples@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)



# Add metadata back to Seurat object
soup_skin_samples@meta.data <- metadata
# Create .RData object to load at any time
save(soup_skin_samples, file="C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/soup_skin_samples.RData")

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=nGene, fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

## filter

soup_skin_samples <- subset(soup_skin_samples, 
                            subset = nGene > 250 & 
                              nGene < 5000 &
                              nUMI >= 500 &
                              nUMI < 40000 & 
                              log10GenesPerUMI > 0.80 &
                              percent.mt < 20)
# and Gene-level filter
counts <- GetAssayData(object = soup_skin_samples, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_soup_skin_samples <- CreateSeuratObject(filtered_counts, meta.data = soup_skin_samples@meta.data)

saveRDS(filtered_soup_skin_samples, file = "filtered_soup_skin_samples.rds")

## perform normalization & UMAP, before finding doublets
# normalizing
filtered_soup_skin_samples <- filtered_soup_skin_samples %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2") %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", reduction.save = "harmony2") %>% 
  RunUMAP(reduction = "harmony2", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony2", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(filtered_soup_skin_samples, reduction = "umap", label = TRUE, label.size = 5)

saveRDS(filtered_soup_skin_samples, file = "fn_soup_skin_samples.rds")

## Doublet correction
# load libraries
BiocManager::install("plger/scDblFinder")
BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

## finding doubletScore 
DefaultAssay(fn_soup_skin_samples) <- "RNA"
sce <- as.SingleCellExperiment(fn_soup_skin_samples)
set.seed(2022)
sce <- scDblFinder(sce, clusters = "seurat_clusters", samples="sample")
sce$scDblFinder.score
sce$scDblFinder.class
fn_soup_skin_samples$scDblFinder.score <- sce$scDblFinder.score
FeaturePlot(fn_soup_skin_samples, features = "scDblFinder.stats")
table(sce$scDblFinder.class)
metadata(sce)$scDblFinder.stats

## add scDblFinder.class column from sce in soup_skin_samples object
fn_soup_skin_samples@meta.data$isDoublet <- sce$scDblFinder.class

## Doublets visualization
DimPlot(fn_soup_skin_samples, reduction = 'umap', group.by = "isDoublet")

## create a new object(fnd_soup_skin_samples) and filter out the doublet
fnd_soup_skin_samples <- subset(fn_soup_skin_samples, isDoublet == "doublet")
DimPlot(fnd_soup_skin_samples, reduction = 'umap', group.by = "isDoublet")
DimPlot(fnd_soup_skin_samples, reduction = 'umap', label = TRUE, label.size = 3)
# normalizing
  fndn_soup_skin_samples <- fnd_soup_skin_samples %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2") %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", reduction.save = "harmony2") %>% 
  RunUMAP(reduction = "harmony2", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony2", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
## save the new object
saveRDS(fndn_soup_skin_samples, file = "C:/Users/~/fnd_soup_skin_samples.rds")

## Cell cycle estimation
## Next, we will use known cell cycle genes to predict in which cell cycle state cells are currently in. We will use this information later on in the clustering analysis to correct for cell cycle state, so that we don't get separate clusters of cell types that correspond to different cycle states.

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

fnd_soup_skin_samples <- CellCycleScoring(fnd_soup_skin_samples, 
                                      s.features = s.genes, 
                                      g2m.features = g2m.genes, 
                                      set.ident = FALSE)

## Let's check how many cells in each Phase we have
fndc_soup_skin_samples <- fnd_soup_skin_samples@meta.data  %>%
  group_by(Phase) %>%
  tally()
ggplot(fndc_soup_skin_samples,aes(Phase,n, fill = Phase)) +
  geom_bar(stat ="identity") +
  scale_fill_brewer(palette = "Dark2")

## load filtered_normalized_doubletCorrected_normalized_soup_skin_samples
fndn_soup_skin_samples <- readRDS("fndn_soup_skin_samples.rds")
DimPlot(fndn_soup_skin_samples, reduction = "umap", label = TRUE, label.size = 5)
## re-perform QC plots to see differences after filter
# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=nGene, fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


## different visualisations
DimPlot(fndn_soup_skin_samples, reduction = "umap", label = TRUE, label.size = 5)
DimPlot(fndn_soup_skin_samples, reduction = "umap", pt.size = .2, split.by = 'batch', label.size = 3)
DimPlot(fndn_soup_skin_samples, reduction = "umap", pt.size = .2, split.by = 'condition', label.size = 3)
DimPlot(fndn_soup_skin_samples, reduction = "umap", pt.size = .1, split.by = 'sample', ncol=5, label=T)
DimPlot(fndn_soup_skin_samples, reduction = "umap", pt.size = .1, split.by = 'age', ncol=3, label=T)

## find markers
fndn_soup_skin_samples.markers <- FindAllMarkers(object = fndn_soup_skin_samples, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

## plot 5 top markers 
markers.to.plot <- head(subset(fndn_soup_skin_samples.markers, cluster==0), n=5)$gene
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==1), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==2), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==3), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==4), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==5), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==6), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==7), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==8), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==9), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==10), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==11), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==12), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==13), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==14), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==15), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==16), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==17), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==18), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==19), n=5)$gene)
markers.to.plot <- append(markers.to.plot, head(subset(fndn_soup_skin_samples.markers, cluster==20), n=5)$gene)

DotPlot(fndn_soup_skin_samples, features = markers.to.plot[!duplicated(markers.to.plot)], cols = c("yellow", "red"), dot.scale = 8) + RotatedAxis() + scale_y_discrete(limits = rev)
ggsave("fndn_soup_skin_samplesDotPlot_top5.png", width=30, height=10, scale=0.7)

## write .csv file of all markers
write.csv(fndn_soup_skin_samples.markers, "C:/Users/pegah/Desktop/skinProject/fromTheTop/fndn_soup_skin_samples.markers.csv")

### The final step would be changing the cluster names (Which are numbers) to the cell type names
### You'll get the name according to you annotation results. 
### Make sure you insert the names in correct order. In this code the name of the cells is according
### to the names that you give to the code. In other word, you always need to name the clusters 
### from cluster 0 to the end. If you skip one cluster, that cluster will lose its name.

## cluster-identification
#to identify the cluster you can use databases such as enrichr (it includes tabula, Pangia, etc. super helpful), Tabula Muris, PangiaoDB, https://www.immunesinglecell.org/atlasList& automated ways 
fndn_soup_skin_samples <- readRDS("fndn_soup_skin_samples.rds")

new.cluster.ids <- c("0-Basal cell of Epidermis", "1-Fibroblasts", "2-Smooth muscle cell","3-Keratinocyte stem cell", "4-Vasculature EC", "5-Epidermal cell", "6-Langerhans cell", "7-Keratinocytes-Epidermal cell", "8-Keratinocyte stem cell", "9-Smoth muscle cell", "10-Epidermal cell-Keratinocyte", "11-basal cell of epidermis", "12-Smooth muscle cell", "13-Keratinocyte stem cell", "14-Fat smooth muscle cell", "15-Epidermal cell", "16-Adiopocyte", "17-Vascular smooth muscle cell", "18-Epidermal cell - Keratinocyte", "19-Melanocyte", "20-Leukocyte")
names(new.cluster.ids) <- levels(fndn_soup_skin_samples)
fndn_soup_skin_samples <- RenameIdents(fndn_soup_skin_samples, new.cluster.ids)
DimPlot(fndn_soup_skin_samples, reduction = "umap", label = TRUE, pt.size = 1.2, label.size = 3)

### Once you figured out the identity of clusters you need to visualize the expression of genes. 
### You can use the following package and codes to make pretty plots.
#to perform high quality pictures
install.packages("viridis")
library(viridis)
cols = c('grey90', inferno(10))
### Using the following codes you will generate high quality images.
tiff('gene name.tiff', units="in", width=12, height=8, res=300, compression = 'lzw')
FeaturePlot(object, features = c("gene_name"), min.cutoff = "q10", max.cutoff = "q90", cols = cols, pt.size = 1.2, order = TRUE)
dev.off()

# Visualize QC metrics as a violin plot
VlnPlot(fndn_soup_skin_samples, group.by = 'seurat_clusters', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#to perform violn plot and umap for specific gene
VlnPlot(soup_skin_samples,  group.by = 'sample', features = "Cdkn2a", assay="RNA")
VlnPlot(soup_skin_samples,  group.by = 'sample', features = "Nav1", assay="RNA")


DimPlot(soup_skin_samples, reduction = "umap", label = TRUE, pt.size = .1, split.by = 'condition', cells.highlight = WhichCells(soup_skin_samples, expression = Cdkn2a > 0))
DimPlot(soup_skin_samples, reduction = "umap", label = TRUE, pt.size = .1, split.by = 'condition', cells.highlight = WhichCells(soup_skin_samples, expression = Cdkn1a > 1))

#feature plot for specific gene

VlnPlot(fndn_soup_skin_samples, features = ("Itgb1")) + NoLegend()
FeaturePlot(fndn_soup_skin_samples, features = ("Itgb1"), order = T)


##to remove some clusters(subclustering)
##to subset target clusters and remove others
#you can subset a cluster to subcluster to
fn_soup_skin_samples_cl1 <- subset(fn_soup_skin_samples, subset = seurat_clusters %in% c('1'))

#to reperfom clustering with subclusters
fn_soup_skin_samples_cl1 <- fn_soup_skin_samples_cl1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(fn_soup_skin_samples_cl1, reduction = "umap", label = TRUE, label.size = 5)
saveRDS(fn_soup_skin_samples_cl1, file = "fn_soup_skin_samples_cl1.rds")

# to find markers of subcluster
fn_soup_skin_samples_cl1 <- FindAllMarkers(object = fn_soup_skin_samples_cl1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(fn_soup_skin_samples_cl1, "C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/fn_soup_skin_samples_cl1.markers.csv")

## finding doublets in subclusters
# load libraries
BiocManager::install("plger/scDblFinder")
BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

DefaultAssay(fn_soup_skin_samples) <- "RNA"
sce <- as.SingleCellExperiment(fn_soup_skin_samples)
set.seed(2022)
sce <- scDblFinder(sce, clusters = "seurat_clusters", samples="sample")
sce$scDblFinder.score
sce$scDblFinder.class
fn_soup_skin_samples$scDblFinder.stats <- sce$scDblFinder.score
FeaturePlot(fn_soup_skin_samples, features = "scDblFinder.stats")
sce$scDblFinder.score %>% hist
table(sce$scDblFinder.class)
metadata(sce)$scDblFinder.stats
saveRDS(fn_soup_skin_samples, file = "fn_soup_skin_samples.rds")

####### _____________________ Different visualisations ___________________________
## Dotplot for spesific genes
cd_genes <- c("Sele", "Cldn5", "Vwf", "Cdh5", "Flt4", "Lyve1", "Prox1", "Pdpn", "Lum", "Dcn", "Vim", "Pdgfra", "Acta2", "Rgs5", "Pdgfrb", "Des", "Lyz2", "Aif1", "Cd68", "Itgax", "Pmel", "Mlana", "Tyrp1", "Dct", "Cald1", "Cnn1", "Cd3d", "Cd3g", "Cd3e", "Lck", "Krt1", "Krt10", "Sbsn", "Krtdap", "Krt5", "Krt14", "Itga6", "Itgb1")
DotPlot(object = fndn_soup_skin_samples, features = cd_genes)

# to identify how many cells are expressing specific gene
df <- data.frame(p16 = (filtered_seurat@assays$RNA@counts["Cdkn2a", ] > 0),
                 sample = filtered_seurat$sample)

p16.pos <- df %>% group_by(sample) %>% count(p16)
p16.pos

df <- data.frame(p21 = (filtered_seurat@assays$RNA@counts["Cdkn1a", ] > 0),
                 sample = filtered_seurat$sample)

p21.pos <- df %>% group_by(sample) %>% count(p21)
p21.pos

# to check genes and getting list of genes in a text file
filtered_seurat@assays[["RNA"]]@data@Dimnames[[1]]["Nav1"]

filtered_seurat@assays[["RNA"]]@data@Dimnames[[1]][991]

filtered_seurat@assays[["RNA"]]@data@Dimnames[[1]][25821]

filtered_seurat@assays[["RNA"]]@data@Dimnames[[1]][25822]

filtered_seurat@assays[["SCT"]]@data@Dimnames[[1]][25821]

filtered_seurat@assays[["SCT"]]@data@Dimnames[[1]][25822]

write.table(filtered_seurat@assays[["SCT"]]@data@Dimnames[[1]],file = "test1.txt")
write.table(filtered_seurat@assays[["RNA"]]@data@Dimnames[[1]],file = "test1.txt")
# to see what cell(based on no. of read ) is expressing a candidate gene
plot(soup_skin_samples@assays$RNA@data["Cdkn2a",], soup_skin_samples$nCount_RNA)
#_______________________________________

# one way of validating doublets in the dataset is to visualize doublet score of subclusters (before filtering out doublets) and find out if the subcluster with high doublet score is expressing nonesense cell-type marker in compare to other subclusters of the same cluster
## subclustering for doublet validation and annotation validation

fn_soup_skin_samples <- readRDS("fn_soup_skin_samples.rds")

fn_soup_skin_samples <- fn_soup_skin_samples %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()
DimPlot(fn_soup_skin_samples, reduction = "umap", label = TRUE, label.size = 5)

## calculating doublet.score in main cluster (before filtering doublets) if not done before
# load libraries
BiocManager::install("plger/scDblFinder")
BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

DefaultAssay(fn_soup_skin_samples) <- "RNA"
sce <- as.SingleCellExperiment(fn_soup_skin_samples)
set.seed(2022)
sce <- scDblFinder(sce, clusters = "seurat_clusters", samples="sample")
sce$scDblFinder.score
sce$scDblFinder.class
fn_soup_skin_samples$scDblFinder.stats <- sce$scDblFinder.score
FeaturePlot(fn_soup_skin_samples, features = "scDblFinder.stats")
sce$scDblFinder.score %>% hist
table(sce$scDblFinder.class)
metadata(sce)$scDblFinder.stats
saveRDS(fn_soup_skin_samples, file = "C:/Users/~/fn_soup_skin_samples_doublet.rds")

## subsetting target clusters and remove others
# you can subset one or more clusters in your data set and remove others 
fn_soup_skin_samples_doublet <- readRDS("fn_soup_skin_samples_doublet.rds")
DimPlot(fn_soup_skin_samples_doublet, reduction = "umap", label = TRUE, label.size = 5)
FeaturePlot(fn_soup_skin_samples_doublet, features = "scDblFinder.stats")
fn_soup_skin_samples_wDoublet_rmKrt <- subset(fn_soup_skin_samples_doublet, subset = seurat_clusters %in% c('1', '2', '5', '6', '7', '8', '12', '13', '14', '15', '16', '17', '18', '19', '20'))

## reperfoming clustering for subclusters
DefaultAssay(fn_soup_skin_samples_wDoublet_rmKrt) <- "SCT"

fn_soup_skin_samples_wDoublet_rmKrt <- fn_soup_skin_samples_wDoublet_rmKrt %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2") %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", reduction.save = "harmony2") %>% 
  RunUMAP(reduction = "harmony2", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony2", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(fn_soup_skin_samples_wDoublet_rmKrt, reduction = "umap", label = TRUE, label.size = 5)
fn_soup_skin_samples_wDoublet_rmKrt.markers <- FindAllMarkers(object = fn_soup_skin_samples_wDoublet_rmKrt, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(fn_soup_skin_samples_wDoublet_rmKrt.markers, "C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/fn_soup_skin_samples_wDoublet_rmKrt.markers.csv")
saveRDS(fn_soup_skin_samples_wDoublet_rmKrt, file = "fn_soup_skin_samples_wDoublet_rmKrt.rds")

fn_soup_skin_samples_wDoublet_rmKrt <- readRDS("fn_soup_skin_samples_wDoublet_rmKrt.rds")

## finding-visiualizing doublets in subclusters
FeaturePlot(fn_soup_skin_samples_wDoublet_rmKrt, features = "scDblFinder.stats", label = TRUE)
FeaturePlot(fn_soup_skin_samples_wDoublet_rmKrt, features = "percent.mt", label = TRUE)

# performing embedding on UMAP with 0.1 resolution to have more generalized clusters to be able to sub-cluster better
fn_soup_skin_samples_wDoublet_rmKrt0.1 <- fn_soup_skin_samples_wDoublet_rmKrt %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()
DimPlot(fn_soup_skin_samples_wDoublet_rmKrt0.1, reduction = "umap", label = TRUE, label.size = 5)
fn_soup_skin_samples_wDoublet_rmKrt0.1.markers0.1 <- FindAllMarkers(object = fn_soup_skin_samples_wDoublet_rmKrt0.1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(fn_soup_skin_samples_wDoublet_rmKrt0.1.markers0.1, "C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/fn_soup_skin_samples_wDoublet_rmKrt0.1.markers0.1.csv")
saveRDS(fn_soup_skin_samples_wDoublet_rmKrt0.1, file = "fn_soup_skin_samples_wDoublet_rmKrt0.1.rds")

FeaturePlot(fn_soup_skin_samples_wDoublet_rmKrt0.1, features = "scDblFinder.stats", label = TRUE, label.size = 5)

# if doublet score inside a cluster is high in a sebcluster and that subcluster or some cells inside it are showing nonsense cell-type markers it means they are doublets and we should correct it

## to rename clusters
new.cluster.ids <- c("0-Muscle cells", "1-vSMC-FB-myoFB", "2-Langerhans cells (dendritic cells)","3-EC", "4-basal cell of epidermis", "5-Epidermal cell-adipocyte", "6-Myocyte", "7-Adipocyte", "8-FB?", "9-SMC-vSMC", "10-Melanocyte", "11-?")
names(new.cluster.ids) <- levels(fn_soup_skin_samples_wDoublet_rmKrt0.1)
fn_soup_skin_samples_wDoublet_rmKrt0.1 <- RenameIdents(fn_soup_skin_samples_wDoublet_rmKrt0.1, new.cluster.ids)
DimPlot(fn_soup_skin_samples_wDoublet_rmKrt0.1, reduction = "umap", label = TRUE, pt.size = 1.2, label.size = 3)

# saving new object
saveRDS(fn_soup_skin_samples_wDoublet_rmKrt0.1, file = "C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/fn_soup_skin_samples_wDoublet_rmKrt0.1.rds")

# subset cluster 0
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0 <- subset(fn_soup_skin_samples_wDoublet_rmKrt0.1, subset = seurat_clusters %in% c('0'))

# performing PCA, UMAP, FindClusters for cl0
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0 <- fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0 %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2") %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", reduction.save = "harmony2") %>% 
  RunUMAP(reduction = "harmony2", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony2", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0, reduction = "umap", label = TRUE, label.size = 5)
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0.markers <- FindAllMarkers(object = fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0.markers, "C:/Users/~/fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0.markers.csv")
saveRDS(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0, file = "fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0.rds")

# doublet score subclusters in cl0
FeaturePlot(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0, features = "scDblFinder.stats", label = TRUE, label.size = 5)

# fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0 <- readRDS("fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0.rds")
## renaming sub-clusters in cl0
new.cluster.ids <- c("0-SMC", "1-Myocyte", "2-Basal cell of epidermis","3-Mayocyte", "4-SMC", "5-SMC", "6-Myocyte", "7-Myocyte", "8-vSMC", "9-Satellite cells", "10-EC-FB", "11-Myocyte-pericyte")
names(new.cluster.ids) <- levels(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0)
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0 <- RenameIdents(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0, new.cluster.ids)
DimPlot(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0, reduction = "umap", label = TRUE, pt.size = 1.2, label.size = 5)
saveRDS(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0, file = "fn_soup_skin_samples_wDoublet_rmKrt0.1_cl0.rds")

# subset cluster 1
fn_soup_skin_samples_wDoublet_rmKrt0.1 <- readRDS("fn_soup_skin_samples_wDoublet_rmKrt0.1.rds")
DimPlot(fn_soup_skin_samples_wDoublet_rmKrt0.1, reduction = "umap", label = TRUE, label.size = 5)

fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1 <- subset(fn_soup_skin_samples_wDoublet_rmKrt0.1, subset = seurat_clusters %in% c('1'))

# perform PCA, UMAP, FindClusters for cl1
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1 <- fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1 %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2") %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", reduction.save = "harmony2") %>% 
  RunUMAP(reduction = "harmony2", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony2", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1, reduction = "umap", label = TRUE, label.size = 5)
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1.markers <- FindAllMarkers(object = fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1.markers, "C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1.markers.csv")
saveRDS(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1, file = "fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1.rds")
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1 <- readRDS("fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1.rds")

# visualizing doublet score subclusters in cl1
FeaturePlot(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1, features = "scDblFinder.stats", label = TRUE, label.size = 5)

## renaming sub-clusters in cl1
new.cluster.ids <- c("0-vSMC", "1-Fibroblasts", "2-Fibroblasts","3-Basal cell of epidermis", "4-vSMC-SMC", "5-vSMC", "6-Myofibroblasts-vSMC", "7-Satellite cells-Myocyte", "8-Myocyte", "9-Schwann Cells", "10-SMC", "11-?", "12-SMC", "13-Keratinocytes")
names(new.cluster.ids) <- levels(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1)
fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1 <- RenameIdents(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1, new.cluster.ids)
DimPlot(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1, reduction = "umap", label = TRUE, pt.size = 1.2, label.size = 5)
saveRDS(fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1, file = "fn_soup_skin_samples_wDoublet_rmKrt0.1_cl1.rds")

## to remove doublets from fn_soup_skin_samples_wDoublet_rmKrt and then validate clusters with key genes
fn_soup_skin_samples_wDoublet_rmKrt <- readRDS("fn_soup_skin_samples_wDoublet_rmKrt.rds")

# add scDblFinder.class column from sce in fn_soup_skin_samples_wDoublet_rmKrt object
DefaultAssay(fn_soup_skin_samples_wDoublet_rmKrt) <- "RNA"

fn_soup_skin_samples_wDoublet_rmKrt@meta.data$isDoublet <- sce$scDblFinder.class

# Doublets visualization
DimPlot(fn_soup_skin_samples, reduction = 'umap', group.by = "isDoublet")

# create a new object(fnd_soup_skin_samples) and delete the doublet
fnd_soup_skin_samples <- subset(fn_soup_skin_samples, isDoublet == "doublet")
DimPlot(fnd_soup_skin_samples, reduction = 'umap', group.by = "isDoublet")
DimPlot(fnd_soup_skin_samples, reduction = 'umap', label = TRUE, label.size = 3)
fnd_soup_skin_samples <- readRDS("fnd_soup_skin_samples.rds")


## DESeq2 differentialy gene expression
# script to perform pseudo-bulk DGA
# setwd("~/Desktop/demo/single_cell_DEG")
library(ExperimentHub)
library(Seurat)
library(DESeq2)
library(tidyverse)
library(pheatmap)

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level
seu.filtered <- readRDS("fndn_soup_skin_samples_rmKrt.0.1.rds ")

View(seu.filtered@meta.data)

# Extract metadata to create new object

metadata <- seu.filtered@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seu.filtered@active.ident)
metadata$samples_id <- paste0(seu.filtered$condition, seu.filtered$sample)

View(metadata)
seu.filtered@meta.data <- metadata

DefaultAssay(seu.filtered)
View(seu.filtered)

cts <- AggregateExpression(seu.filtered, 
                           group.by = c("cluster_id", "samples_id"),
                           assays = 'SCT',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$SCT

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split

splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame

cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

#gsub('.*_(.*)', '\\1', '2-EC_ntrs39Q')

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})




# running DE analysis with EC 
# 1. Get counts matrix
counts_EC <- cts.split.modified$`2-EC`


# 2. generate sample level metadat
colData <- data.frame(samples = colnames(counts_EC), 
                      condition = c('young','young','young','young','tr','tr','tr','veh','veh','veh')) %>%
  column_to_rownames(var = 'samples')

colData

# perform DESeq2 --------

#make sure the row names in colData matches to column names in counts data (here count_EC) 

all(colnames(counts_EC) %in% rownames(colData))
#if it is FALSE run
counts_EC <- counts_EC[, rownames(colData)]

# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_EC,
                              colData = colData,
                              design = ~ condition)


# filter

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

### plots
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "condition")


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = colData[, c("condition", "batch"), drop=F])

# specifying the reference level

dds$condition <- relevel(dds$condition, ref = "young")

# likelihood ratio test ()

dds <- DESeq(dds, test="LRT", reduced=~1)

# export normalized read counts

normCountsDDS <- counts(dds, normalized = T)
write.csv(normCountsDDS, "normCountsDDS.csv")

# Plot dispersion estimates
plotDispEsts(dds)

# results to compare each two condition
res_tr_veh <- results(dds, contrast=c("condition","tr","veh"))
res_veh_young <- results(dds, contrast=c("condition","veh","young"))
res_tr_young <- results(dds, contrast=c("condition","tr","young"))

# output tables
write.table(res_tr_veh, file="deseq2_out_tr_veh", quote=FALSE, sep="\t")
write.table(res_veh_young, file="deseq2_out_veh_young", quote=FALSE, sep="\t")
write.table(res_tr_young, file="deseq2_out_tr_young", quote=FALSE, sep="\t")

## Explore more
# explore results
summary(res_tr_veh)
summary(res_veh_young)
summary(res_tr_young)

# MA plot
plotMA(res_tr_veh)
plotMA(res_veh_young)
plotMA(res_tr_young)

## checking batch effect on deseq 
# Extract metadata to create new object

metadata <- seu.filtered@meta.data

# Setting up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seu.filtered@active.ident)
metadata$samples_id <- paste0(seu.filtered$condition, seu.filtered$sample)
#metadata$batch_id <- factor(seu.filtered$condition, seu.filtered$batch)

View(metadata)
seu.filtered@meta.data <- metadata

DefaultAssay(seu.filtered)
View(seu.filtered)

cts <- AggregateExpression(seu.filtered, 
                           group.by = c("cluster_id", "samples_id"),
                           assays = 'SCT',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$SCT

# transpose

cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split

splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame

cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

#gsub('.*_(.*)', '\\1', '2-EC_ntrs39Q')

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})


# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_EC <- cts.split.modified$`2-EC`


# 2. generate sample level metadat

#to add batch
colData <- data.frame(samples = colnames(counts_EC), 
                      condition = c('young','young','young','young','tr','tr','tr','veh','veh','veh'),
                      batch = c('batch1', 'batch2', 'batch2', 'batch2', 'batch3', 'batch3', 'batch3', 'batch2', 'batch3', 'batch3')) %>%
  column_to_rownames(var = 'samples')
# get more information from metadata



# perform DESeq2 --------

#make sure the row names in colData matches to column names in counts data (here count_EC) 

all(colnames(counts_EC) %in% rownames(colData))
#if it is FALSE run
counts_EC <- counts_EC[, rownames(colData)]

# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_EC,
                              colData = colData,
                              design = ~ condition)


# filter

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
###plots
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","batch")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


## Heatmap of the sample-to-sample distances
# apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
sampleDists <- dist(t(assay(vsd)))

#this heatmap function would calculate a clustering based on the distances between the rows/columns of the distance matrix

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
##Principal component plot of the samples
#This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.
plotPCA(vsd, intgroup=c("condition", "batch"))
plotPCA(vsd, intgroup=c("condition"))
plotPCA(vsd, intgroup=c("batch"))

#naming plots inside the PCA plot using geom_label
library(DESeq2)
## after running example("plotPCA")
z <- plotPCA(vsd)
z + geom_label(aes(label = name))

w <- plotPCA(rld, intgroup=c("condition", "batch"))
w + geom_label(aes(label = name))

w <- plotPCA(rld, intgroup=c("batch"))
w + geom_label(aes(label = name))

w <- plotPCA(rld, intgroup=c("condition"))
w + geom_label(aes(label = name))

## That will only be able to use what's in the 
## or for arbitrary extra things
zz <- plotPCA(rld, returnData = TRUE)
## add some random column 
zz$whatever <- c("FOO","BAR")[sample(1:2, 12, TRUE)]
z +  geom_label(data = zz, aes(label = whatever))

## or if you want it to be cool
library(ggrepel)
z + geom_label_repel(data = zz, aes(label = whatever))



### Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "batch")
DESeq2::plotPCA(rld, ntop = 500, intgroup = "condition")


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# specifying the reference level

dds$condition <- relevel(dds$condition, ref = "young")

# likelihood ratio test ()

dds <- DESeq(dds, test="LRT", reduced=~1)

# Plot dispersion estimates
plotDispEsts(dds)

#results to compare each two condition
res_tr_veh <- results(dds, contrast=c("condition","tr","veh"))
res_veh_young <- results(dds, contrast=c("condition","veh","young"))
res_tr_young <- results(dds, contrast=c("condition","tr","young"))

# output tables
write.table(res_tr_veh, file="deseq2_out_tr_veh", quote=FALSE, sep="\t")
write.table(res_veh_young, file="deseq2_out_veh_young", quote=FALSE, sep="\t")
write.table(res_tr_young, file="deseq2_out_tr_young", quote=FALSE, sep="\t")

##Explore more
#explore results
summary(res_tr_veh)
summary(res_veh_young)
summary(res_tr_young)

#MA plot
plotMA(res_tr_veh)
plotMA(res_veh_young)
plotMA(res_tr_young)

###to subset target sub-clusters from EC cluster
wwwECclean <- subset(wwwEC, subset = seurat_clusters %in% c('0', '1', '2', '3', '5'))

##to reperfom clustering for subclusters
#in case that needs to change assay use below code before clustering (here we have two assay--> RNA and SCT)
DefaultAssay(wwwECclean) <- "SCT"

wwwECclean <- wwwECclean %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2") %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", reduction.save = "harmony2") %>% 
  RunUMAP(reduction = "harmony2", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony2", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(wwwECclean, reduction = "umap", label = TRUE, label.size = 5)
wwwECclean.markers <- FindAllMarkers(object = wwwECclean, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(wwwECclean.markers, "C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/wwwECclean.markers.csv")
saveRDS(wwwECclean, file = "wwwECclean.rds")

##to rename clusters
new.cluster.ids <- c("0-venous EC", "1-vEC", "2-arterial EC","3-lEC", "4-venous EC", "5-EC")
names(new.cluster.ids) <- levels(wwwECclean)
wwwECclean <- RenameIdents(wwwECclean, new.cluster.ids)
DimPlot(wwwECclean, reduction = "umap", label = TRUE, pt.size = 1.2, label.size = 3)

FeaturePlot(wwwECclean, features = "scDblFinder.stats", label = TRUE, label.size = 5)
saveRDS(wwwECclean, file = "wwwECclean.rds")


## to perform DESeq only on clean EC cluster
# counts aggregate to sample level
seu.filtered <- readRDS("wwwECclean.rds")

#write.csv(df_with_special_characters, "first.csv", row.names=FALSE)

View(seu.filtered@meta.data)
# Extract metadata to create new object

metadata <- seu.filtered@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seu.filtered@active.ident)
metadata$samples_id <- paste0(seu.filtered$condition, seu.filtered$sample)

View(metadata)
seu.filtered@meta.data <- metadata

DefaultAssay(seu.filtered)
View(seu.filtered)

cts <- AggregateExpression(seu.filtered, 
                           group.by = c("cluster_id", "samples_id"),
                           assays = 'SCT',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$SCT

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split

splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame

cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

#gsub('.*_(.*)', '\\1', '2-EC_ntrs39Q')

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})




# Let's run DE analysis with EC 
# 1. Get counts matrix


# 2. generate sample level metadat
colData <- data.frame(samples = colnames(cts.split.modified), 
                      condition = c('young','young','young','young','tr','tr','tr','veh','veh','veh'),
                      batch = c('batch1', 'batch2', 'batch2', 'batch2', 'batch3', 'batch3', 'batch3', 'batch2', 'batch3', 'batch3'))

column_to_rownames(var = 'samples')
colData

# get more information from metadata



# perform DESeq2 --------

#make sure the row names in colData matches to column names in counts data (here count_EC) if it is not TRUE whatch the tutorial https://www.youtube.com/watch?v=nks7ibkBud8

all(colnames(counts_EC) %in% rownames(colData))
#if it is FALSE run #counts_EC <- counts_EC[, rownames(colData)]

# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_EC,
                              colData = colData,
                              design = ~ condition)


# filter

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
###plots
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "condition")


# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = colData[, c("condition", "batch"), drop=F])

# specifying the reference level

dds$condition <- relevel(dds$condition, ref = "young")

# likelihood ratio test ()

dds <- DESeq(dds, test="LRT", reduced=~1)

# Plot dispersion estimates
plotDispEsts(dds)

#results to compare each two condition
res_tr_veh <- results(dds, contrast=c("condition","tr","veh"))
res_veh_young <- results(dds, contrast=c("condition","veh","young"))
res_tr_young <- results(dds, contrast=c("condition","tr","young"))

# output tables
write.table(res_tr_veh, file="deseq2_out_tr_veh", quote=FALSE, sep="\t")
write.table(res_veh_young, file="deseq2_out_veh_young", quote=FALSE, sep="\t")
write.table(res_tr_young, file="deseq2_out_tr_young", quote=FALSE, sep="\t")

##Explore more
#explore results
summary(res_tr_veh)
summary(res_veh_young)
summary(res_tr_young)

#MA plot
plotMA(res_tr_veh)
plotMA(res_veh_young)
plotMA(res_tr_young)

##in case you want to make csv file of your counts and import it 
#####to make df for running DESeq for cleanECcluster
# counts aggregate to sample level
CluterECclean <- readRDS("CluterECclean.rds")

View(CluterECclean@meta.data)
# Extract metadata to create new object

metadata <- CluterECclean@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(CluterECclean@active.ident)
metadata$samples_id <- paste0(CluterECclean$condition, CluterECclean$sample)

View(metadata)
CluterECclean@meta.data <- metadata

DefaultAssay(CluterECclean)
View(CluterECclean)

#aggregate counts across the cells (make a count matrix for DESeq2)
cts <- AggregateExpression(CluterECclean, 
                           group.by =  "samples_id",
                           assays = 'SCT',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$SCT
#to export count as a CSV file https://www.biostars.org/p/437713/

write.csv(cts, append = F, 'C:/Users/pegah/Desktop/skinProject/~/ctsCounts.csv', sep = ',', row.names = T, col.names = T, quote = T)



# 2. generate sample level metadata
#I made an excel file name colData


# perform DESeq2 --------
#load data

ctsCounts <- read.csv(file = "ctsCounts.csv" , header = T, row.names=1, sep=",")


colData <- read.csv(file = "colData.csv" , header = T, row.names=1, sep=",")

#to see that row names in colData matches to column names in counts data (here count_EC) if it is not TRUE whatch the tutorial https://www.youtube.com/watch?v=nks7ibkBud8

all(colnames(ctsCounts) %in% rownames(colData))
#if it is FALSE run #counts_EC <- counts_EC[, rownames(colData)]

# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = ctsCounts,
                              colData = colData,
                              design = ~ condition)


# filter

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
###plots
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)

select <- order(rowMeans(counts(dds,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","batch")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


# Plot PCA
DESeq2::plotPCA(rld, ntop = 500, intgroup = "condition")

## Plot PCA 
z <- plotPCA(vsd)
z + geom_label(aes(label = name))

#w <- plotPCA(rld, intgroup=c("condition", "batch"))
#w + geom_label(aes(label = name))

w <- plotPCA(rld, intgroup=c("batch"))
w + geom_label(aes(label = name))

w <- plotPCA(rld, intgroup=c("condition"))
w + geom_label(aes(label = name))

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = colData[, c("condition", "batch"), drop=F])

# specifying the reference level

dds$condition <- relevel(dds$condition, ref = "young")

# likelihood ratio test ()

dds <- DESeq(dds, test="LRT", reduced=~1)

# Plot dispersion estimates
plotDispEsts(dds)

#results to compare each two condition
ECclean_res_tr_veh <- results(dds, contrast=c("condition","tr","veh"))
ECclean_res_veh_young <- results(dds, contrast=c("condition","veh","young"))
ECclean_res_tr_young <- results(dds, contrast=c("condition","tr","young"))

# output tables
write.table(ECclean_res_tr_veh, file="ECclean_out_tr_veh", quote=FALSE, sep="\t")
write.table(ECclean_res_veh_young, file="ECclean_out_veh_young", quote=FALSE, sep="\t")
write.table(ECclean_res_tr_young, file="ECclean_out_tr_young", quote=FALSE, sep="\t")

##Explore more
#explore results
summary(ECclean_res_tr_veh)
summary(ECclean_res_veh_young)
summary(ECclean_res_tr_young)

#MA plot
plotMA(ECclean_res_tr_veh)
plotMA(ECclean_res_veh_young)
plotMA(ECclean_res_tr_young)

#to manipulate dataset and perform ploting
#load data(norm count )
normCountsDDS <- read.csv("normCountsDDS.csv", header = T, row.names = 1)

ECclean_DESeq2_out_tr_young$sig <- ifelse(ECclean_DESeq2_out_tr_young$pvalue <= 0.05, "yes", "no")

ECclean_DESeq2_out_tr_young <- na.omit(ECclean_DESeq2_out_tr_young)
library(ggplot2)

#plotMA
ggplot2::ggplot(ECclean_DESeq2_out_tr_young, aes(x = log10(baseMean), y = log2FoldChange, color = sig)) +
  geom_point()


#res<-res[order(res$padj),]

###Plot counts### http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
#d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
#returnData=TRUE)
#library("ggplot2")
#ggplot(d, aes(x=condition, y=count)) + 
#  geom_point(position=position_jitter(w=0.1,h=0)) + 
#  scale_y_log10(breaks=c(25,100,400))



### pathway analysis

#pathway analysis
library(Seurat)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(data.table)

##load data
DEG_tr_veh <- read.csv(file = "DEG_tr_veh.csv" , header = T, sep=",")
DEG_tr_young <- read.csv(file = "DEG_tr_young.csv" , header = T, sep=",")
DEG_veh_young <- read.csv(file = "DEG_veh_young.csv" , header = T, sep=",")


##Enrichment##
#get gene database
#all gene
all_gene_sets = msigdbr(species = "Mus musculus")
head(all_gene_sets)

#halmark
h_gene_sets = msigdbr(species = "mouse", category = "H")
head(h_gene_sets)

#to check if it is a data.frame
class(h_gene_sets)

#to define significant genes
# to look at padj

ggplot(DEG_tr_veh, aes(x=pvalue)) +
  geom_histogram(bins = 10) +
  theme_classic() 

#to see if we have padj==0
table(DEG_tr_veh$pvalue == 0)


#cutoff non-significant padj bc padj should not be zero (arbitrary step)
signif <- ECclean_DESeq2_out_tr_veh %>%
  filter(padj <= 1E-16)

table(ECclean_DESeq2_out_tr_veh$padj == 0)

###to get matching columns between my dataset and databese

gene.name <- unique(DEG_tr_veh$gene_id)
H.gene.symbol <- select(h_gene_sets, gs_name, gene_symbol)

#run enrichment
enrich.H<- enricher(gene = gene.name, TERM2GENE = H.gene.symbol)

##Extract results
class(enrich.H)

head(enrich.H@result)

class(enrich.H@result$GeneRatio)

##format the results to be able to plot
#separate ratios into 2 columns of data (to see how big the overlap between the gene set and halmark database)

enrich.H.df <- enrich.H@result %>% 
    separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
#convert to numeric
         mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 

#Calculate k/K to see how many significant genes are in the gene set relative to haw big the geneset is
            mutate("k.K"=size.overlap.term/size.term) 



##visualize results##
enrich.H.df %>% 
  filter(p.adjust <= 0.05) %>%
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Mtb significant genes \nenriched in Hallmark gene sets")


################################
######GSEA########
#this method could not be used for more than two condition comparison, if your data has more than two groups or condition separate group of interest before

#load data
DEG_ECclean_DESeq2_trVSveh <- read.csv(file = "DEG_ECclean_DESeq2_trVSveh.csv" , header = T, sep=",", row.names = 1)


#Get gene set database
H = msigdbr(species = "mouse", category = "H")
head(h_gene_sets)

class(H)
#format gene set database
head(H)

H.gene.name.ls <- H %>% 
  select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()

# Calculate fold change if your data is not coming from DESeq analysis it does not have LFC and following commonds help to cakculate it


H.gene.name2 <- read.csv(file = "H.gene.name2.csv" , header = T, sep=",")

#Extract expression data
FC <- as.data.frame(DEG$log2FoldChange) %>% 
  as.data.frame(DEG$gene_id)%>% 
  #Move gene IDs from rownames to a column
  rownames_to_column("gene_id")

#format for gsea 
FC.vec <- FC$"DEG$log2FoldChange"
names(FC.vec) <- FC$gene_id

#set score type 
min(FC.vec)
max(FC.vec)

#bc we have + & - LFC standardiz it 

scoreType <- "std"

#Run GSEA#####
#the number of nperm=1000 can increase till 1milion depending on the pvalue that you have
gsea.H <- fgseaSimple(pathways = H.gene.name.ls,
                      stats = FC.vec,
                      scoreType = scoreType,
                      nperm=1000)

#### plot gsea results ####
class(gsea.H)

gsea.H %>% 
  filter(pval <= 0.05) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(pathway = gsub("HALLMARK_","", pathway),
         pathway = gsub("_"," ", pathway)) %>% 
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  #Force equal max min
  lims(y=c(-3.2,3.2)) +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES), pval <= 0.05",
       x="Gene set",
       title = "Hallmark & senescence GSEA \nDown in +Mtb <--         --> Up in +Mtb")

###


gsea.H %>% 
   #Beautify descriptions by removing _ and HALLMARK
  mutate(pathway = gsub("HALLMARK_","", pathway),
         pathway = gsub("_"," ", pathway)) %>% 
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  #Force equal max min
  lims(y=c(-3.2,3.2)) +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES)",
       x="Gene set",
       title = "Hallmark & senescence GSEA \nDown in +Mtb <--         --> Up in +Mtb")

##output##
#to write a table
gsea.H <- apply(gsea.H,2,as.character)
write.table(gsea.H, file="GSEA2_only_tr_veh_EC.csv", sep=",", row.names = F)

# Save a single object to a file
saveRDS(gsea.H, "GSEA2_only_tr_veh_EC.rds")

## convert gene ID to Symbol
edox <- setReadable(enrich.H, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)









