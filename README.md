# seuratHelpR
#### A group of R functions to help manipulate and visualize scRNA-seq data (seurat objects) in R 

### description/origin
These functions are built for seurat versions >= 5. Many functions in this package are influenced by or originate entirely from Shekhar, K. et al. (https://doi.org/10.1016/j.cell.2016.07.054) and Loo L. et al. (https://doi.org/10.1038/s41467-018-08079-9)

### installation
```
library(devtools)

install_github('https://github.com/jkniehaus/seuratHelpR.git')
library(seuratHelpR)
```

Each R function is listed within the 'man' directory. Information regarding utility and parameters can be found using `?function` in R.
e.g. `?binomcount.test`

# Functions Overview

## 1. `binomcount.test`
**Description:** Performs a binomial test to determine the probability of differences in the proportion of non-zero values between groups in a Seurat object.
- **Parameters:** 
  - `object`: Seurat object.
  - `cells.1`, `cells.2`: Groups for comparison.
  - `effect.size`: Minimum effect size to filter results.
- **Examples:** `binomcount.test(object, cells.1, cells.2, effect.size)`

---

## 2. `markers.binom`
**Description:** Identifies differentially expressed genes (DEGs) between groups using a binomial test.
- **Parameters:**
  - `object`: Seurat object.
  - `clust.1`, `clust.2`: Clusters or groups for comparison.
  - `effect.size`, `assay`: Effect size threshold and Seurat assay.
  - `posFrac`: Minimum proportion of expressing cells to include genes.
- **Examples:** `markers.binom(object, clust.1, clust.2, effect.size, assay)`

---

## 3. `ClusterCentroids`
**Description:** Computes cluster centroids in principal component (PC) space.
- **Parameters:**
  - `object`: Seurat object.
  - `reduction.use`: Dimensionality reduction technique (default: "pca").
  - `pcs.use`: Principal components to use.
- **Examples:** `ClusterCentroids(object, reduction.use = "pca", pcs.use = 1:10)`

---

## 4. `ComputeClusterDistances`
**Description:** Calculates distances between clusters based on centroids or nearest neighbors.
- **Parameters:**
  - `object`: Seurat object.
  - `reduction.use`: Dimensionality reduction technique.
  - `dist.type`: Distance calculation type ('centroid' or 'nn').
  - `pcs.use`: Principal components to use.
- **Examples:** `ComputeClusterDistances(object, reduction.use = "pca", dist.type = "centroid")`

---

## 5. `merge.clusters.DE`
**Description:** Merges clusters without enough DEGs to distinguish them.
- **Parameters:**
  - `object`: Seurat object.
  - `min.de.genes`: Minimum DEGs to retain clusters.
  - `effect.size`, `pval.cutoff`: Thresholds for DEG calling.
- **Examples:** `merge.clusters.DE(object, min.de.genes, effect.size, pval.cutoff)`

---

## 6. `findDoublets`
**Description:** Identifies and removes doublets using the DoubletFinder package.
- **Parameters:**
  - `object`: Seurat object.
  - `cores`, `pcs`: Number of cores and principal components.
- **Examples:** `findDoublets(object, cores = 1, pcs = 30)`

---

## 7. `modify_vlnplot`
**Description:** Customizes violin plots generated from Seurat VlnPlot.
- **Parameters:**
  - `object`: Seurat object.
  - `feature`: Gene to plot.
  - `pt.size`, `plot.margin`: Customization parameters.
- **Examples:** `modify_vlnplot(obj, feature = "GeneX")`

---

## 8. `StackedVlnPlot`
**Description:** Creates stacked violin plots for multiple features.
- **Parameters:**
  - `obj`: Seurat object.
  - `features`: List of features to plot.
- **Examples:** `StackedVlnPlot(obj, features = c("Gene1", "Gene2"))`

---

## 9. `PercExp`
**Description:** Calculates the percentage of cells expressing specific genes across groups.
- **Parameters:**
  - `object`: Seurat object.
  - `features`: Genes of interest.
  - `group_by`, `split_by`: Metadata variables for grouping.
- **Examples:** `PercExp(object, features = c("GeneX"))`

---

## 10. `sexDif`
**Description:** Analyzes sex-based differences in gene expression proportions.
- **Parameters:**
  - `object`: Seurat object.
  - `gene`: Gene of interest.
  - `plot`: Whether to generate plots.
- **Examples:** `sexDif(object, gene = "GeneX", plot = TRUE)`

---

## 11. `euler_plot`
**Description:** Generates Euler diagrams to visualize gene expression overlaps.
- **Parameters:**
  - `object`: Seurat object.
  - `genes`: List of genes.
  - `metacol`, `metavar`: Metadata columns for subsetting.
- **Examples:** `euler_plot(object, genes = c("Gene1", "Gene2"))`

---

## 12. `upset_plot`
**Description:** Creates UpSet plots for gene expression overlaps.
- **Parameters:**
  - `object`: Seurat object.
  - `genes`: List of genes.
- **Examples:** `upset_plot(object, genes = c("GeneX", "GeneY"))`

---

## 13. `difAlignment`
**Description:** Tests alignment effects on gene expression proportions.
- **Parameters:**
  - `object`: Seurat object.
  - `gene`: Gene of interest.
  - `assay1`, `assay2`: Assays to compare.
- **Examples:** `difAlignment(object, gene = "GeneX", assay1 = "RNA", assay2 = "SCT")`

---

## 14. `IncMorphDif`
**Description:** Analyzes pain and morphine data for expression and effect sizes.
- **Parameters:**
  - `object`: Seurat object.
  - `module`: Module score metadata column.
  - `gene_list`: List of genes in the module.
- **Examples:** `IncMorphDif(object, module = "signature_1", gene_list = c("Gene1", "Gene2"))`
