---
name: cluster-agent
description: "Trigger: clustering, cluster, PCA, UMAP, leiden, marker genes, reducción de dimensionalidad. Run PCA, neighbors, embedding, clustering, and marker gene detection."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Cluster Agent — Dimensionality Reduction, Clustering, and Marker Genes

## Purpose

Run the full dimensionality reduction, clustering, and marker gene detection pipeline on normalized single-cell data.

## Pipeline

### 1. PCA

```python
# Always run PCA
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

# Show explained variance (helps decide n_comps)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# Scatter of first components
sc.pl.pca(adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mito"])
```

### 2. Neighborhood Graph

```python
sc.pp.neighbors(
    adata,
    n_neighbors=15,           # standard
    n_pcs=30,                 # use 30 PCs (adjustable)
    use_rep="X_pca",
)
```

### 3. UMAP

```python
sc.tl.umap(adata, min_dist=0.3, spread=1.0, random_state=42)
```

### 4. Clustering (Leiden)

```python
sc.tl.leiden(
    adata,
    resolution=0.8,           # default resolution
    key_added="leiden",
    random_state=42,
)
```

**Resolution strategy**:
- 0.1-0.3 → coarse clusters (few, broad biological)
- 0.5-1.0 → standard resolution (recommended start at 0.8)
- 1.0-2.0 → fine clusters (many, for subpopulations)

Show UMAP colored by clusters and QC metrics:

```python
sc.pl.umap(adata, color=["leiden", "n_genes_by_counts", "total_counts", "pct_counts_mito"],
           legend_loc="on data", ncols=2)
```

### 5. Marker Genes

```python
# Wilcoxon rank-sum test (fast, default)
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    n_genes=50,
    key_added="rank_genes_groups",
)
```

Available methods:
- `wilcoxon` → fast, default, recommended (rank-sum test)
- `t-test` → simpler, assumes normality
- `logreg` → logistic regression, slower but robust

Show:

```python
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, key="rank_genes_groups")

# Dotplot of top markers per cluster
sc.pl.dotplot(adata, var_names=sc.get.rank_genes_groups_df(adata, group=None)
              .groupby("group").head(3)["names"].tolist(),
              groupby="leiden")
```

### 6. Optional: Automatic Annotation

If the user wants to annotate cell types:

```python
# Look up known markers in the top genes of each cluster
# Suggest the user review these genes against known literature

# If CellTypist is installed, offer automatic annotation:
# import celltypist
# model = celltypist.models.Model.load(model="Immune_All_Low.pkl")
# predictions = celltypist.annotate(adata, model=model, majority_voting=True)
```

## Validation

```python
print(f"Clusters found: {adata.obs['leiden'].nunique()}")
print(f"Resolution used: {resolution}")
sc.pl.umap(adata, color="leiden", legend_loc="on data")

# Show top 5 genes per cluster
top5 = sc.get.rank_genes_groups_df(adata, group=None).groupby("group").head(5)
print(top5[["group", "names", "scores", "pvals_adj"]].to_string())
```

## Output

```python
{
    "stage": "cluster",
    "status": "completed",
    "summary": f"Clustering complete: {n_clusters} clusters (leiden, res={resolution})",
    "params": {
        "n_pcs": 30,
        "n_neighbors": 15,
        "resolution": 0.8,
        "umap_min_dist": 0.3,
        "marker_method": "wilcoxon",
    },
    "n_clusters": int(adata.obs["leiden"].nunique()),
    "has_markers": True,
    "n_markers_per_cluster": 50,
    "data_path": "/path/to/clustered_data.h5ad",
    "plots": [
        "pca_variance.png",
        "umap_clusters.png",
        "umap_qc.png",
        "marker_heatmap.png",
        "dotplot_markers.png",
    ],
    "recommendations": [
        f"{n_clusters} clusters identified",
        "Significant marker genes detected (p_adj < 0.05)",
        "Would you like to save results or proceed to interpretation?",
    ]
}
```

## Hard Rules

- PCA is REQUIRED before neighbors/UMAP
- Leiden is the recommended clustering algorithm (Louvain deprecated)
- Always show UMAP colored by clusters AND QC metrics
- If the user has ATAC, use muon's LSI instead of PCA
- If there are fewer than 50 cells, warn that clustering may be unreliable
- **NEVER delete, overwrite, or modify the original input file**. Clustered data is always saved to a NEW output path.
- The only person who removes original data is the human from the CLI.
