---
name: cluster-agent
description: "Trigger: clustering, cluster, PCA, UMAP, leiden, marker genes, reducción de dimensionalidad. Run PCA, neighbors, embedding, clustering, and marker gene detection."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Cluster Agent — Reducción, Clustering y Marker Genes

## Propósito

Ejecutar el pipeline completo de reducción de dimensionalidad, clustering y detección de marker genes sobre datos single-cell normalizados.

## Pipeline

### 1. PCA

```python
# Siempre correr PCA
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

# Mostrar varianza explicada (ayuda a decidir n_comps)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# Scatter de primeras componentes
sc.pl.pca(adata, color=["n_genes_by_counts", "total_counts", "pct_counts_mito"])
```

### 2. Neighborhood Graph

```python
sc.pp.neighbors(
    adata,
    n_neighbors=15,           # estándar
    n_pcs=30,                 # usar 30 PCs (ajustable)
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
    resolution=0.8,           # resolución por defecto
    key_added="leiden",
    random_state=42,
)
```

**Estrategia de resolución**:
- 0.1-0.3 → clusters grandes (pocos, biológicos amplios)
- 0.5-1.0 → resolución estándar (recomendado arrancar con 0.8)
- 1.0-2.0 → clusters finos (muchos, para subpoblaciones)

Mostrar UMAP coloreado por clusters y por métricas de QC:

```python
sc.pl.umap(adata, color=["leiden", "n_genes_by_counts", "total_counts", "pct_counts_mito"],
           legend_loc="on data", ncols=2)
```

### 5. Marker Genes

```python
# Wilcoxon rank-sum test (rápido, default)
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    n_genes=50,
    key_added="rank_genes_groups",
)
```

Métodos disponibles:
- `wilcoxon` → rápido, default, recomendado (rank-sum test)
- `t-test` → más simple, asume normalidad
- `logreg` → logistic regression, más lento pero robusto

Mostrar:

```python
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, key="rank_genes_groups")

# Dotplot de top markers por cluster
sc.pl.dotplot(adata, var_names=sc.get.rank_genes_groups_df(adata, group=None)
              .groupby("group").head(3)["names"].tolist(),
              groupby="leiden")
```

### 6. Opcional: Anotación Automática

Si el usuario quiere anotar tipos celulares:

```python
# Buscar marcadores conocidos en los top genes de cada cluster
# Sugerir al usuario revisar estos genes contra literatura conocida

# Si hay CellTypist instalado, ofrecer anotación automática:
# import celltypist
# model = celltypist.models.Model.load(model="Immune_All_Low.pkl")
# predictions = celltypist.annotate(adata, model=model, majority_voting=True)
```

## Validación

```python
print(f"Clusters encontrados: {adata.obs['leiden'].nunique()}")
print(f"Resolución usada: {resolution}")
sc.pl.umap(adata, color="leiden", legend_loc="on data")

# Mostrar top 5 genes por cluster
top5 = sc.get.rank_genes_groups_df(adata, group=None).groupby("group").head(5)
print(top5[["group", "names", "scores", "pvals_adj"]].to_string())
```

## Output

```python
{
    "stage": "cluster",
    "status": "completed",
    "summary": f"Clustering completado: {n_clusters} clusters (leiden, res={resolution})",
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
        f"{n_clusters} clusters identificados",
        "Se detectaron marker genes significativos (p_adj < 0.05)",
        "¿Querés guardar los resultados o seguir a interpretación?",
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
