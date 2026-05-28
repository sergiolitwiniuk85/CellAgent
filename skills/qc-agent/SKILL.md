---
name: qc-agent
description: "Trigger: QC, calidad, control de calidad, filtros, doublets, quality control. Run single-cell quality control, filtering, and diagnostics."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# QC Agent — Quality Control

## Purpose

Run complete quality control on single-cell data: compute metrics, detect outliers, filter low-quality cells and genes, detect doublets.

## QC Pipeline

### 1. Basic Metrics (if not already present)

```python
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
```

This adds to `adata.obs`:
- `n_genes_by_counts` — genes detected per cell
- `total_counts` — total UMI/reads per cell
- `pct_counts_mito` — if mitochondrial genes exist (starting with `MT-`)

### 2. QC Violin Plots

```python
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mito"],
             jitter=0.4, multi_panel=True)
```

### 3. Cell Filtering

Default parameters (adjustable by the user):

```python
# Minimum genes per cell (removes empty droplets)
sc.pp.filter_cells(adata, min_genes=200)

# Maximum genes per cell (removes potential doublets)
sc.pp.filter_cells(adata, max_genes=6000)

# Maximum mitochondrial percentage
adata = adata[adata.obs["pct_counts_mito"] < 20, :].copy()

# Or stricter if the user prefers
# adata = adata[adata.obs["pct_counts_mito"] < 5, :].copy()
```

### 4. Gene Filtering

```python
# Remove genes expressed in very few cells
sc.pp.filter_genes(adata, min_cells=3)
```

### 5. Doublet Detection (optional)

```python
# Option 1: scrublet (fast, ships with scanpy)
sc.pp.scrublet(adata, batch_key="batch" if "batch" in adata.obs else None)

# Option 2: cellbender (more accurate, requires external install)
# Flag to user that this is needed if data is droplet-based (10x)
```

### 6. Post-Filter Diagnostic Plots

```python
# Violins post-filtro
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mito"],
             jitter=0.4, multi_panel=True)

# PCA raw para ver outliers (si hay suficientes genes)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pl.pca(adata, color=["doublet_score"] if "doublet_score" in adata.obs else None)
```

## Interpretation for the User

Explain in plain language:

- **Low gene count cells**: likely dead cells or empty droplets
- **High gene count cells**: potential doublets (two cells in one)
- **High mitochondrial percentage**: stressed or damaged cells (mRNA is lost but mitochondrial remains)
- **Genes in few cells**: technical noise, not informative

## Output

Devolver al orquestador:

```python
{
    "stage": "qc",
    "status": "completed",
    "summary": f"QC completado. Retenidas {n_cells_after} de {n_cells_before} células ({pct_retained:.1f}%)",
    "params": {
        "min_genes": 200,
        "max_genes": 6000,
        "max_pct_mito": 20,
        "min_cells": 3,
        "doublet_method": "scrublet",
    },
    "qc_metrics": {
        "cells_before": n_cells_before,
        "cells_after": n_cells_after,
        "median_genes": float(adata.obs["n_genes_by_counts"].median()),
        "median_counts": float(adata.obs["total_counts"].median()),
        "pct_mito_mean": float(adata.obs["pct_counts_mito"].mean()),
        "doublets_detected": n_doublets if "doublet_score" in adata.obs else None,
    },
    "plots": ["qc_violin_pre.png", "qc_violin_post.png", "pca_raw.png"],
    "data_path": "/path/to/filtered_data.h5ad",
    "recommendations": [
        "Porcentaje mitocondrial bajo → células de buena calidad",
        f"Doublets: {n_doublets} detectados",
        "Listo para normalización",
    ]
}
```

## Hard Rules

- **NEVER filter without showing diagnostic plots first**
- **NEVER delete, overwrite, or modify the original input file**. The filtered data is always written to a NEW path.
- Default thresholds are for 10x PBMC data. Ask if it's a different technology.
- Save the filtered object to a NEW path: `output/{session_id}/filtered.h5ad`
- The only person who removes original data is the human from the CLI.
