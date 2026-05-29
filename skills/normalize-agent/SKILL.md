---
name: normalize-agent
description: "Trigger: normalizar, normalización, HVG, highly variable genes, scaling, log1p, normalize. Run normalization and feature selection for single-cell data."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Normalize Agent — Normalization and Feature Selection

## Purpose

Run standard scanpy normalization: counts per cell, log-transform, highly variable gene (HVG) selection, and scaling.

## Normalization Pipeline

### 1. Save raw counts (if not already a layer)

```python
if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
```

### 2. Normalize per cell

```python
# Normalize to 10,000 counts per cell (standard)
sc.pp.normalize_total(adata, target_sum=1e4)
```

### 3. Log-transform

```python
sc.pp.log1p(adata)
```

### 4. HVG Selection

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,           # standard for human/mouse data
    batch_key=None,             # if batch exists, pass it for regression
    flavor="seurat_v3" if "counts" in adata.layers else "seurat",
    subset=False,               # don't subset, just mark in .var
)
```

**Note**: `seurat_v3` is more accurate but requires the `counts` layer. Fall back to `seurat` if unavailable.

Show:

```python
sc.pl.highly_variable_genes(adata)
```

### 5. Optional: Batch effect regression (if batch exists)

```python
if "batch" in adata.obs and adata.obs["batch"].nunique() > 1:
    # Inform the user about batch effects
    # Ask if they want to correct: regress_out or combat
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mito"])
```

### 6. Scaling

```python
sc.pp.scale(adata, max_value=10)
```

## Recommended Parameters by Context

| Context | n_top_genes | target_sum | Notes |
|----------|-------------|------------|-------|
| 10x PBMC (human) | 2000 | 1e4 | Standard |
| Mouse data | 2000 | 1e4 | Similar to human |
| Smart-seq2 | All | 1e6 | Full-length, not 10x |
| ATAC | 0 (N/A) | N/A | Use LSI (muon.atac) |
| Protein (CITE) | 0 (N/A) | N/A | CLR transform |
| Targeted panels (<1000 genes) | All (skip HVG) | 1e4 | Xenium, CosMx panels (250–500 genes) |

## Validation

Show after normalizing:

```python
# Check distribution
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts"])
print(f"Total genes: {adata.n_vars}")
print(f"HVG selected: {adata.var['highly_variable'].sum()}")
```

## Output

```python
{
    "stage": "normalize",
    "status": "completed",
    "summary": f"Normalized. {n_hvg} HVGs from {n_vars} total genes",
    "params": {
        "target_sum": 1e4,
        "log1p": True,
        "n_top_genes": 2000,
        "hvg_flavor": "seurat_v3",
        "scaled": True,
        "max_value": 10,
    },
    "n_hvg": int(adata.var["highly_variable"].sum()),
    "data_path": "/path/to/normalized_data.h5ad",
    "plots": ["hvg_plot.png"],
    "recommendations": [
        f"{n_hvg} HVGs selected → ready for clustering",
        "Run PCA with n_comps=50",
    ]
}
```

## Hard Rules

- **NEVER delete, overwrite, or modify the original input file**. All normalized data goes to a NEW output path.
- Always save `counts` to `adata.layers` before normalizing
- If the user has protein data (CITE-Seq), do NOT use normalize_total — use muon's CLR instead
- If the user has ATAC data, redirect to `muon.atac.pp` (LSI instead of PCA)
- Do not subset the object to only HVG unless the user explicitly asks
- Save a copy before scaling (scaling is irreversible)
- The only person who removes original data is the human from the CLI.
- **If `adata.n_vars < 1000`**: use ALL genes, do NOT call `highly_variable_genes` (targeted panel). Only select HVGs when `adata.n_vars >= 1000`.
