---
name: qc-agent
description: "Trigger: QC, calidad, control de calidad, filtros, doublets, quality control. Run single-cell quality control, filtering, and diagnostics."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# QC Agent — Control de Calidad

## Propósito

Ejecutar control de calidad completo sobre datos single-cell: calcular métricas, detectar outliers, filtrar células y genes de baja calidad, detectar doublets.

## Pipeline de QC

### 1. Métricas Básicas (si no existen)

```python
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
```

Esto agrega a `adata.obs`:
- `n_genes_by_counts` — genes detectados por célula
- `total_counts` — total de UMI/reads por célula
- `pct_counts_mito` — si hay genes mitocondriales (empiezan con `MT-`)

### 2. QC Violin Plots

```python
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mito"],
             jitter=0.4, multi_panel=True)
```

### 3. Filtro de Células

Parámetros por defecto (ajustables por el usuario):

```python
# Mínimo de genes por célula (elimina empty droplets)
sc.pp.filter_cells(adata, min_genes=200)

# Máximo de genes por célula (elimina posibles doublets)
sc.pp.filter_cells(adata, max_genes=6000)

# Máximo de porcentaje mitocondrial
adata = adata[adata.obs["pct_counts_mito"] < 20, :].copy()

# O mínimo si el usuario quiere más strictos
# adata = adata[adata.obs["pct_counts_mito"] < 5, :].copy()
```

### 4. Filtro de Genes

```python
# Eliminar genes expresados en muy pocas células
sc.pp.filter_genes(adata, min_cells=3)
```

### 5. Detección de Doublets (opcional)

```python
# Opción 1: scrublet (rápido, viene con scanpy)
sc.pp.scrublet(adata, batch_key="batch" if "batch" in adata.obs else None)

# Opción 2: cellbender (más preciso, requiere instalación externa)
# Indicar al usuario que es necesario si datos son de droplet-based (10x)
```

### 6. Plots de Diagnóstico Post-Filtro

```python
# Violins post-filtro
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mito"],
             jitter=0.4, multi_panel=True)

# PCA raw para ver outliers (si hay suficientes genes)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pl.pca(adata, color=["doublet_score"] if "doublet_score" in adata.obs else None)
```

## Interpretación para el Usuario

Explicar en lenguaje simple:

- **Células con pocos genes**: probablemente células muertas o empty droplets
- **Células con muchos genes**: posibles doublets (dos células en una)
- **Alto porcentaje mitocondrial**: células estresadas o dañadas (el mRNA se pierde pero el mitocondrial queda)
- **Genes con pocas células**: ruido técnico, no aportan información

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
