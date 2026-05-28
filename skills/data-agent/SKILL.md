---
name: data-agent
description: "Trigger: cargar datos, h5ad, h5mu, zarr, AnnData, MuData, SpatialData. Load and describe single-cell datasets for analysis."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Data Agent — Data Loading and Description

## Purpose

Load single-cell data from files and generate a complete dataset description: dimensions, layers, metadata, basic quality.

## Supported Formats

| Formato | Librería | Objeto | Uso |
|---------|----------|--------|-----|
| .h5ad | `scanpy` / `anndata` | `AnnData` | RNA-seq unimodal |
| .h5mu | `mudata` | `MuData` | Multimodal (RNA+ATAC+proteína) |
| .zarr | `spatialdata` | `SpatialData` | Datos espaciales (Xenium, Visium) |
| .h5ad + imágenes | `scanpy` + `squidpy` | `AnnData` + `sdata` | Espacial legacy |

## APIs Clave

```python
import scanpy as sc
import mudata as md
import spatialdata as sd

# AnnData
adata = sc.read_h5ad("path/to/data.h5ad")

# MuData
mdata = md.read("path/to/data.h5mu")

# SpatialData
sdata = sd.read_zarr("path/to/data.zarr")
```

## Dataset Description

Always report:

- **Dimensions**: `n_obs × n_vars` (cells × genes)
- **Layers**: what's in `.layers` (counts, normalized, etc.)
- **Obs columns**: available cell metadata
- **Var columns**: available gene metadata (if applicable)
- **Embeddings**: existing reductions (PCA, UMAP, etc.)
- **Uns**: content of `.uns` (marker genes, colors, etc.)
- **Raw**: whether `.raw` is saved
- **For MuData**: modalities, dimensions of each
- **For SpatialData**: elements (Images, Labels, Points, Shapes, Tables)

## Output

Devolver al orquestador:

```python
{
    "status": "ok",
    "data_type": "AnnData | MuData | SpatialData",
    "data_path": "/path/to/saved/data.h5ad",
    "summary": {
        "n_obs": 10000,
        "n_vars": 20000,
        "layers": ["counts", "normalized"],
        "obs_columns": ["n_genes", "total_counts", "pct_mito"],
        "has_raw": True,
        "modalities": None,  # o ["rna", "atac"] para MuData
    },
    "recommendations": [
        "Tiene pct_mito → QC con filtro mitocondrial recomendado",
        "No tiene HVG → correr normalize-agent",
    ]
}
```

## Hard Rules

- **NEVER modify, overwrite, or delete the original input file**. Open it in read-only mode by default.
- If the file doesn't exist, report immediately.
- If the format is unrecognized, suggest conversion.
- Save a copy of the loaded object to `output/{session_id}/` for the next agents to read.
- The only person who removes original data is the human from the CLI.
