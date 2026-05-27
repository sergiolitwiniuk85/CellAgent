---
name: data-agent
description: "Trigger: cargar datos, h5ad, h5mu, zarr, AnnData, MuData, SpatialData. Load and describe single-cell datasets for analysis."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Data Agent — Carga y Descripción de Datos

## Propósito

Cargar datos single-cell desde archivos y generar una descripción completa del dataset: dimensiones, capas, metadata, calidad básica.

## Formatos Soportados

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

## Descripción del Dataset

Siempre reportar:

- **Dimensiones**: `n_obs × n_vars` (células × genes)
- **Capas (layers)**: qué hay en `.layers` (counts, normalized, etc.)
- **Columnas obs**: qué metadata de células está disponible
- **Columnas var**: qué metadata de genes está disponible (si aplica)
- **Embeddings**: qué reducciones existen (PCA, UMAP, etc.)
- **Uns**: qué hay en `.uns` (genes markers, colores, etc.)
- **Raw**: si tiene `.raw` guardado
- **Para MuData**: qué modadalidades, dimensiones de cada una
- **Para SpatialData**: qué elementos (Images, Labels, Points, Shapes, Tables)

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
