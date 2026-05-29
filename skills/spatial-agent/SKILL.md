---
name: spatial-agent
description: "Trigger: espacial, spatial, squidpy, spatialdata, xenium, visium, imágenes, coordenadas. Run spatial omics analysis with Squidpy and SpatialData."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Spatial Agent — Spatial Analysis

## Purpose

Analyze spatial data (Xenium, Visium, MERFISH) using SpatialData for structure and Squidpy for analysis.

## Input Format

```python
import spatialdata as sd
import squidpy as sq

# SpatialData (formato moderno)
sdata = sd.read_zarr("path/to/data.zarr")

# AnnData legacy (formato squidpy clásico)
adata = sc.read_h5ad("path/to/spatial_data.h5ad")
```

## SpatialData Elements

| Element | Access | Description |
|----------|--------|-------------|
| Images | `sdata["image_name"]` | H&E, DAPI, staining |
| Labels | `sdata["label_name"]` | Pixel-level segmentation |
| Points | `sdata["transcripts"]` | Transcript coordinates |
| Shapes | `sdata["boundaries"]` | Cell/nuclei polygons |
| Tables | `sdata["table"]` | AnnData with gene expression |

## Spatial Analysis Pipeline

### 1. Dataset Exploration

```python
print(sdata)

# Mostrar elementos
sdata.pl.render_images("raw_image").pl.show()
sdata.pl.render_shapes("nucleus_boundaries").pl.show()
sdata.pl.render_points("transcripts", color="gene", groups="Vwf").pl.show()
```

### 2. Spatial Neighbors

```python
# Calcular vecinos espaciales
sq.gr.spatial_neighbors(sdata["table"], coord_type="grid")  # Visium
# o
sq.gr.spatial_neighbors(sdata["table"], coord_type="generic")  # Xenium
```

### 3. Análisis de Enriquecimiento de Vecindad

```python
# Cluster enrichment en espacio
sq.gr.nhood_enrichment(sdata["table"], cluster_key="leiden")
sq.pl.nhood_enrichment(sdata["table"], cluster_key="leiden")
```

### 4. Visualización Espacial

```python
# Scatter espacial de clusters
sq.pl.spatial_scatter(sdata["table"], shape=None, color="leiden")

# O usando spatialdata-plot
sdata.pl.render_shapes("nucleus_boundaries", color="leiden").pl.show()
```

### 5. Análisis de Interacción (opcional)

```python
# Co-ocurrencia de tipos celulares
sq.gr.co_occurrence(sdata["table"], cluster_key="leiden")
sq.pl.co_occurrence(sdata["table"], cluster_key="leiden")

# Ligando-receptor (con datos de tissue)
# sq.gr.ligand_receptor(sdata["table"], ...)
```

### 6. Spatial Autocorrelation (after clustering)

```python
# Prerequisite: compute spatial neighbors graph from coordinates
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

# Moran's I — global spatial autocorrelation
sq.gr.spatial_autocorr(adata, mode="moran")

# Geary's C — local spatial variation (more sensitive)
sq.gr.spatial_autocorr(adata, mode="geary")

# Results stored in adata.uns
print(adata.uns['moranI'].head())
print(adata.uns['gearyC'].head())
```

**Output structure** (`adata.uns['moranI']` and `adata.uns['gearyC']`):

| Column | Description |
|--------|-------------|
| `I` / `C` | Moran's I / Geary's C statistic |
| `pval_norm` | p-value under normality assumption |
| `var_norm` | Variance under normality |
| `pval_norm_fdr_bh` | Benjamini-Hochberg corrected p-value |

**Interpretation**:
- **Moran's I**: Positive values → clustering (similar values near each other). Values near zero → random spatial distribution. Negative → dispersion.
- **Geary's C**: Values < 1 → positive spatial autocorrelation. Values > 1 → negative autocorrelation. More sensitive to local variation than Moran's I.

## Validación

```python
print(f"Elementos: {list(sdata.attr_keys())}")
n_cells = sdata["table"].n_obs
print(f"Células en tabla: {n_cells}")
print(f"Clusters: {sdata['table'].obs['leiden'].nunique()}")
```

## Output

```python
{
    "stage": "spatial",
    "status": "completed",
    "summary": f"Análisis espacial completado: {n_cells} células, {n_clusters} clusters",
    "n_elements": len(sdata.attr_keys()),
    "has_images": True,
    "n_cells": n_cells,
    "data_path": "/path/to/spatial_analysis.zarr",
    "plots": ["tissue_image.png", "spatial_clusters.png", "nhood_enrichment.png"],
    "recommendations": [
        "Se observa segregación espacial de clusters",
        "Revisar interacciones ligando-receptor entre regiones",
    ]
}
```

## Hard Rules

- SpatialData is the modern format; Squidpy is for analysis, not storage
- Spatial data is LARGE — recommend backed mode when possible
- Do not process images (segmentation, etc.) unless the user explicitly asks
- For Visium: remember resolution is spot-level, not single-cell
- For Xenium: data is already single-cell, segmentation is included
- **NEVER delete, overwrite, or modify the original input file**. All results go to a NEW output path.
- The only person who removes original data is the human from the CLI.
