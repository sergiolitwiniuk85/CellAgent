---
name: report-agent
description: "Trigger: reporte, report, resumen, summary, informe, executive summary. Generate an executive report of completed single-cell analysis."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Report Agent — Generación de Reporte

## Propósito

Sintetizar todo el análisis single-cell en un reporte ejecutivo que el científico pueda leer, compartir y usar en publicaciones.

## Formato de Salida

El reporte se genera como:

1. **Markdown** → informe legible, con plots embedidos (base64 o paths)
2. **HTML** (opcional) → más lindo, con tabs interactivos
3. **Resumen en chat** → el orquestador lo muestra al usuario

## Estructura del Reporte

```markdown
# Reporte de Análisis Single-Cell
## Proyecto: {nombre}
## Fecha: {fecha}

### 1. Resumen Ejecutivo
- {n_células} células analizadas
- {n_genes} genes detectados
- {n_clusters} clusters identificados
- {modadalidades} (si aplica)

### 2. Control de Calidad
- Células antes del filtro: {n}
- Células después del filtro: {n} ({pct}% retenidas)
- Dobletes detectados: {n}
- ![QC Violins](qc_violin_post.png)

### 3. Normalización
- Método: scanpy (normalize_total + log1p)
- HVG seleccionados: {n} de {n} totales
- ![HVG Plot](hvg_plot.png)

### 4. Reducción de Dimensionalidad
- PCA: {n} componentes
- ![PCA Variance](pca_variance.png)

### 5. Clustering
- Algoritmo: Leiden, resolución {res}
- Clusters encontrados: {n}
- ![UMAP Clusters](umap_clusters.png)

### 6. Marker Genes
- Método: {method}
- Top genes por cluster:
  | Cluster | Gene 1 | Gene 2 | Gene 3 |
  |---------|--------|--------|--------|
  | 0       | CD3D   | IL7R   | CCR7   |
  | 1       | CD14   | LYZ    | CST3   |
  | ...     | ...    | ...    | ...    |
- ![Marker Heatmap](marker_heatmap.png)

### 7. Interpretación Preliminar
- Cluster 0: parece ser {cell_type} (expresa {markers})
- Cluster 1: parece ser {cell_type}
- {recomendaciones para el próximo paso}

### 8. Parámetros Usados
| Etapa | Parámetro | Valor |
|-------|-----------|-------|
| QC | min_genes | 200 |
| QC | max_pct_mito | 20 |
| Normalización | n_top_genes | 2000 |
| Clustering | resolution | 0.8 |
| ... | ... | ... |
```

## Cómo Construir el Reporte

El report-agent recibe el `PipelineState` completo del orquestador con toda la historia y resultados. Simplemente:

1. Iterar `state["history"]` para cada etapa completada
2. Extraer métricas, parámetros y paths de plots
3. Generar el markdown
4. Si hay marker genes, convertir a tabla
5. Incluir interpretación básica de tipos celulares si es posible

## Interpretación de Tipos Celulares (Básica)

Usar una lookup table de marcadores clásicos (la más común):

```python
# Marcadores clásicos para PBMC humanos
cell_type_markers = {
    "T cell": ["CD3D", "CD3E", "CD7"],
    "CD4+ T cell": ["CD3D", "CD4", "IL7R"],
    "CD8+ T cell": ["CD3D", "CD8A", "CD8B"],
    "NK cell": ["NKG7", "GNLY", "KLRD1"],
    "Monocyte": ["CD14", "LYZ", "CST3"],
    "Macrophage": ["CD68", "CD163", "CSF1R"],
    "B cell": ["MS4A1", "CD79A", "CD19"],
    "Plasma cell": ["MZB1", "SDC1", "JCHAIN"],
    "Dendritic cell (DC)": ["FCER1A", "CLEC10A", "CST3"],
    "pDC": ["IL3RA", "CLEC4C", "TCF4"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8"],
}
```

Para cada cluster: intersectar top markers con esta tabla y sugerir tipo celular.

## Output

```python
{
    "stage": "report",
    "status": "completed",
    "report_path": "/path/to/report.md",
    "summary": "Reporte generado con {n_clusters} clusters, {n_plots} plots",
    "cell_type_annotations": {
        "0": "CD4+ T cell (CD3D+, CD4+)",
        "1": "Monocyte (CD14+, LYZ+)",
        ...
    },
    "plots_in_report": n_plots,
    "recommendations": [
        "Revisar anotación de cluster 3 — no coincide claramente con marcadores conocidos",
        "Considerar análisis diferencial entre condiciones si hay metadata de grupos",
    ]
}
```

## Hard Rules

- Do not fabricate biological interpretations. If there aren't enough markers, say so.
- Always include the parameter table (reproducibility)
- If the pipeline has fewer than 3 stages, report as "partial analysis"
- Do not mention "significant" without valid statistical tests
- **NEVER delete, overwrite, or modify any data file**. Reports are write-only — they don't touch original data.
