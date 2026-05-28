---
name: report-agent
description: "Trigger: reporte, report, resumen, summary, informe, executive summary. Generate an executive report of completed single-cell analysis."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Report Agent — Report Generation

## Purpose

Synthesize the entire single-cell analysis into an executive report that the scientist can read, share, and use in publications.

## Output Format

The report is generated as:

1. **Markdown** → readable report with embedded plots (base64 or paths)
2. **PDF** → publication-ready via Typst (lightweight, modern typesetting)
3. **HTML** (optional) → prettier, with interactive tabs
4. **Chat summary** → the orchestrator shows it to the user

### PDF Generation (Typst)

After writing the Markdown report, convert to PDF using the bundled script:

```bash
# Basic usage
python skills/report-agent/scripts/md_to_pdf.py report.md -o report.pdf

# With custom margins
python skills/report-agent/scripts/md_to_pdf.py report.md -o report.pdf --margin 1.5in

# Keep intermediate .typ file (debugging)
python skills/report-agent/scripts/md_to_pdf.py report.md -o report.pdf --keep-typ
```

The pipeline is: **Markdown → Pandoc → Typst → PDF**

- **Pandoc** converts MD → Typst format
- **Typst** compiles to PDF (lightweight, ~15 MB binary vs ~500 MB for LaTeX)

**Requirements**: `pandoc >= 3.1` (with typst writer) and `typst >= 0.12`.
Both are included in the Docker image or can be installed via:

```bash
# Via conda
conda install -c conda-forge pandoc
# Typst binary (Linux x86_64)
curl -sSfL https://github.com/typst/typst/releases/latest/download/typst-x86_64-unknown-linux-musl.tgz | tar xz
sudo mv typst /usr/local/bin/
```

> **Why Typst over LaTeX?** Typst is ~30x smaller, compiles 5-10x faster, has cleaner syntax, and produces equally professional PDFs.

### Script reference

- `scripts/md_to_pdf.py` — Markdown → PDF converter (uses pandoc + typst)
- `scripts/scai_template.typ` — Typst report template (for standalone use)

## Report Structure

```markdown
# Single-Cell Analysis Report
## Project: {name}
## Date: {date}

### 1. Executive Summary
- {n_cells} cells analyzed
- {n_genes} genes detected
- {n_clusters} clusters identified
- {modalities} (if applicable)

### 2. Quality Control
- Cells before filtering: {n}
- Cells after filtering: {n} ({pct}% retained)
- Doublets detected: {n}
- ![QC Violins](qc_violin_post.png)

### 3. Normalization
- Method: scanpy (normalize_total + log1p)
- HVGs selected: {n} of {n} total
- ![HVG Plot](hvg_plot.png)

### 4. Dimensionality Reduction
- PCA: {n} components
- ![PCA Variance](pca_variance.png)

### 5. Clustering
- Algorithm: Leiden, resolution {res}
- Clusters found: {n}
- ![UMAP Clusters](umap_clusters.png)

### 6. Marker Genes
- Method: {method}
- Top genes per cluster:
  | Cluster | Gene 1 | Gene 2 | Gene 3 |
  |---------|--------|--------|--------|
  | 0       | CD3D   | IL7R   | CCR7   |
  | 1       | CD14   | LYZ    | CST3   |
  | ...     | ...    | ...    | ...    |
- ![Marker Heatmap](marker_heatmap.png)

### 7. Preliminary Interpretation
- Cluster 0: appears to be {cell_type} (expresses {markers})
- Cluster 1: appears to be {cell_type}
- {recommendations for next steps}

### 8. Image-Derived Features

- Cells with image features: {n_cells}
- Images processed: {images}
- Total features extracted: {n_features}
- Levels completed: {levels}
- Tissue compartments found: {compartments}

**Image Feature Summary**:
| Feature | Mean | Std | Min | Max |
|---------|------|-----|-----|-----|
| {feature_name} | {mean} | {std} | {min} | {max} |
| ... | ... | ... | ... | ... |

**Top-5 Most Variable Features** (by coefficient of variation, with spatial scatter):
- ![Spatial Scatter](image_top5_spatial_scatter.png)

**PCA Variance Explained**:
- ![PCA Variance](image_pca_variance.png)

**PC1 Spatial Overlay** (first principal component mapped onto tissue coordinates):
- ![PC1 Spatial](image_pc1_spatial.png)

### 9. Parameters Used
| Stage | Parameter | Value |
|-------|-----------|-------|
| QC | min_genes | 200 |
| QC | max_pct_mito | 20 |
| Normalization | n_top_genes | 2000 |
| Clustering | resolution | 0.8 |
| ... | ... | ... |
```

## How to Build the Report

The report-agent receives the full `PipelineState` from the orchestrator with all history and results. Simply:

1. Iterate `state["history"]` for each completed stage
2. Extract metrics, parameters, and plot paths
3. Generate the markdown
4. If marker genes exist, convert to table
5. Include basic cell type interpretation if possible

## Cell Type Interpretation (Basic)

Use a lookup table of classic markers (the most common):

```python
# Classic markers for human PBMC
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

For each cluster: intersect top markers with this table and suggest cell type.

## Output

```python
{
    "stage": "report",
    "status": "completed",
    "report_path": "/path/to/report.md",
    "summary": "Report generated with {n_clusters} clusters, {n_plots} plots",
    "cell_type_annotations": {
        "0": "CD4+ T cell (CD3D+, CD4+)",
        "1": "Monocyte (CD14+, LYZ+)",
        ...
    },
    "plots_in_report": n_plots,
    "recommendations": [
        "Review cluster 3 annotation — doesn't clearly match known markers",
        "Consider differential analysis between conditions if group metadata is available",
    ]
}
```

## Hard Rules

- Do not fabricate biological interpretations. If there aren't enough markers, say so.
- Always include the parameter table (reproducibility)
- If the pipeline has fewer than 3 stages, report as "partial analysis"
- Do not mention "significant" without valid statistical tests
- **NEVER delete, overwrite, or modify any data file**. Reports are write-only — they don't touch original data.
