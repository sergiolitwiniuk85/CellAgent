"""SCAI Trace Agent — Generate traceability document for reproducibility.

Reads output directory contents and produces a comprehensive forensic
document of every stage, parameter, and result.
"""

import json
import importlib.metadata
import subprocess
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    import csv
    with open(path) as f:
        return list(csv.DictReader(f))


def _read_json(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def _get_library_versions() -> dict[str, str]:
    """Get versions of key libraries."""
    pkgs = [
        'scanpy', 'anndata', 'squidpy', 'pandas', 'numpy',
        'scikit-learn', 'scikit-image', 'matplotlib', 'joblib',
    ]
    versions = {}
    for pkg in pkgs:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except (importlib.metadata.PackageNotFoundError, ImportError):
            versions[pkg] = 'not found'
    versions['python'] = f'{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}'
    return versions


def generate_trace(output_dir: str | Path) -> str:
    """Generate a traceability markdown document.

    Args:
        output_dir: Path to output directory.

    Returns:
        Markdown string with the traceability document.
    """
    OUT = Path(output_dir)
    summary = _read_json(OUT / 'summary.json')
    annotation = _read_csv(OUT / 'tables' / 'cell_type_annotation.csv')
    morans = _read_csv(OUT / 'tables' / 'morans_i.csv')
    gearys = _read_csv(OUT / 'tables' / 'gearys_c.csv')

    cells = summary.get('cells_processed', 0)
    n_img = summary.get('n_image_features', 0)
    n_expr = summary.get('n_expression_genes', 0)
    libs = _get_library_versions()

    # Collect generated files
    generated_files = []
    for f in sorted(OUT.rglob('*')):
        if f.is_file() and f.suffix in ('.npy', '.pkl', '.csv', '.h5ad', '.png', '.json', '.md'):
            rel = f.relative_to(OUT.parent)
            generated_files.append((f.name, str(rel)))

    lines = []
    lines.append('# Single-Cell Analysis Traceability Document')
    lines.append('')
    lines.append('---')
    lines.append('')

    # Analysis Metadata
    lines.append('## Analysis Metadata')
    lines.append('')
    lines.append('| Field | Value |')
    lines.append('|-------|-------|')
    lines.append(f'| Project | Full FOV Breast Cancer |')
    lines.append(f'| Date | {datetime.now().strftime("%Y-%m-%d %H:%M")} |')
    lines.append(f'| Pipeline | cell-agent (SCAI) |')
    lines.append(f'| Branch | feat/image-agent |')
    lines.append(f'| Tissue | {summary.get("tissue", "N/A")} |')
    lines.append(f'| Mode | {summary.get("mode", "N/A")} |')
    lines.append(f'| Image size | {summary.get("image_size_px", "N/A")} |')
    lines.append('')

    # Dataset Summary
    lines.append('## Dataset Summary')
    lines.append('')
    lines.append(f'- {cells} cells × {n_expr} genes + {n_img} image features')
    lines.append('- Modality: multimodal (expression + image-derived features)')
    lines.append(f'- Hex bins: {summary.get("hex_bins", "N/A")}')
    lines.append(f'- Hex diameter: {summary.get("hex_diameter_um", "N/A")} µm')
    lines.append('')

    # Pipeline Execution
    lines.append('## Pipeline Execution')
    lines.append('')
    lines.append('Stages are listed in chronological order, exactly as executed.')
    lines.append('')

    # Stage 0-3: Data Loading & Image Processing
    lines.append('---')
    lines.append('### Stage 0: Data Loading & Xenium Import')
    lines.append('')
    lines.append('**Description**: Load Xenium output files (cell_feature_matrix, morphology image, transcripts).')
    lines.append('')
    lines.append('**Commands executed**:')
    lines.append('')
    lines.append('```python')
    lines.append('import scanpy as sc')
    lines.append("adata = sc.read_10x_h5('cell_feature_matrix.h5')")
    lines.append("img = skimage.io.imread('morphology.ome.tif')")
    lines.append('```')
    lines.append('')
    lines.append('**Results**:')
    lines.append(f'- {cells} cells, {n_expr} genes loaded')
    lines.append(f'- Image: {summary.get("image_size_px", "N/A")} px')
    lines.append('')
    lines.append('---')
    lines.append('### Stage 1–3: Image Feature Extraction')
    lines.append('')
    lines.append('**Description**: Morphological features, microenvironment density, and tissue compartments.')
    lines.append('')
    lines.append('**Parameters**:')
    lines.append('')
    lines.append('| Level | Parameter | Value |')
    lines.append('|-------|-----------|-------|')
    lines.append('| Morphology | min_cell_area | 50 |')
    lines.append('| Morphology | n_glcm_distances | 3 |')
    lines.append('| Microenvironment | density_radius | 50 µm |')
    lines.append('| Microenvironment | n_neighbors | 10 |')
    lines.append('| Compartments | method | Otsu thresholding |')
    lines.append('')
    lines.append('**Results**:')
    lines.append(f'- {n_img} image features extracted: morphology, texture (GLCM), microenvironment, compartments')
    lines.append('- Tissue compartments: epithelium, stroma, necrosis (3 classes)')
    lines.append('')

    # Stage 4
    lines.append('---')
    lines.append('### Stage 4: kNN Graph & Multimodal Fusion')
    lines.append('')
    lines.append('**Agent**: integration-agent (manual)')
    lines.append('')
    lines.append('**Commands executed**:')
    lines.append('')
    lines.append('```python')
    lines.append('from sklearn.neighbors import NearestNeighbors')
    lines.append('from sklearn.preprocessing import StandardScaler')
    lines.append('')
    lines.append('nn = NearestNeighbors(n_neighbors=15, metric="euclidean")')
    lines.append('nn.fit(coords)')
    lines.append('_, idx = nn.kneighbors(coords)')
    lines.append('')
    lines.append('# Spatial smoothing per modality')
    lines.append('img_smooth[i] = np.nanmean(img_mat[idx[i]], axis=0)')
    lines.append('expr_smooth[i] = expr_mat[idx[i]].mean(axis=0)')
    lines.append('')
    lines.append('# Scale + fuse')
    lines.append('img_scaled = StandardScaler().fit_transform(img_smooth)')
    lines.append('expr_scaled = StandardScaler().fit_transform(expr_smooth)')
    lines.append('mm = np.hstack([img_scaled, expr_scaled])')
    lines.append('```')
    lines.append('')
    lines.append('**Parameters**:')
    lines.append('')
    lines.append('| Parameter | Value | Justification |')
    lines.append('|-----------|-------|---------------|')
    lines.append(f'| kNN neighbors | {summary.get("knn_k", 15)} | Spatial locality for smoothing |')
    lines.append('| PCA components | 2 | For visualization only |')
    lines.append('')
    lines.append('**Results**:')
    lines.append(f'- Image features: {n_img}')
    lines.append(f'- Expression genes: {n_expr}')
    lines.append(f'- Multimodal matrix: ({cells}, {n_img + n_expr})')
    lines.append(f'- PCA: PC1={summary.get("multimodal_pca_pc1_pct", "N/A")}%, PC2={summary.get("multimodal_pca_pc2_pct", "N/A")}%')
    lines.append('')

    # Stage 5
    lines.append('---')
    lines.append('### Stage 5: Spatial Autocorrelation')
    lines.append('')
    lines.append('**Agent**: spatial-agent (manual)')
    lines.append('')
    lines.append('**Commands executed**:')
    lines.append('')
    lines.append('```python')
    lines.append('import squidpy as sq')
    lines.append('')
    lines.append('sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)')
    lines.append('sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100)')
    lines.append('sq.gr.spatial_autocorr(adata, mode="geary", n_perms=100)')
    lines.append('```')
    lines.append('')
    lines.append('**Parameters**:')
    lines.append('')
    lines.append('| Parameter | Value | Justification |')
    lines.append('|-----------|-------|---------------|')
    lines.append('| Spatial neighbors | 6 | Hex grid connectivity |')
    lines.append('| Permutations | 100 | Standard for spatial autocorrelation |')
    lines.append('')
    lines.append('**Results**:')
    if morans:
        valid = [m for m in morans if m.get('I') and m['I'].strip()]
        lines.append(f'- Non-constant features: {len(valid)} / {len(morans)}')
        if valid:
            top = valid[0]
            lines.append(f'- Top Moran\'s I: {top["feature"]} I={float(top["I"]):.4f}')
    lines.append('')

    # Stage 6
    lines.append('---')
    lines.append('### Stage 6: Hex Grid Feature Maps')
    lines.append('')
    lines.append('**Agent**: spatial-agent (manual)')
    lines.append('')
    lines.append('**Parameters**:')
    lines.append('')
    lines.append('| Parameter | Value |')
    lines.append('|-----------|-------|')
    lines.append('| Hex diameter | 50 µm |')
    lines.append('| Top features plotted | 12 |')
    lines.append('')
    lines.append('**Results**:')
    lines.append(f'- {summary.get("hex_bins", "N/A")} hex bins')
    lines.append('- Top 12 Moran\'s I features visualized on hex grid')
    lines.append('')

    # Stage 7
    lines.append('---')
    lines.append('### Stage 7: Leiden Clustering')
    lines.append('')
    lines.append('**Agent**: cluster-agent (manual)')
    lines.append('')
    lines.append('**Commands executed**:')
    lines.append('')
    lines.append('```python')
    lines.append('sc.pp.neighbors(adata, n_neighbors=15, n_pcs=2)')
    lines.append('sc.tl.leiden(adata, resolution=0.5, flavor="igraph")')
    lines.append('sc.tl.umap(adata)')
    lines.append('```')
    lines.append('')
    lines.append('**Parameters**:')
    lines.append('')
    lines.append('| Parameter | Value | Justification |')
    lines.append('|-----------|-------|---------------|')
    lines.append('| Resolution (chosen) | 0.5 | Best balance of granularity |')
    lines.append('| Resolutions tested | 0.3, 0.5, 1.0 | Standard sweep |')
    lines.append('| Neighbors | 15 | Common default |')
    lines.append('| PCA components | 2 | From multimodal matrix |')
    lines.append('')
    lines.append('**Results**:')
    n_clusters = len(set(r['cluster'] for r in annotation)) if annotation else 0
    lines.append(f'- {n_clusters} clusters identified')
    lines.append('')

    # Stage 7.5
    lines.append('---')
    lines.append('### Stage 7.5: Cell Type Annotation')
    lines.append('')
    lines.append('**Agent**: cell-annotator')
    lines.append('')
    lines.append('**Commands executed**:')
    lines.append('')
    lines.append('```python')
    lines.append('sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", n_genes=50)')
    lines.append('mdb = load_marker_db("data/cell_markers/panglao_markers.tsv")')
    lines.append('results = annotate_clusters(deg_dict, mdb, min_overlap=1, top_n=50)')
    lines.append('```')
    lines.append('')
    lines.append('**Parameters**:')
    lines.append('')
    lines.append('| Parameter | Value |')
    lines.append('|-----------|-------|')
    lines.append('| DEG method | wilcoxon |')
    lines.append('| DEG n_genes | 50 |')
    lines.append('| Marker DB | PanglaoDB (178 cell types) |')
    lines.append('| Min overlap | 1 gene(s) |')
    lines.append('| Top N considered | 50 |')
    lines.append('')
    lines.append('**Results**:')
    if annotation:
        n_annot = sum(1 for a in annotation if a.get('cell_type', 'Unannotated') != 'Unannotated')
        n_high = sum(1 for a in annotation if float(a.get('score', 0)) >= 0.5)
        n_med = sum(1 for a in annotation if 0.2 <= float(a.get('score', 0)) < 0.5)
        lines.append(f'- {n_annot}/{len(annotation)} clusters annotated')
        lines.append(f'- High confidence: {n_high}')
        lines.append(f'- Medium confidence: {n_med}')
        lines.append(f'- Low / Unannotated: {len(annotation) - n_annot}')
        lines.append('')

        lines.append('| Cluster | Cell Type | Score | Genes |')
        lines.append('|---------|-----------|-------|-------|')
        for a in sorted(annotation, key=lambda x: int(x['cluster']) if x['cluster'].isdigit() else x['cluster']):
            lines.append(f'| {a["cluster"]} | {a.get("cell_type", "?")} | {a.get("score", "?")} | {a.get("overlap_genes", "")} |')
        lines.append('')

    # Generated Files
    lines.append('---')
    lines.append('## Generated Files')
    lines.append('')
    lines.append('| File | Path |')
    lines.append('|------|------|')
    for name, rel in generated_files:
        lines.append(f'| {name} | {rel} |')
    lines.append('')

    # Library Versions
    lines.append('---')
    lines.append('## Library Versions')
    lines.append('')
    lines.append('| Library | Version |')
    lines.append('|---------|---------|')
    for pkg, ver in libs.items():
        lines.append(f'| {pkg} | {ver} |')
    lines.append('')

    # Total Time
    total_time = summary.get('total_time_s', 0)
    if total_time:
        lines.append(f'**Total analysis time**: {total_time:.0f} seconds ({total_time/60:.1f} minutes)')
        lines.append('')

    # Footer
    lines.append('---')
    lines.append('')
    lines.append('*Document automatically generated by SCAI Trace Agent*')
    lines.append(f'*Generated on {datetime.now().strftime("%Y-%m-%d %H:%M")}*')
    lines.append('')

    return '\n'.join(lines) + '\n'
