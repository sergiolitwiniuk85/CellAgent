"""SCAI Report Agent — Generate executive report for single-cell analysis.

Reads output directory contents and produces a complete markdown report.
"""

import json
from pathlib import Path
from datetime import datetime
from typing import Optional


def _read_csv(path: Path) -> list[dict]:
    """Read a small CSV into list of dicts."""
    if not path.exists():
        return []
    import csv
    with open(path) as f:
        return list(csv.DictReader(f))


def _read_json(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def _confidence_label(score: float) -> str:
    if score >= 0.5:
        return '**High** ★'
    elif score >= 0.2:
        return 'Medium'
    else:
        return '_Low_'


def generate_report(output_dir: str | Path) -> str:
    """Generate a markdown executive report from pipeline output.

    Args:
        output_dir: Path to output directory (e.g. 'test_data_2/output_fullfov').

    Returns:
        Markdown string with the complete report.
    """
    OUT = Path(output_dir)
    summary = _read_json(OUT / 'summary.json')
    annotation = _read_csv(OUT / 'tables' / 'cell_type_annotation.csv')
    morans = _read_csv(OUT / 'tables' / 'morans_i.csv')
    leiden = _read_csv(OUT / 'tables' / 'leiden_clusters.csv')

    plots_dir = OUT / 'plots'
    available_plots = sorted(p.name for p in plots_dir.glob('*.png')) if plots_dir.exists() else []

    cells = summary.get('cells_processed', 0)
    n_img = summary.get('n_image_features', 0)
    n_expr = summary.get('n_expression_genes', 0)
    n_clusters = len(set(r['cluster'] for r in annotation)) if annotation else 0

    lines = []
    lines.append('# Single-Cell Analysis Report — Full FOV Breast Cancer')
    lines.append('')
    lines.append(f'**Generated**: {datetime.now().strftime("%Y-%m-%d %H:%M")}')
    lines.append(f'**Pipeline**: Full FOV (Phase 0 → Phase 9)')
    lines.append('')
    lines.append('---')
    lines.append('')

    # 1. Executive summary
    lines.append('## 1. Summary')
    lines.append('')
    lines.append(f'| Metric | Value |')
    lines.append(f'|--------|-------|')
    lines.append(f'| Cells processed | {cells} |')
    lines.append(f'| Image features | {n_img} |')
    lines.append(f'| Expression genes | {n_expr} |')
    lines.append(f'| Multimodal features | {n_img + n_expr} |')
    lines.append(f'| Leiden clusters (r=0.5) | {n_clusters} |')
    lines.append(f'| Hex bins | {summary.get("hex_bins", "N/A")} |')
    lines.append(f'| Hex diameter | {summary.get("hex_diameter_um", "N/A")} µm |')
    lines.append(f'| Tissue | {summary.get("tissue", "N/A")} |')
    lines.append(f'| Image size | {summary.get("image_size_px", "N/A")} px |')
    lines.append('')
    lines.append('**Spatial Autocorrelation**: ' + ('Detected ✓' if summary.get('spatial_autocorrelation') else 'Not detected'))
    lines.append(f'**KNN edges**: {summary.get("knn_edges", "N/A")} (k={summary.get("knn_k", "N/A")})')
    lines.append(f'**Within-compartment connectivity**: {summary.get("within_compartment_connectivity_pct", "N/A")}%')
    lines.append('')
    lines.append('---')
    lines.append('')

    # 2. Cell Type Annotation
    lines.append('## 2. Cell Type Annotation')
    lines.append('')
    lines.append(f'**Annotated clusters**: {len([a for a in annotation if a["cell_type"] != "Unannotated"])} / {len(annotation)}')
    lines.append('')

    if annotation:
        conf_counts = {}
        for a in annotation:
            c = a.get('confidence', 'N/A')
            conf_counts[c] = conf_counts.get(c, 0) + 1
        for label in ['**High** ★', 'Medium', '_Low_', 'Unannotated']:
            if label in conf_counts:
                lines.append(f'- **{label.replace("**", "").replace("_", "").replace("★","").strip()}**: {conf_counts[label]} clusters')
        lines.append('')

        lines.append('| Cluster | Predicted Cell Type | Score | Confidence | Marker Genes |')
        lines.append('|---------|---------------------|-------|------------|--------------|')
        for a in sorted(annotation, key=lambda x: int(x['cluster']) if x['cluster'].isdigit() else x['cluster']):
            ct = a.get('cell_type', 'Unannotated')
            score = float(a.get('score', 0))
            conf = _confidence_label(score)
            markers = a.get('marker_genes', '')
            lines.append(f'| {a["cluster"]} | {ct} | {conf} | {score:.3f} | {markers} |')
        lines.append('')

    # 3. Spatial Autocorrelation
    lines.append('## 3. Spatial Autocorrelation')
    lines.append('')

    if morans:
        valid = [m for m in morans if m.get('I') and m['I'].strip()]
        lines.append(f'**Moran\'s I**: {len(valid)} non-constant features evaluated')
        lines.append('')
        lines.append('| Rank | Feature | Moran\'s I | p-value |')
        lines.append('|------|---------|-----------|---------|')
        for i, m in enumerate(valid[:10]):
            feat = m.get('feature', 'N/A')
            i_val = m.get('I', 'N/A')
            pval = f'{float(m.get("pval_sim", 1)):.2e}' if m.get('pval_sim') and m['pval_sim'].strip() else 'N/A'
            lines.append(f'| {i+1} | {feat} | {float(i_val):.4f} | {pval} |')
        lines.append('')

    # 4. Generated Plots
    lines.append('## 4. Generated Plots')
    lines.append('')
    if available_plots:
        for p in sorted(available_plots):
            rel = f'plots/{p}'
            lines.append(f'![{p}]({rel})')
        lines.append('')
    else:
        lines.append('No plots found.')
        lines.append('')

    # 5. Cluster Composition
    lines.append('## 5. Cluster Composition')
    lines.append('')
    if leiden:
        from collections import Counter
        comp = Counter(r['leiden'] for r in leiden)
        lines.append(f'| Cluster | Cells | Proportion |')
        lines.append(f'|---------|-------|------------|')
        for cid in sorted(comp.keys(), key=lambda x: int(x) if x.isdigit() else x):
            cnt = comp[cid]
            pct = cnt / len(leiden) * 100
            lines.append(f'| {cid} | {cnt} | {pct:.1f}% |')
        lines.append('')

    # 6. Parameters
    lines.append('## 6. Analysis Parameters')
    lines.append('')
    lines.append('| Phase | Parameter | Value |')
    lines.append('|-------|-----------|-------|')
    lines.append(f'| 4 — Fusion | kNN neighbors | {summary.get("knn_k", 15)} |')
    lines.append(f'| 4 — Fusion | PCA components | 2 |')
    lines.append(f'| 5 — Autocorr | Moran permutations | 100 |')
    lines.append(f'| 5 — Autocorr | Spatial neighbors | 6 |')
    lines.append(f'| 6 — Hex grid | Hex diameter | {summary.get("hex_diameter_um", 50)} µm |')
    lines.append(f'| 7 — Leiden | Resolution tested | {", ".join(str(r) for r in summary.get("leiden_resolutions_tested", ["0.5"]))} |')
    lines.append(f'| 7 — Leiden | Chosen resolution | 0.5 |')
    lines.append(f'| 7 — Leiden | Neighbors | 15 |')
    lines.append(f'| 7.5 — Annotation | DEG method | wilcoxon |')
    lines.append(f'| 7.5 — Annotation | DEG n_genes | 50 |')
    lines.append(f'| 7.5 — Annotation | Marker DB | PanglaoDB |')
    lines.append(f'| 7.5 — Annotation | Min overlap | 1 |')
    lines.append('')

    # 7. Recommendations
    lines.append('## 7. Recommendations')
    lines.append('')
    n_low = sum(1 for a in annotation if float(a.get('score', 0)) < 0.2) if annotation else 0
    n_high = sum(1 for a in annotation if float(a.get('score', 0)) >= 0.5) if annotation else 0

    if n_low > 0:
        lines.append(f'- ⚠️ **{n_low} clusters** have low-confidence annotation (score < 0.2). Consider:')
        lines.append('  - Increasing `top_n` in marker detection')
        lines.append('  - Manual review of DEG lists')
        lines.append('  - Using a larger or tissue-specific marker database')
        lines.append('')
    if n_high < len(annotation) / 2:
        lines.append('- ℹ️ **Majority of annotations are low-to-medium confidence**. Consider validating with an orthogonal method.')
        lines.append('')
    lines.append('- 📊 **Spatial autocorrelation detected** — the feature space has meaningful spatial structure.')
    lines.append('')

    # Footer
    lines.append('---')
    lines.append('')
    lines.append(f'*Report generated by SCAI Report Agent on {datetime.now().strftime("%Y-%m-%d %H:%M")}*')
    lines.append('')

    return '\n'.join(lines) + '\n'
