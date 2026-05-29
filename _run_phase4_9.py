"""Run pipeline from Phase 4 through Phase 9 using saved Phase 0-3 output."""
import sys, time, importlib
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

sc.settings.verbosity = 0

OUT = Path('test_data_2/output_fullfov')
OUT.joinpath('plots').mkdir(parents=True, exist_ok=True)
OUT.joinpath('tables').mkdir(parents=True, exist_ok=True)

m = pd.read_csv(OUT / 'tables' / 'aligned_cell_data.csv')
print(f'Loaded: {len(m)} cells, {len(m.columns)} cols')

# Reconstruct compartment label from one-hot columns
cd = [c for c in m.columns
      if c.startswith('compartment_')
      and c not in ('compartment_boundary_distance',
                    'compartment_intensity_gradient',
                    'compartment_border')]
cl = m[cd].idxmax(axis=1).str.replace('compartment_', '')

# ═══════════════════════════════════════════════════════════════
# PHASE 4: Shared kNN Graph & Multimodal Fusion
# ═══════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PHASE 4: Shared kNN Graph & Multimodal Fusion')
print('=' * 60)
t0 = time.time()

N_NEIGHBORS = 15
coords = m[['cx_microns', 'cy_microns']].values
nn = NearestNeighbors(n_neighbors=N_NEIGHBORS, metric='euclidean')
nn.fit(coords)
_, idx = nn.kneighbors(coords)

# Image features
feat_img = [c for c in m.columns
            if c.startswith(('morph_', 'texture_', 'microenv_'))
            or (c.startswith('compartment_') and c != 'compartment_label')]
img_mat = m[feat_img].values.astype(float)
img_smooth = np.zeros_like(img_mat)
for i in range(len(m)):
    img_smooth[i] = np.nanmean(img_mat[idx[i]], axis=0)

# Expression features — exclude summary stats
feat_expr = [c for c in m.columns
             if c.startswith('expr_')
             and c not in ('expr_n_genes', 'expr_total_counts', 'expr_mean')]
expr_mat = m[feat_expr].values
expr_smooth = np.zeros_like(expr_mat)
for i in range(len(m)):
    expr_smooth[i] = expr_mat[idx[i]].mean(axis=0)

# Scale each modality independently
si = StandardScaler()
se = StandardScaler()
img_scaled = si.fit_transform(img_smooth)
expr_scaled = se.fit_transform(expr_smooth)

joblib.dump(si, OUT / 'scaler_img.pkl')
joblib.dump(se, OUT / 'scaler_expr.pkl')

# Fuse
mm = np.hstack([img_scaled, expr_scaled])
mm = np.nan_to_num(np.asarray(mm, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
pca = PCA(n_components=2)
pca_r = pca.fit_transform(mm)

print(f'  Image features: {img_smooth.shape[1]}')
print(f'  Expression genes: {expr_smooth.shape[1]}')
print(f'  Multimodal matrix: {mm.shape}')
print(f'  PCA: PC1={pca.explained_variance_ratio_[0]:.1%}, '
      f'PC2={pca.explained_variance_ratio_[1]:.1%}')

np.save(OUT / 'multimodal_matrix.npy', mm)
print(f'  Phase 4 done in {time.time()-t0:.1f}s')

# ═══════════════════════════════════════════════════════════════
# PHASE 5: Spatial Autocorrelation (Moran's I + Geary's C)
# ═══════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PHASE 5: Spatial Autocorrelation')
print('=' * 60)
t0 = time.time()

# Build var_names
gene_names = [c.replace('expr_', '') for c in feat_expr]
var_names = list(feat_img) + [f'expr_{g}' for g in gene_names]

adata = sc.AnnData(X=mm)
adata.var_names = var_names
print(f'  var_names assigned: {len(var_names)}')

adata.obs['compartment'] = cl.values
adata.obs['cell_id'] = m['cell_id'].values
adata.obsm['spatial'] = m[['cx_microns', 'cy_microns']].values
adata.obsm['X_pca'] = pca_r

# Spatial neighbors
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

# Moran's I
print('  Computing Moran\'s I...')
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100)
mr = adata.uns['moranI']
md = pd.DataFrame({
    'feature': mr['I'].index,
    'I': mr['I'].values,
    'pval_sim': mr['pval_sim'].values
}).sort_values('I', ascending=False)
md.to_csv(OUT / 'tables' / 'morans_i.csv', index=False)
top = md.iloc[0]
print(f'  Top Moran I: {top["feature"]} I={top["I"]:.4f} p={top["pval_sim"]:.4e}')
valid = md.dropna(subset=['I'])
print(f'  Non-constant features: {len(valid)} / {len(md)}')

# Geary's C
print('  Computing Geary\'s C...')
sq.gr.spatial_autocorr(adata, mode='geary', n_perms=100)
gd = pd.DataFrame({
    'feature': adata.uns['gearyC']['C'].index,
    'C': adata.uns['gearyC']['C'].values,
    'pval_sim': adata.uns['gearyC']['pval_sim'].values
}).sort_values('C')
gd.to_csv(OUT / 'tables' / 'gearys_c.csv', index=False)
print(f'  Phase 5 done in {time.time()-t0:.1f}s')

# ═══════════════════════════════════════════════════════════════
# PHASE 6: Hex Grid Feature Maps
# ═══════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PHASE 6: Hex Grid Feature Maps')
print('=' * 60)

ha = pd.read_csv(OUT / 'tables' / 'hex_grid_aggregated.csv')
top_feats = md.dropna(subset=['I']).head(12)['feature'].tolist()

fig, axes = plt.subplots(3, 4, figsize=(16, 12))
axes = axes.flatten()
for i, feat in enumerate(top_feats):
    if feat in ha.columns:
        sv = ha[feat].values
        valid = ~(np.isnan(sv) | np.isinf(sv))
        if valid.sum() > 0:
            lo, hi = np.percentile(sv[valid], [2, 98])
            axes[i].scatter(
                ha['cx_microns'].values[valid],
                ha['cy_microns'].values[valid],
                c=sv[valid], cmap='viridis', s=20,
                vmin=lo, vmax=hi
            )
            axes[i].set_title(feat[:30], fontsize=9)
    axes[i].set_xlabel('X (um)')
    axes[i].set_ylabel('Y (um)')
for j in range(len(top_feats), 12):
    axes[j].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / 'plots' / 'hex_grid_features.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print('  hex_grid_features.png saved')

# Save top features for report
md.to_csv(OUT / 'tables' / 'morans_i_top100.csv', index=False)

# ═══════════════════════════════════════════════════════════════
# PHASE 7: Leiden Clustering
# ═══════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PHASE 7: Leiden Clustering')
print('=' * 60)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(10, pca_r.shape[1]))
for res in [0.2, 0.3, 0.5, 1.0]:
    sc.tl.leiden(adata, resolution=res, flavor='igraph',
                 n_iterations=2, key_added=f'leiden_r{res}')
adata.obs['leiden'] = adata.obs['leiden_r0.3']
n_clusters = adata.obs['leiden'].nunique()
print(f'  Clusters (r=0.5): {n_clusters}')

# UMAP & spatial
sc.tl.pca(adata, n_comps=2)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=2)
sc.tl.umap(adata)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='Leiden (UMAP)')
sc.pl.spatial(adata, color='leiden', spot_size=40, ax=axes[1],
              show=False, title='Leiden (Spatial)')
fig.tight_layout()
fig.savefig(OUT / 'plots' / 'leiden_clusters.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print('  leiden_clusters.png saved')

pd.DataFrame({
    'cell_id': adata.obs['cell_id'],
    'leiden': adata.obs['leiden'],
    'compartment': adata.obs['compartment']
}).to_csv(OUT / 'tables' / 'leiden_clusters.csv', index=False)

# ═══════════════════════════════════════════════════════════════
# PHASE 7.5: Marker Gene Detection + Cell Type Annotation
# ═══════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PHASE 7.5: Cell Type Annotation')
print('=' * 60)

sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', n_genes=50)
print('  rank_genes_groups OK')

# Load annotation module (importlib due to hyphen in path)
am = importlib.import_module('skills.cell-annotator.annotate')
mdb = am.load_marker_db('data/cell_markers/panglao_markers.tsv', tissue_context='breast')
n_entries = sum(len(v) for v in mdb.values())
print(f'  Marker DB: {n_entries} entries, {len(mdb)} types')

# Extract DEGs per cluster (scanpy >=1.12 uses recarray)
rg = adata.uns['rank_genes_groups']
cluster_ids = rg['names'].dtype.names
deg_dict = {}
for g in cluster_ids:
    raw = rg['names'][g].tolist()
    clean = [x.replace('expr_', '') for x in raw if x.startswith('expr_')]
    deg_dict[g] = clean
    print(f'  Cluster {g}: {len(clean)} DEGs')

# Annotate
results = am.annotate_clusters(deg_dict, mdb, min_overlap=1, top_n=50)
print('\nCell Type Annotation Results:')
print(am.format_annotation_table(results))

n_annot = sum(1 for r in results.values() if r['cell_type'] != 'Unannotated')
n_high = sum(1 for r in results.values() if r['score'] >= 0.5)
n_med = sum(1 for r in results.values() if 0.2 <= r['score'] < 0.5)
print(f'\nSummary: {n_annot}/{len(results)} annotated')
print(f'  High confidence (score >= 0.5):  {n_high}')
print(f'  Medium confidence (0.2-0.5):     {n_med}')
print(f'  Low / Unannotated:               {len(results) - n_annot}')

# Save annotation table
adf = pd.DataFrame([{
    'cluster': k,
    'cell_type': v['cell_type'],
    'score': v['score'],
    'confidence': am.confidence_label(v['score']),
    'overlap_genes': ', '.join(v['matching_markers']),
    'marker_genes': ', '.join(v['matching_markers'][:5])
} for k, v in results.items()])
adf.to_csv(OUT / 'tables' / 'cell_type_annotation.csv', index=False)

# Map to adata
ct_map = {k: v['cell_type'] for k, v in results.items()}
adata.obs['cell_type'] = adata.obs['leiden'].map(ct_map)

# Plot
fig, axes = plt.subplots(1, 2, figsize=(16, 5))
sc.pl.umap(adata, color='cell_type', ax=axes[0], show=False,
           title='Cell Types (UMAP)')
sc.pl.spatial(adata, color='cell_type', spot_size=40, ax=axes[1],
              show=False, title='Cell Types (Spatial)')
fig.tight_layout()
fig.savefig(OUT / 'plots' / 'cell_type_annotation.png', dpi=150,
            bbox_inches='tight')
plt.close(fig)
print('  cell_type_annotation.png saved')

# Save full adata
adata.write(OUT / 'adata_annotated.h5ad')
print('  adata_annotated.h5ad saved')

# ═══════════════════════════════════════════════════════════════
# PHASE 8: Executive Report
# ═══════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PHASE 8: Executive Report')
print('=' * 60)

rm = importlib.import_module('skills.report-agent.report')
report_md = rm.generate_report(OUT)
(OUT / 'report.md').write_text(report_md)
print('  report.md saved')

# ═══════════════════════════════════════════════════════════════
# PHASE 9: Traceability Document
# ═══════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PHASE 9: Traceability Document')
print('=' * 60)

tm = importlib.import_module('skills.trace-agent.trace')
trace_md = tm.generate_trace(OUT)
(OUT / 'trace.md').write_text(trace_md)
print('  trace.md saved')

# ───────────────────────────────────────────────────────────────
print(f'\n{"=" * 60}')
print(f'PIPELINE COMPLETE')
print(f'{time.time()-t0:.0f}s total for Phases 4-9')
print(f'{"=" * 60}')
