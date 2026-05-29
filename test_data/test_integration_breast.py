#!/usr/bin/env python3
"""
Full Integration Test: Image Features + RNA Expression (Breast Cancer)
========================================================================

Pipeline:
  Phase 1 — Extract image features (Level 1-3) for a tissue crop
  Phase 2 — Load RNA expression matrix, align by cell_id
  Phase 3 — Hexagonal grid: aggregate both modalities per spatial bin
  Phase 4 — Shared kNN: propagate features across neighborhood graph
  Phase 5 — Tissue compartments as spatial scaffold
  Phase 6 — Validation + PDF report generation

Data: Xenium Human Breast 2 FOV (test dataset)
"""

import numpy as np
import pandas as pd
import tifffile
import h5py
import skimage.measure
import skimage.feature
import skimage.filters
import skimage.morphology
import scipy.spatial
import scipy.ndimage
import scipy.sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from pathlib import Path
import joblib
from skimage.draw import polygon as sk_polygon
import warnings
import time
import json
import math


warnings.filterwarnings("ignore", category=UserWarning)

BASE = Path("/home/sergiolitwiniuk/gentle-projects/scai/test_data_2")
OUT = BASE / "output_integration"
OUT.mkdir(exist_ok=True)
(OUT / "plots").mkdir(exist_ok=True)
(OUT / "tables").mkdir(exist_ok=True)

MICRONS_TO_PX = 1.0 / 0.2125  # Xenium resolution

# ────────────────────────────────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────────────────────────────────
CROP_Y0, CROP_X0 = 2000, 1000     # Top-left of crop (pixels) — densest region
CROP_SIZE = 1500                    # Crop size (pixels) — larger for breast density
HEX_DIAMETER_MICRONS = 50.0         # Hex bin diameter in microns
HEX_DIAMETER_PX = HEX_DIAMETER_MICRONS * MICRONS_TO_PX
N_NEIGHBORS_KNN = 10                # k for shared kNN
MIN_CELL_AREA = 16                  # px² — skip tiny artifacts
DENSITY_RADIUS_PX = 50              # px — microenvironment radius
N_NND = 5                           # Number of nearest neighbors for NND

print("=" * 70)
print("INTEGRATION TEST: Image Features + RNA Expression (Breast Cancer)")
print("=" * 70)
t_start = time.time()

# ════════════════════════════════════════════════════════════════════════════
# PHASE 0: LOAD DATA
# ════════════════════════════════════════════════════════════════════════════

print("\n" + "─" * 70)
print("PHASE 0: Load Data")
print("─" * 70)

# --- 0a. Image ---
print("\nLoading morphology image...")
img = tifffile.imread(str(BASE / "morphology.ome.tif"))
print(f"  Raw image shape: {img.shape}, dtype={img.dtype}")
# Breast has 11 Z-planes — do maximum intensity projection
if img.ndim == 3:
    image_full = img.max(axis=0)  # MIP over Z
    print(f"  MIP over {img.shape[0]} Z-planes -> {image_full.shape}")
elif img.ndim == 2:
    image_full = img
print(f"  Full image: {image_full.shape}, dtype={image_full.dtype}")

# Crop
y1, x1 = CROP_Y0, CROP_X0
y2, x2 = y1 + CROP_SIZE, x1 + CROP_SIZE
image = image_full[y1:y2, x1:x2]
print(f"  Crop: ({CROP_SIZE}×{CROP_SIZE}) at ({y1},{x1})")
print(f"  Crop range: [{image.min()}, {image.max()}]")

# --- 0b. Cell boundaries ---
print("\nLoading nucleus boundaries...")
nuclei_df = pd.read_parquet(BASE / "nucleus_boundaries.parquet")
print(f"  Total cells: {nuclei_df['cell_id'].nunique()}")

# --- 0c. Cell metadata ---
print("\nLoading cell metadata...")
cells_df = pd.read_csv(BASE / "cells.csv.gz", compression="gzip")
print(f"  Cells: {len(cells_df)}")

# --- 0d. Expression matrix ---
print("\nLoading expression matrix...")
with h5py.File(str(BASE / "cell_feature_matrix.h5"), "r") as f:
    matrix = f["matrix"]
    barcodes_h5 = [b.decode() for b in matrix["barcodes"][:]]
    gene_names = [g.decode() for g in matrix["features"]["name"][:]]
    data = matrix["data"][:]
    indices = matrix["indices"][:]
    indptr = matrix["indptr"][:]
    n_genes, n_cells = matrix["shape"][:]

print(f"  Genes: {n_genes}, Cells: {n_cells}")
print(f"  Non-zero entries: {len(data)}")
print(f"  First 5 genes: {gene_names[:5]}")
print(f"  First 5 barcodes (h5): {barcodes_h5[:5]}")
print(f"  First 5 cell_ids (cells.csv): {cells_df['cell_id'].head(5).tolist()}")

# Verify ordering matches
assert barcodes_h5 == cells_df["cell_id"].tolist(), "Cell ordering mismatch between h5 and cells.csv!"
print("  ✅ Cell ordering: h5 ↔ cells.csv match")

# h5 stores in CSC format (10X standard): genes × cells
# Use CSC directly then transpose to cells × genes
expr_csc = scipy.sparse.csc_matrix((data, indices, indptr),
                                   shape=(int(n_genes), int(n_cells)))
expr = expr_csc.T.toarray()  # cells × genes
print(f"  Expression matrix: {expr.shape} (cells × genes)")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 1: EXTRACT IMAGE FEATURES
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("PHASE 1: Image Feature Extraction")
print("─" * 70)
t1 = time.time()

# --- 1a. Filter cells in crop & rasterize ---
print("\nFiltering cells in crop region...")
xs_all = (nuclei_df["vertex_x"].values * MICRONS_TO_PX).astype(int)
ys_all = (nuclei_df["vertex_y"].values * MICRONS_TO_PX).astype(int)

in_crop = ((xs_all >= x1) & (xs_all < x2) &
           (ys_all >= y1) & (ys_all < y2))
crop_cell_ids = nuclei_df["cell_id"].iloc[np.where(in_crop)[0]].unique()
print(f"  Cells in crop: {len(crop_cell_ids)}")

if len(crop_cell_ids) > 3000:
    print(f"  ⚠ Large crop ({len(crop_cell_ids)} cells), may be slow.")

nuclei_crop = nuclei_df[nuclei_df["cell_id"].isin(crop_cell_ids)].copy()

print("Rasterizing cell masks...")
mask = np.zeros(image.shape, dtype=np.int32)

xs_crop = (nuclei_crop["vertex_x"].values * MICRONS_TO_PX).astype(int) - x1
ys_crop = (nuclei_crop["vertex_y"].values * MICRONS_TO_PX).astype(int) - y1
nuclei_crop = nuclei_crop.assign(px_x=xs_crop, px_y=ys_crop)

grouped = nuclei_crop.groupby("cell_id")
n_cells_crop = len(grouped)
cell_centroids_px = {}
cell_labels = {}  # cell_id -> label int
label = 1
skipped = 0

for cell_id, group in grouped:
    if len(group) < 3:
        skipped += 1
        continue
    xs = np.clip(group["px_x"].values, 0, CROP_SIZE - 1)
    ys = np.clip(group["px_y"].values, 0, CROP_SIZE - 1)
    rr, cc = sk_polygon(ys.astype(int), xs.astype(int))
    if len(rr) == 0:
        skipped += 1
        continue
    mask[rr, cc] = label
    cell_centroids_px[cell_id] = (
        float(ys.astype(int).mean()), float(xs.astype(int).mean()))
    cell_labels[cell_id] = label
    label += 1

n_cells_rasterized = len(cell_labels)
print(f"  Rasterized: {n_cells_rasterized} cells ({skipped} skipped)")
print(f"  Mask unique labels: {np.unique(mask).size - 1}")

label_to_cell_id = {v: k for k, v in cell_labels.items()}

# --- 1b. Level 1: Morphology ---
print("\nLevel 1: Morphology + GLCM...")
image_8bit = ((image - image.min()) /
              (image.max() - image.min()) * 255).astype(np.uint8)

props = skimage.measure.regionprops_table(
    mask, intensity_image=image,
    properties=('label', 'area', 'eccentricity', 'solidity',
                'perimeter', 'extent', 'equivalent_diameter', 'centroid')
)

morph_df = pd.DataFrame(props)
morph_df.rename(columns={
    'label': 'cell_label',
    'area': 'morph_area',
    'eccentricity': 'morph_eccentricity',
    'solidity': 'morph_solidity',
    'perimeter': 'morph_perimeter',
    'extent': 'morph_extent',
    'equivalent_diameter': 'morph_equiv_diameter',
}, inplace=True)

# Map labels back to cell_ids
morph_df["cell_id"] = morph_df["cell_label"].map(label_to_cell_id)

# Per-cell GLCM
print("  GLCM textures...")
glcm_data = {p: [] for p in
             ['contrast', 'dissimilarity', 'homogeneity', 'energy', 'correlation']}
for _, row in morph_df.iterrows():
    lbl = int(row['cell_label'])
    cell_mask = mask == lbl
    if cell_mask.sum() < MIN_CELL_AREA:
        for p in glcm_data:
            glcm_data[p].append(np.nan)
        continue
    ys, xs = np.where(cell_mask)
    y1_, x1_, y2_, x2_ = ys.min(), xs.min(), ys.max() + 1, xs.max() + 1
    cell_intensity = image_8bit[y1_:y2_, x1_:x2_]
    glcm = skimage.feature.graycomatrix(
        cell_intensity, distances=[1],
        angles=[0, np.pi/4, np.pi/2, 3*np.pi/4],
        symmetric=True, normed=True
    )
    for p in glcm_data:
        glcm_data[p].append(skimage.feature.graycoprops(glcm, p).mean())

for p, vals in glcm_data.items():
    morph_df[f"texture_{p}"] = vals

print(f"  → {len(morph_df)} cells, morphology features: "
      f"{[c for c in morph_df.columns if c.startswith(('morph_', 'texture_'))][:6]}...")

# --- 1c. Level 2: Microenvironment ---
print("\nLevel 2: Microenvironment...")
centroids = morph_df[['centroid-0', 'centroid-1']].values  # (y, x) from regionprops
centroids_xy = centroids[:, [1, 0]]  # (x, y)

tree = scipy.spatial.cKDTree(centroids_xy)

# Cell density
density = np.array([len(tree.query_ball_point(c, r=DENSITY_RADIUS_PX)) - 1
                    for c in centroids_xy])
morph_df["microenv_cell_density"] = density

# Edge distance
dist_from_boundary = scipy.ndimage.distance_transform_edt(mask > 0)
edge_dist = np.array([
    dist_from_boundary[int(row['centroid-0']), int(row['centroid-1'])]
    for _, row in morph_df.iterrows()
])
morph_df["microenv_edge_distance"] = edge_dist

# NND
nndist, _ = tree.query(centroids_xy, k=N_NND + 1)
for i in range(1, N_NND + 1):
    morph_df[f"microenv_nndist_{i}"] = nndist[:, i]

print(f"  → density: [{density.min()}, {density.max()}], "
      f"edge: [{edge_dist.min():.1f}, {edge_dist.max():.1f}], "
      f"NND1 mean: {nndist[:, 1].mean():.1f} px")

# --- 1d. Level 3: Tissue Compartments ---
print("\nLevel 3: Tissue Compartments...")
thresh = skimage.filters.threshold_otsu(image)
print(f"  Otsu threshold: {thresh:.1f}")

epithelium = skimage.morphology.closing(image > thresh,
                                        skimage.morphology.disk(15))
stroma = skimage.morphology.closing(image <= thresh,
                                     skimage.morphology.disk(15))
border = skimage.morphology.dilation(epithelium,
                                     skimage.morphology.disk(5)) & stroma
necrosis = image < (thresh * 0.3)

labels = []
comp_boundary_dist = []
comp_intensity_grad = []
for _, row in morph_df.iterrows():
    cy = int(min(row['centroid-0'], CROP_SIZE - 1))
    cx = int(min(row['centroid-1'], CROP_SIZE - 1))
    if necrosis[cy, cx]:
        labels.append("necrosis")
    elif epithelium[cy, cx] and not border[cy, cx]:
        labels.append("epithelium")
    elif border[cy, cx]:
        labels.append("border")
    elif stroma[cy, cx]:
        labels.append("stroma")
    else:
        labels.append("unclassified")
    ep_edge = skimage.morphology.dilation(epithelium,
                                          skimage.morphology.disk(1)) ^ epithelium
    if ep_edge.any():
        comp_boundary_dist.append(float(
            scipy.ndimage.distance_transform_edt(~ep_edge)[cy, cx]))
    else:
        comp_boundary_dist.append(np.nan)
    gy, gx = np.gradient(image.astype(float))
    comp_intensity_grad.append(np.sqrt(gy[cy, cx]**2 + gx[cy, cx]**2))

morph_df["compartment_label"] = labels
morph_df["compartment_boundary_distance"] = comp_boundary_dist
morph_df["compartment_intensity_gradient"] = comp_intensity_grad

comp_counts = morph_df["compartment_label"].value_counts()
print("  Compartment distribution:")
for comp, count in comp_counts.items():
    print(f"    {comp}: {count} ({count/len(morph_df)*100:.1f}%)")

t1_end = time.time()
print(f"\n⏱ Image features: {t1_end - t1:.1f}s")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 2: ALIGN EXPRESSION + IMAGE FEATURES
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("PHASE 2: Align Expression + Image Features by cell_id")
print("─" * 70)
t2 = time.time()

# Build feature matrix
feature_cols = [c for c in morph_df.columns
                if c.startswith(('morph_', 'texture_', 'microenv_', 'compartment_'))]
image_features = morph_df[['cell_id'] + feature_cols].copy()

# Convert compartment_label to one-hot
comp_dummies = pd.get_dummies(image_features['compartment_label'],
                               prefix='compartment')
image_features = pd.concat([image_features.drop(columns=['compartment_label']),
                            comp_dummies], axis=1)

feature_cols_final = [c for c in image_features.columns if c != 'cell_id']
print(f"  Image features: {len(feature_cols_final)} columns")

# Now align with expression: both go by cell_id
# Get expression for the cells in our crop
crop_idx = np.where(cells_df["cell_id"].isin(image_features["cell_id"]))[0]
expr_crop = expr[crop_idx]  # (n_cells_crop, n_genes)

# Expression summary stats per cell
expr_per_cell = pd.DataFrame({
    "cell_id": cells_df["cell_id"].iloc[crop_idx].values,
    "expr_n_genes": (expr_crop > 0).sum(axis=1),       # genes detected per cell
    "expr_total_counts": expr_crop.sum(axis=1),          # total UMI per cell
    "expr_mean": expr_crop.mean(axis=1),                 # mean expression
})

# Top variable genes — adaptive: all genes for targeted panels (<1000), top 2000 HVGs otherwise
n_genes_total = int(n_genes)
if n_genes_total >= 1000:
    n_top_var = 2000
    gene_var = np.var(expr_crop, axis=0)
    top_var_genes_idx = np.argsort(gene_var)[-n_top_var:][::-1]
    print(f"  Whole transcriptome ({n_genes_total} genes) → selecting top {n_top_var} HVGs")
else:
    n_top_var = n_genes_total
    top_var_genes_idx = np.arange(n_top_var)
    print(f"  Targeted panel ({n_genes_total} genes) → using ALL genes")
top_var_genes = [gene_names[i] for i in top_var_genes_idx]
print(f"  Top 5 variable genes: {top_var_genes[:5]}")

# Merge image features + expression summary
merged = image_features.merge(expr_per_cell, on="cell_id", how="inner")
print(f"  Merged cells: {len(merged)}")
print(f"  Total columns: {len(merged.columns)}")

# Full expression for top genes (for hex grid aggregation)
expr_top_genes = pd.DataFrame(
    expr_crop[:, top_var_genes_idx],
    columns=[f"expr_{g}" for g in top_var_genes]
)
expr_top_genes["cell_id"] = cells_df["cell_id"].iloc[crop_idx].values

# Merge everything
merged_full = merged.merge(expr_top_genes, on="cell_id", how="inner")
print(f"  After adding top genes: {len(merged_full.columns)} columns, "
      f"{len(merged_full)} cells")

# Centroids (in microns for the hexagonal grid)
cell_to_centroid_microns = {
    row["cell_id"]: (row["x_centroid"], row["y_centroid"])
    for _, row in cells_df.iterrows()
}
microns_df = pd.DataFrame([
    {"cell_id": cid,
     "cx_microns": cell_to_centroid_microns[cid][0],
     "cy_microns": cell_to_centroid_microns[cid][1]}
    for cid in merged_full["cell_id"] if cid in cell_to_centroid_microns
])
merged_full = merged_full.merge(microns_df, on="cell_id", how="inner")
print(f"  With spatial coordinates: {len(merged_full)} cells")

t2_end = time.time()
print(f"⏱ Alignment: {t2_end - t2:.1f}s")

# Save aligned data
merged_full.to_csv(OUT / "tables" / "aligned_cell_data.csv", index=False)
print(f"  Saved: aligned_cell_data.csv")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 3: HEXAGONAL GRID — Aggregate by Spatial Bin
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print(f"PHASE 3: Hexagonal Grid (diameter = {HEX_DIAMETER_MICRONS} µm)")
print("─" * 70)
t3 = time.time()

def hex_grid_bin(cx, cy, diameter):
    """Assign a point (cx, cy) to its hex bin center.
    Uses pointy-top hex grid with axial coordinates.
    """
    # Hex width = diameter, height = diameter * sqrt(3)/2
    # Axial coordinates
    q = (np.sqrt(3) / 3 * cx - 1.0 / 3 * cy) / (diameter / 2)
    r = (2.0 / 3 * cy) / (diameter / 2)
    # Round to nearest hex center
    q_round = np.round(q)
    r_round = np.round(r)
    s_round = np.round(-q - r)
    # Correct rounding
    q_diff = np.abs(q - q_round)
    r_diff = np.abs(r - r_round)
    s_diff = np.abs(-q - r - s_round)
    if q_diff > r_diff and q_diff > s_diff:
        q_round = -r_round - s_round
    elif r_diff > s_diff:
        r_round = -q_round - s_round
    return q_round, r_round

# Assign each cell to a hex bin
cx_all = merged_full["cx_microns"].values
cy_all = merged_full["cy_microns"].values
hex_q = np.zeros(len(merged_full))
hex_r = np.zeros(len(merged_full))
for i in range(len(merged_full)):
    q, r = hex_grid_bin(cx_all[i], cy_all[i], HEX_DIAMETER_MICRONS)
    hex_q[i] = q
    hex_r[i] = r

merged_full["hex_q"] = hex_q
merged_full["hex_r"] = hex_r
merged_full["hex_bin"] = [f"hex_{int(q)}_{int(r)}" for q, r in zip(hex_q, hex_r)]

n_bins = merged_full["hex_bin"].nunique()
print(f"  Cells in grid: {len(merged_full)}")
print(f"  Unique hex bins: {n_bins}")
print(f"  Cells per bin: mean={merged_full.groupby('hex_bin').size().mean():.1f}, "
      f"max={merged_full.groupby('hex_bin').size().max()}")

# --- Aggregate per hex bin ---
# Numeric columns to aggregate
image_feat_to_agg = [c for c in merged_full.columns if c.startswith('morph_') or
                     c.startswith('texture_') or c.startswith('microenv_') or
                     c.startswith('compartment_')]
expr_agg_cols = [c for c in merged_full.columns if c.startswith('expr_')]
coord_cols = ["cx_microns", "cy_microns"]

agg_dict = {}
for col in image_feat_to_agg + expr_agg_cols:
    agg_dict[col] = ["mean", "std"]
agg_dict["cx_microns"] = "mean"
agg_dict["cy_microns"] = "mean"

hex_agg = merged_full.groupby("hex_bin").agg(agg_dict).reset_index()
hex_agg.columns = ['_'.join(c).strip('_') if c[1] else c[0]
                   for c in hex_agg.columns.tolist()]

# Rename properly
new_cols = []
for c in hex_agg.columns:
    if c == 'hex_bin':
        new_cols.append(c)
    elif c.endswith('_mean'):
        new_cols.append(c.replace('_mean', ''))
    elif c.endswith('_std'):
        new_cols.append(c)
    else:
        new_cols.append(c)
hex_agg.columns = new_cols

hex_agg["n_cells"] = merged_full.groupby("hex_bin").size().values
print(f"  Hex grid aggregated: {len(hex_agg)} bins × {len(hex_agg.columns)} columns")

# Hex bin centroids in microns
hex_bin_centers = {}
for q_val, r_val in zip(merged_full.groupby("hex_bin")["hex_q"].first(),
                         merged_full.groupby("hex_bin")["hex_r"].first()):
    # Axial to pixel conversion
    x_center_microns = (HEX_DIAMETER_MICRONS / 2) * (np.sqrt(3) * q_val + np.sqrt(3)/2 * r_val)
    y_center_microns = (HEX_DIAMETER_MICRONS / 2) * (3.0 / 2 * r_val)
    hex_bin_centers[f"hex_{int(q_val)}_{int(r_val)}"] = (x_center_microns, y_center_microns)

hex_agg["cx_microns"] = hex_agg["hex_bin"].map(
    lambda b: hex_bin_centers.get(b, (np.nan, np.nan))[0])
hex_agg["cy_microns"] = hex_agg["hex_bin"].map(
    lambda b: hex_bin_centers.get(b, (np.nan, np.nan))[1])

hex_agg.to_csv(OUT / "tables" / "hex_grid_aggregated.csv", index=False)
print(f"  Saved: hex_grid_aggregated.csv")

t3_end = time.time()
print(f"⏱ Hex grid: {t3_end - t3:.1f}s")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 4: SHARED kNN — Multimodal Propagation
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print(f"PHASE 4: Shared kNN Graph (k={N_NEIGHBORS_KNN})")
print("─" * 70)
t4 = time.time()

# Use ALL spatial coordinates (image centroids in microns)
coords_knn = merged_full[["cx_microns", "cy_microns"]].values

# Build shared kNN graph
nn = NearestNeighbors(n_neighbors=N_NEIGHBORS_KNN, metric="euclidean")
nn.fit(coords_knn)
distances, indices = nn.kneighbors(coords_knn)

# --- Propagate image features through the graph ---
# For each cell, smooth its image features by averaging over its kNN neighborhood
image_feat_to_agg_final = [c for c in merged_full.columns
                           if c.startswith(('morph_', 'texture_', 'microenv_'))
                           or (c.startswith('compartment_') and c != 'compartment_label')]
image_feat_matrix = merged_full[image_feat_to_agg_final].values.astype(float)
print(f"  Image feature matrix: {image_feat_matrix.shape}, NaN: {np.isnan(image_feat_matrix).sum()}")

# Propagate: for each cell, average the features of its neighbors
image_feat_smoothed = np.zeros_like(image_feat_matrix)
for i in range(len(merged_full)):
    neighbors = indices[i]
    image_feat_smoothed[i] = np.nanmean(image_feat_matrix[neighbors], axis=0)

# --- Propagate expression through the same graph ---
expr_matrix = merged_full[[c for c in merged_full.columns
                          if c.startswith('expr_') and c != 'expr_n_genes'
                          and c != 'expr_total_counts']].values

expr_smoothed = np.zeros_like(expr_matrix)
for i in range(len(merged_full)):
    neighbors = indices[i]
    expr_smoothed[i] = expr_matrix[neighbors].mean(axis=0)

# --- Scale each modality independently before fusion ---
scaler_img = StandardScaler()
img_scaled = scaler_img.fit_transform(image_feat_smoothed)
scaler_expr = StandardScaler()
expr_scaled = scaler_expr.fit_transform(expr_smoothed)

joblib.dump(scaler_img, OUT / "scaler_img.pkl")
joblib.dump(scaler_expr, OUT / "scaler_expr.pkl")
print(f"  Saved scalers: scaler_img.pkl, scaler_expr.pkl")

multimodal_matrix = np.hstack([img_scaled, expr_scaled])

print(f"  Graph edges: {len(merged_full)} × {N_NEIGHBORS_KNN} = "
      f"{len(merged_full) * N_NEIGHBORS_KNN}")
print(f"  Multimodal matrix shape: {multimodal_matrix.shape}")
print(f"    Image features: {image_feat_smoothed.shape[1]}")
print(f"    Expression (smoothed): {expr_smoothed.shape[1]}")

# PCA on multimodal matrix — ensure clean numeric
print(f"  multimodal_matrix dtype: {multimodal_matrix.dtype}")
print(f"  multimodal_matrix has NaN: {np.isnan(multimodal_matrix).any()}")
print(f"  multimodal_matrix has Inf: {np.isinf(multimodal_matrix).any()}")
multimodal_matrix = np.asarray(multimodal_matrix, dtype=float)
multimodal_matrix_clean = np.nan_to_num(multimodal_matrix, nan=0.0, posinf=0.0, neginf=0.0)
print(f"  Clean shape: {multimodal_matrix_clean.shape}, NaN after: {np.isnan(multimodal_matrix_clean).sum()}")
pca_mm = PCA(n_components=2)
pca_mm_result = pca_mm.fit_transform(multimodal_matrix_clean)
print(f"  PCA variance explained: "
      f"{pca_mm.explained_variance_ratio_[0]*100:.1f}% + "
      f"{pca_mm.explained_variance_ratio_[1]*100:.1f}% = "
      f"{pca_mm.explained_variance_ratio_[:2].sum()*100:.1f}%")

t4_end = time.time()
print(f"⏱ Shared kNN: {t4_end - t4:.1f}s")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 5: Compartments as Spatial Scaffold
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("PHASE 5: Compartments as Spatial Scaffold")
print("─" * 70)
t5 = time.time()

# Use compartment labels to stratify the kNN
# For each compartment, build a separate kNN and compare within vs across

# Use compartment labels from morph_df (before dummy encoding)
comp_label_map = morph_df.set_index("cell_id")["compartment_label"].to_dict()
merged_full["compartment_label"] = merged_full["cell_id"].map(comp_label_map)

# Per-compartment feature profiles
compartment_profiles = merged_full.groupby("compartment_label")[
    image_feat_to_agg_final + [c for c in merged_full.columns
                         if c.startswith('expr_') and c != 'expr_n_genes'
                         and c != 'expr_total_counts']
].mean()

print("  Per-compartment feature profiles:")
for comp in compartment_profiles.index:
    top_expr = compartment_profiles.loc[comp].filter(like='expr_').sort_values(ascending=False)
    top_genes = top_expr.index[:3].tolist()
    print(f"    {comp}: {len(merged_full[merged_full['compartment_label'] == comp])} cells")

# --- Cross-compartment connectivity ---
# What fraction of kNN edges connect cells in the same compartment?
same_comp = 0
total_edges = 0
comp_labels = merged_full["compartment_label"].values
for i in range(len(merged_full)):
    for j in indices[i]:
        if i == j:
            continue
        total_edges += 1
        if comp_labels[i] == comp_labels[j]:
            same_comp += 1

print(f"  Within-compartment edges: {same_comp}/{total_edges} "
      f"({same_comp/total_edges*100:.1f}%)")
print(f"  Cross-compartment edges: {total_edges - same_comp}/{total_edges} "
      f"({(total_edges - same_comp)/total_edges*100:.1f}%)")

t5_end = time.time()
print(f"⏱ Compartment scaffold: {t5_end - t5:.1f}s")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 6: PLOTS
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("PHASE 6: Plots")
print("─" * 70)

plots_dir = OUT / "plots"

# 6a. Hex grid spatial map
fig, ax = plt.subplots(figsize=(10, 8))
cmap = plt.cm.viridis
bin_sizes = merged_full.groupby("hex_bin").size()
bin_centers_microns = {
    bin_id: (merged_full[merged_full["hex_bin"] == bin_id]["cx_microns"].mean(),
             merged_full[merged_full["hex_bin"] == bin_id]["cy_microns"].mean())
    for bin_id in bin_sizes.index
}
for bin_id in bin_sizes.index:
    cx, cy = bin_centers_microns[bin_id]
    size = bin_sizes[bin_id]
    hex_radius = HEX_DIAMETER_MICRONS / 2 * MICRONS_TO_PX * 0.7
    hex = RegularPolygon((cx, cy), numVertices=6,
                         radius=HEX_DIAMETER_MICRONS/2 * MICRONS_TO_PX * 0.7,
                         orientation=np.radians(30),
                         color=cmap(size / bin_sizes.max()),
                         ec='none', alpha=0.8)
    ax.add_patch(hex)

# Plot cell centroids on top
sc = ax.scatter(cx_all * MICRONS_TO_PX, cy_all * MICRONS_TO_PX,
                c='white', s=1, alpha=0.3)
ax.set_aspect('equal')
ax.set_xlabel("X (pixels)")
ax.set_ylabel("Y (pixels)")
ax.set_title(f"Hex Grid (Ø={HEX_DIAMETER_MICRONS}µm): "
             f"{n_bins} bins, {len(merged_full)} cells")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, bin_sizes.max()))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, label="Cells per bin")
fig.savefig(plots_dir / "hex_grid_map.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  hex_grid_map.png ✓")

# 6b. Multimodal PCA (image + expression combined)
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# PCA on image features only
pca_img = PCA(n_components=2)
img_feat_clean = np.nan_to_num(
    merged_full[[c for c in merged_full.columns
                 if c.startswith(('morph_', 'texture_', 'microenv_'))
                 and c != 'compartment_label']].values,
    nan=0.0)
pca_img_result = pca_img.fit_transform(img_feat_clean)
axes[0].scatter(pca_img_result[:, 0], pca_img_result[:, 1],
                c=merged_full["compartment_label"].map(
                    {"epithelium": "#e41a1c", "stroma": "#377eb8",
                     "necrosis": "#4daf4a", "border": "#ff7f00",
                     "unclassified": "#999999"}),
                s=10, alpha=0.6)
axes[0].set_title(f"Image Features Only\n"
                  f"({pca_img.explained_variance_ratio_[0]*100:.1f}% + "
                  f"{pca_img.explained_variance_ratio_[1]*100:.1f}%)")

# PCA on expression only
pca_expr = PCA(n_components=2)
pca_expr_result = pca_expr.fit_transform(
    np.nan_to_num(expr_smoothed, nan=0.0))
axes[1].scatter(pca_expr_result[:, 0], pca_expr_result[:, 1],
                c=merged_full["compartment_label"].map(
                    {"epithelium": "#e41a1c", "stroma": "#377eb8",
                     "necrosis": "#4daf4a", "border": "#ff7f00",
                     "unclassified": "#999999"}),
                s=10, alpha=0.6)
axes[1].set_title(f"Expression (smoothed)\n"
                  f"({pca_expr.explained_variance_ratio_[0]*100:.1f}% + "
                  f"{pca_expr.explained_variance_ratio_[1]*100:.1f}%)")

# PCA on multimodal
axes[2].scatter(pca_mm_result[:, 0], pca_mm_result[:, 1],
                c=merged_full["compartment_label"].map(
                    {"epithelium": "#e41a1c", "stroma": "#377eb8",
                     "necrosis": "#4daf4a", "border": "#ff7f00",
                     "unclassified": "#999999"}),
                s=10, alpha=0.6)
axes[2].set_title(f"Multimodal (Image + Expression)\n"
                  f"({pca_mm.explained_variance_ratio_[0]*100:.1f}% + "
                  f"{pca_mm.explained_variance_ratio_[1]*100:.1f}%)")

for ax in axes:
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

# Add compartment legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#e41a1c", label="Epithelium"),
    Patch(facecolor="#377eb8", label="Stroma"),
    Patch(facecolor="#4daf4a", label="Necrosis"),
    Patch(facecolor="#ff7f00", label="Border"),
]
fig.legend(handles=legend_elements, loc="lower center",
           ncol=5, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(plots_dir / "multimodal_pca_comparison.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  multimodal_pca_comparison.png ✓")

# 6c. Compartment spatial map
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Image with compartments overlay
comp_colors = {
    "epithelium": "#e41a1c", "stroma": "#377eb8",
    "necrosis": "#4daf4a", "border": "#ff7f00", "unclassified": "#999999"
}
comp_overlay = np.zeros((*image.shape, 3), dtype=np.uint8)
for comp, color in comp_colors.items():
    if comp == "epithelium":
        comp_overlay[epithelium] = [228, 26, 28]
    elif comp == "stroma":
        comp_overlay[stroma] = [55, 126, 184]
    elif comp == "necrosis":
        comp_overlay[necrosis] = [77, 175, 74]

axes[0].imshow(image, cmap="gray", alpha=0.6)
axes[0].imshow(comp_overlay, alpha=0.4)
axes[0].set_title("Tissue Compartments (Otsu)")
axes[0].set_xlabel("X (pixels)")
axes[0].set_ylabel("Y (pixels)")

# Cells colored by compartment
for comp, color in comp_colors.items():
    comp_cells = merged_full[merged_full["compartment_label"] == comp]
    if len(comp_cells) == 0:
        continue
    # Convert microns to pixels for plotting
    px_x = comp_cells["cx_microns"].values * MICRONS_TO_PX - x1
    px_y = comp_cells["cy_microns"].values * MICRONS_TO_PX - y1
    axes[1].scatter(px_x, px_y, c=color, s=8, alpha=0.5, label=comp, edgecolors='none')
axes[1].set_title(f"Cells by Compartment ({len(merged_full)} cells)")
axes[1].set_xlabel("X (pixels)")
axes[1].set_ylabel("Y (pixels)")
axes[1].legend()
axes[1].set_xlim(0, CROP_SIZE)
axes[1].set_ylim(CROP_SIZE, 0)

fig.tight_layout()
fig.savefig(plots_dir / "compartment_spatial_map.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  compartment_spatial_map.png ✓")

# 6d. kNN connectivity by compartment
fig, ax = plt.subplots(figsize=(8, 6))
compartment_labels_sorted = ["epithelium", "stroma", "border", "necrosis", "unclassified"]
connectivity_matrix = np.zeros((5, 5))
for i, ci in enumerate(compartment_labels_sorted):
    for j, cj in enumerate(compartment_labels_sorted):
        if ci not in comp_counts.index or cj not in comp_counts.index:
            continue
        cells_i = merged_full[merged_full["compartment_label"] == ci].index
        cells_j = merged_full[merged_full["compartment_label"] == cj].index
        if len(cells_i) == 0 or len(cells_j) == 0:
            continue
        # Count edges from i to j
        n_edges = 0
        for idx in cells_i:
            for nbr in indices[idx]:
                if nbr in cells_j:
                    n_edges += 1
        connectivity_matrix[i, j] = n_edges

# Normalize by rows
row_sums = connectivity_matrix.sum(axis=1, keepdims=True)
row_sums[row_sums == 0] = 1
connectivity_matrix_norm = connectivity_matrix / row_sums

im = ax.imshow(connectivity_matrix_norm, cmap="YlOrRd", vmin=0, vmax=1)
ax.set_xticks(range(5))
ax.set_yticks(range(5))
ax.set_xticklabels(compartment_labels_sorted, rotation=45)
ax.set_yticklabels(compartment_labels_sorted)
ax.set_xlabel("To")
ax.set_ylabel("From")
ax.set_title("kNN Connectivity Between Compartments\n(normalized by row)")
plt.colorbar(im, ax=ax, label="Fraction of edges")
fig.tight_layout()
fig.savefig(plots_dir / "compartment_knn_connectivity.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  compartment_knn_connectivity.png ✓")

# 6e. Hex grid feature summary
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Morph area per hex bin
hex_agg_sorted = hex_agg.sort_values("hex_bin")
im1 = axes[0, 0].scatter(hex_agg_sorted["cx_microns"],
                          hex_agg_sorted["cy_microns"],
                          c=hex_agg_sorted["morph_area"],
                          s=hex_agg_sorted["n_cells"] * 5,
                          cmap="viridis", alpha=0.7)
axes[0, 0].set_title("morph_area per bin")
plt.colorbar(im1, ax=axes[0, 0])

# expr_total_counts per hex bin
if "expr_total_counts" in hex_agg_sorted.columns:
    im2 = axes[0, 1].scatter(hex_agg_sorted["cx_microns"],
                              hex_agg_sorted["cy_microns"],
                              c=hex_agg_sorted["expr_total_counts"],
                              s=hex_agg_sorted["n_cells"] * 5,
                              cmap="plasma", alpha=0.7)
    axes[0, 1].set_title("expr_total_counts per bin")
    plt.colorbar(im2, ax=axes[0, 1])

# N_cells per bin
im3 = axes[1, 0].scatter(hex_agg_sorted["cx_microns"],
                          hex_agg_sorted["cy_microns"],
                          c=hex_agg_sorted["n_cells"],
                          s=hex_agg_sorted["n_cells"] * 5,
                          cmap="Blues", alpha=0.7)
axes[1, 0].set_title("Cells per bin")
plt.colorbar(im3, ax=axes[1, 0])

# Microenvironment density per bin
if "microenv_cell_density" in hex_agg_sorted.columns:
    im4 = axes[1, 1].scatter(hex_agg_sorted["cx_microns"],
                              hex_agg_sorted["cy_microns"],
                              c=hex_agg_sorted["microenv_cell_density"],
                              s=hex_agg_sorted["n_cells"] * 5,
                              cmap="RdYlBu_r", alpha=0.7)
    axes[1, 1].set_title("Cell density per bin (microenv)")
    plt.colorbar(im4, ax=axes[1, 1])

for ax in axes.flatten():
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
fig.tight_layout()
fig.savefig(plots_dir / "hex_grid_features.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  hex_grid_features.png ✓")

# ════════════════════════════════════════════════════════════════════════════
# PHASE 7: SUMMARY
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 70)
print("SUMMARY")
print("─" * 70)

t_total = time.time() - t_start
print(f"""
┌──────────────────────────────────────────────────────────────┐
│                 INTEGRATION TEST RESULTS                      │
├──────────────────────────────────────────────────────────────┤
│  Crop region:          ({CROP_SIZE}×{CROP_SIZE}) px at ({CROP_Y0},{CROP_X0})     │
│  Cells processed:      {n_cells_rasterized}                                   │
│  Features extracted:   {len(feature_cols_final)} (image) + {n_top_var} (RNA top genes)    │
│  Hex bins:             {n_bins} (Ø={HEX_DIAMETER_MICRONS} µm)                    │
│  kNN graph:            k={N_NEIGHBORS_KNN}, {
len(merged_full) * N_NEIGHBORS_KNN} edges                     │
│  Within-compartment    {same_comp/total_edges*100:.1f}%                           │
│    connectivity:                                             │
│  Multimodal PCA:       {pca_mm.explained_variance_ratio_[0]*100:.1f}% + {
pca_mm.explained_variance_ratio_[1]*100:.1f}%            │
│                                                              │
│  Output:              {OUT}                │
│  Total time:          {t_total:.1f}s                                    │
└──────────────────────────────────────────────────────────────┘
""")

# Save summary JSON
summary = {
    "crop": f"{CROP_SIZE}×{CROP_SIZE} at ({CROP_Y0},{CROP_X0})",
    "cells_processed": n_cells_rasterized,
    "n_image_features": len(feature_cols_final),
    "n_expression_genes": n_top_var,
    "hex_bins": n_bins,
    "hex_diameter_um": HEX_DIAMETER_MICRONS,
    "knn_k": N_NEIGHBORS_KNN,
    "knn_edges": len(merged_full) * N_NEIGHBORS_KNN,
    "within_compartment_connectivity_pct": round(same_comp/total_edges*100, 1),
    "multimodal_pca_pc1_pct": round(pca_mm.explained_variance_ratio_[0]*100, 1),
    "multimodal_pca_pc2_pct": round(pca_mm.explained_variance_ratio_[1]*100, 1),
    "scaler_img": "scaler_img.pkl",
    "scaler_expr": "scaler_expr.pkl",
    "total_time_s": round(t_total, 1),
}
with open(OUT / "summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print(f"✅ Test complete. Results saved to {OUT}/")
