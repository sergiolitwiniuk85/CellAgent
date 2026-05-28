"""
Test script: image-agent feature extraction on Xenium Lung 2 FOV.
Uses a cropped region for fast testing of the full pipeline.
"""

import numpy as np
import pandas as pd
import tifffile
import skimage.measure
import skimage.feature
import skimage.filters
import skimage.morphology
import scipy.spatial
import scipy.ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from skimage.draw import polygon as sk_polygon
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

BASE = Path("/home/sergiolitwiniuk/gentle-projects/scai/test_data/lung_2fov")
OUT = BASE / "output"
OUT.mkdir(exist_ok=True)
(OUT / "plots").mkdir(exist_ok=True)

# Xenium scale: 0.2125 µm/pixel
MICRONS_TO_PX = 1.0 / 0.2125

# ── 1. Load image ──────────────────────────────────────────────────────
print("=" * 60)
print("Loading morphology image...")
img = tifffile.imread(str(BASE / "morphology.ome.tif"))
image = img[0] if img.ndim == 3 else img  # First Z-slice (DAPI)
print(f"  Full image: {image.shape}, dtype={image.dtype}")

# Crop to a small region with tissue (rough center of the FOV)
# Physical coords ~600x400 µm → ~2823x1882 px
# But a 500x500 px crop is enough for testing
y0, x0 = 1500, 2500
crop_size = 600
image_crop = image[y0:y0+crop_size, x0:x0+crop_size]
print(f"  Crop: ({crop_size}×{crop_size}) at ({y0},{x0})")
print(f"  Crop range: [{image_crop.min()}, {image_crop.max()}]")

# ── 2. Load cell boundaries ────────────────────────────────────────────
print("\nLoading nucleus boundaries...")
nuclei_df = pd.read_parquet(BASE / "nucleus_boundaries.parquet")
print(f"  Total cells in full data: {nuclei_df['cell_id'].nunique()}")

# ── 3. Filter cells in crop region ─────────────────────────────────────
# Convert all vertices to pixels and filter to crop
xs_all = (nuclei_df["vertex_x"].values * MICRONS_TO_PX).astype(int)
ys_all = (nuclei_df["vertex_y"].values * MICRONS_TO_PX).astype(int)

# A cell is "in crop" if its centroid would be in the crop region
# Approximate by checking if any vertex falls in crop
in_crop = ((xs_all >= x0) & (xs_all < x0 + crop_size) &
           (ys_all >= y0) & (ys_all < y0 + crop_size))
crop_cell_ids = nuclei_df["cell_id"].iloc[np.where(in_crop)[0]].unique()
print(f"  Cells in crop region: {len(crop_cell_ids)}")

# Filter to these cells
nuclei_crop = nuclei_df[nuclei_df["cell_id"].isin(crop_cell_ids)].copy()

# ── 4. Build crop mask ─────────────────────────────────────────────────
print("\nBuilding segmentation mask for crop...")
mask = np.zeros(image_crop.shape, dtype=np.int32)

# Pre-compute pixel coords for crop cells
xs_crop = (nuclei_crop["vertex_x"].values * MICRONS_TO_PX).astype(int) - x0
ys_crop = (nuclei_crop["vertex_y"].values * MICRONS_TO_PX).astype(int) - y0
nuclei_crop = nuclei_crop.assign(px_x=xs_crop, px_y=ys_crop)

# Groupby with progress
grouped = nuclei_crop.groupby("cell_id")
n_cells = len(grouped)
cell_centroids = {}
skipped = 0

for i, (cell_id, group) in enumerate(grouped):
    if len(group) < 3:
        skipped += 1
        continue
    xs = np.clip(group["px_x"].values, 0, crop_size - 1)
    ys = np.clip(group["px_y"].values, 0, crop_size - 1)
    rr, cc = sk_polygon(ys, xs)
    if len(rr) == 0:
        skipped += 1
        continue
    label = i + 1 - skipped
    mask[rr, cc] = label
    cell_centroids[cell_id] = (float(ys.mean()), float(xs.mean()))

print(f"  Rasterized {len(cell_centroids)} cells ({skipped} skipped)")
print(f"  Unique labels in mask: {np.unique(mask).size - 1}")

# Map for reverse lookup
cell_ids_list = list(cell_centroids.keys())
label_to_cell_id = {i+1: cid for i, cid in enumerate(cell_ids_list)}

if len(cell_centroids) == 0:
    print("❌ No cells in crop region. Try different coordinates.")
    exit(1)

# ── 5. LEVEL 1: Cell Morphology ────────────────────────────────────────
print("\n" + "=" * 60)
print("LEVEL 1: Cell Morphology — regionprops + GLCM")
print("-" * 60)

image_8bit = ((image_crop - image_crop.min()) /
              (image_crop.max() - image_crop.min()) * 255).astype(np.uint8)

props = skimage.measure.regionprops_table(
    mask, intensity_image=image_crop,
    properties=('label', 'area', 'eccentricity', 'solidity',
                'perimeter', 'extent', 'equivalent_diameter', 'centroid')
)

morph_df = pd.DataFrame(props)
morph_df.rename(columns={
    'label': 'cell_id',
    'area': 'morph_area',
    'eccentricity': 'morph_eccentricity',
    'solidity': 'morph_solidity',
    'perimeter': 'morph_perimeter',
    'extent': 'morph_extent',
    'equivalent_diameter': 'morph_equiv_diameter',
}, inplace=True)

print(f"  regionprops done: {len(morph_df)} cells")
print(f"  morph_area: [{morph_df['morph_area'].min():.1f}, {morph_df['morph_area'].max():.1f}] px²")
print(f"  morph_eccentricity: [{morph_df['morph_eccentricity'].min():.3f}, {morph_df['morph_eccentricity'].max():.3f}]")
print(f"  morph_solidity: [{morph_df['morph_solidity'].min():.3f}, {morph_df['morph_solidity'].max():.3f}]")

# Per-cell GLCM
print("  Computing per-cell GLCM textures...")
min_cell_area = 16
glcm_data = {p: [] for p in ['contrast', 'dissimilarity', 'homogeneity', 'energy', 'correlation']}

for _, row in morph_df.iterrows():
    label = int(row['cell_id'])
    cell_mask = mask == label
    if cell_mask.sum() < min_cell_area:
        for p in glcm_data:
            glcm_data[p].append(np.nan)
        continue
    ys, xs = np.where(cell_mask)
    y1, x1, y2, x2 = ys.min(), xs.min(), ys.max() + 1, xs.max() + 1
    cell_intensity = image_8bit[y1:y2, x1:x2]
    glcm = skimage.feature.graycomatrix(
        cell_intensity,
        distances=[1], angles=[0, np.pi/4, np.pi/2, 3*np.pi/4],
        symmetric=True, normed=True
    )
    for p in glcm_data:
        glcm_data[p].append(skimage.feature.graycoprops(glcm, p).mean())

for p, vals in glcm_data.items():
    morph_df[f"texture_{p}"] = vals
    valid = [v for v in vals if not np.isnan(v)]
    if valid:
        print(f"  texture_{p}: [{min(valid):.4f}, {max(valid):.4f}]")

# Map back to string cell_id
morph_df["cell_id_str"] = morph_df["cell_id"].map(label_to_cell_id)

# ── 6. LEVEL 2: Microenvironment ───────────────────────────────────────
print("\n" + "=" * 60)
print("LEVEL 2: Microenvironment — cKDTree")
print("-" * 60)

centroids = morph_df[['centroid-0', 'centroid-1']].values
centroids_xy = centroids[:, [1, 0]]
tree = scipy.spatial.cKDTree(centroids_xy)

density_radius = 50  # pixels
density = np.array([len(tree.query_ball_point(c, r=density_radius)) - 1
                    for c in centroids_xy])
morph_df["microenv_cell_density"] = density
print(f"  Cell density (r={density_radius}px): "
      f"[{density.min()}, {density.max()}] mean={density.mean():.1f}")

dist_from_boundary = scipy.ndimage.distance_transform_edt(mask > 0)
edge_dist = np.array([
    dist_from_boundary[int(row['centroid-0']), int(row['centroid-1'])]
    for _, row in morph_df.iterrows()
])
morph_df["microenv_edge_distance"] = edge_dist
print(f"  Edge distance: [{edge_dist.min():.1f}, {edge_dist.max():.1f}] "
      f"mean={edge_dist.mean():.1f} px")

n_neighbors = 5
nndist, _ = tree.query(centroids_xy, k=n_neighbors + 1)
for i in range(1, n_neighbors + 1):
    morph_df[f"microenv_nndist_{i}"] = nndist[:, i]
print(f"  NND 1: mean={nndist[:, 1].mean():.1f}")

# ── 7. LEVEL 3: Tissue Compartments ────────────────────────────────────
print("\n" + "=" * 60)
print("LEVEL 3: Tissue Compartments — Otsu")
print("-" * 60)

thresh = skimage.filters.threshold_otsu(image_crop)
print(f"  Otsu threshold: {thresh:.1f}")

epithelium = skimage.morphology.closing(image_crop > thresh,
                                        skimage.morphology.disk(15))
stroma = skimage.morphology.closing(image_crop <= thresh,
                                     skimage.morphology.disk(15))
border = skimage.morphology.dilation(epithelium,
                                     skimage.morphology.disk(5)) & stroma
necrosis = image_crop < (thresh * 0.3)

labels = []
comp_boundary_dist = []
comp_intensity_grad = []
for _, row in morph_df.iterrows():
    cy = int(min(row['centroid-0'], crop_size - 1))
    cx = int(min(row['centroid-1'], crop_size - 1))
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
    gy, gx = np.gradient(image_crop.astype(float))
    comp_intensity_grad.append(np.sqrt(gy[cy, cx]**2 + gx[cy, cx]**2))

morph_df["compartment_label"] = labels
morph_df["compartment_boundary_distance"] = comp_boundary_dist
morph_df["compartment_intensity_gradient"] = comp_intensity_grad

comp_counts = morph_df["compartment_label"].value_counts()
print("  Compartments:")
for comp, count in comp_counts.items():
    print(f"    {comp}: {count} ({count/len(morph_df)*100:.1f}%)")

# ── 8. Save ────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("SAVING OUTPUT")
print("-" * 60)

feature_cols = [c for c in morph_df.columns
                if c.startswith(('morph_', 'texture_', 'microenv_', 'compartment_'))]
feature_df = morph_df[['cell_id_str'] + feature_cols].copy()
feature_df.rename(columns={'cell_id_str': 'cell_id'}, inplace=True)
feature_df.to_csv(OUT / "features.csv", index=False)
print(f"  features.csv: {len(feature_df)} rows × {len(feature_cols)} features")
print(f"  Feature columns: {feature_cols}")

obs_metrics = pd.DataFrame({
    "n_image_features": len(feature_cols),
    "mean_cell_area": feature_df.filter(like="morph_area").mean(axis=1),
    "mean_eccentricity": feature_df.filter(like="morph_eccentricity").mean(axis=1),
    "tissue_compartment": feature_df["compartment_label"],
})
obs_metrics.to_csv(OUT / "obs_metrics.csv", index=False)

# ── 9. QC Plots ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("QC PLOTS")
print("-" * 60)

plots_dir = OUT / "plots"
numeric_cols = feature_df.select_dtypes(include=np.number)
numeric_cols = numeric_cols.dropna(axis=1, how='any')

pca = PCA(n_components=2)
pca_result = pca.fit_transform(numeric_cols)

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(pca_result[:, 0], pca_result[:, 1], s=10, alpha=0.6, c='steelblue')
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title(f"Image Feature PCA ({len(feature_df)} cells)")
fig.savefig(plots_dir / "features_pca.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  features_pca.png ✓")

feature_subset = feature_cols[:8]
ncols, nrows = 4, (len(feature_subset) + 3) // 4
fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3 * nrows))
axes = axes.flatten()
for i, col in enumerate(feature_subset):
    axes[i].hist(feature_df[col].dropna(), bins=50, color='steelblue', alpha=0.7)
    axes[i].set_title(col)
for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)
fig.tight_layout()
fig.savefig(plots_dir / "feature_distributions.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  feature_distributions.png ✓")

fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(morph_df['centroid-1'], morph_df['centroid-0'],
                c=pca_result[:, 0], s=15, cmap="RdBu_r", alpha=0.7)
ax.set_xlabel("X (pixels)")
ax.set_ylabel("Y (pixels)")
ax.set_title("PC1 Spatial Overlay")
plt.colorbar(sc, ax=ax, label="PC1")
fig.savefig(plots_dir / "spatial_pc1.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("  spatial_pc1.png ✓")

# ── 10. Summary ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("SUMMARY ✅")
print("-" * 60)
print(f"  Cells processed: {len(feature_df)}")
print(f"  Features extracted: {len(feature_cols)}")
print(f"  Levels completed: morphology, microenvironment, compartments")
print(f"  Morph area: {feature_df['morph_area'].min():.1f} – {feature_df['morph_area'].max():.1f} px²")
print(f"  Compartments: {len(comp_counts)} tissue regions found")
print(f"  Output: {OUT}/")
