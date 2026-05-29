---
name: image-agent
description: "Trigger: image, imagen, morphology, morfología, feature extraction, cell shape, GLCM, texture, tissue compartments. Extract morphological features from tissue images using classical computer vision (scikit-image)."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Image Agent — Image Feature Extraction

## Purpose

Extract quantitative morphological features from tissue images in spatial transcriptomics data (Xenium, Visium, MERFISH, H&E, DAPI) using **classical computer vision** (scikit-image). Three levels of extraction: cell morphology, microenvironment, and tissue compartments. CPU-only by default with optional cellpose for GPU segmentation.

## Input Format

```python
import spatialdata as sd
sdata = sd.read_zarr("path/to/data.zarr")
```

## Image Validation

Before processing, validate image and mask keys exist:

```python
valid_images = {}
for key in sdata.images:
    img = sdata[key]
    if img.ndim not in (2, 3):
        print(f"WARNING: {key} has shape {img.shape}, skipping")
        continue
    valid_images[key] = img
    print(f"Loaded: {key}, shape={img.shape}, dtype={img.dtype}")
```

Fail with descriptive error if a referenced image key is missing or empty. Infer image type from metadata — nucleic acid (DAPI) vs brightfield (H&E) vs fluorescence — for appropriate feature extraction paths.

## Processing Levels

### Level 1 — Cell Morphology

Use `skimage.measure.regionprops_table` on segmented masks. Use `skimage.feature.graycomatrix` for per-cell texture.

```python
from skimage.measure import regionprops_table
from skimage.feature import graycomatrix, graycoprops
import numpy as np

props = regionprops_table(
    mask, intensity_image=image,
    properties=('label', 'area', 'eccentricity', 'solidity',
                'perimeter', 'extent', 'equivalent_diameter',
                'centroid', 'bbox')
)

glcm_features = {f"texture_{p}": [] for p in
    ['contrast', 'dissimilarity', 'homogeneity', 'energy', 'correlation']}

for cell_id in range(1, mask.max() + 1):
    cell_mask = mask == cell_id
    if cell_mask.sum() < min_cell_area:  # default 16 px²
        continue
    rows, cols = np.where(cell_mask)
    y1, x1, y2, x2 = rows.min(), cols.min(), rows.max() + 1, cols.max() + 1
    cell_intensity = image[y1:y2, x1:x2]
    glcm = graycomatrix(
        (cell_intensity * 255).astype(np.uint8),
        distances=[1], angles=[0, np.pi/4, np.pi/2, 3*np.pi/4],
        symmetric=True, normed=True
    )
    for prop in ['contrast', 'dissimilarity', 'homogeneity', 'energy', 'correlation']:
        glcm_features[f"texture_{prop}"].append(graycoprops(glcm, prop).mean())
```

Columns: `morph_area`, `morph_eccentricity`, `morph_solidity`, `morph_perimeter`, `morph_extent`, `morph_equiv_diameter`, `texture_contrast`, `texture_dissimilarity`, `texture_homogeneity`, `texture_energy`, `texture_correlation`.

### Level 2 — Microenvironment

Use `scipy.spatial.cKDTree` for spatial statistics.

```python
from scipy.spatial import cKDTree
from scipy.ndimage import distance_transform_edt

tree = cKDTree(centroids)

# Cell density within configurable radius (default 50 µm)
density = np.array([len(tree.query_ball_point(c, r=density_radius)) - 1
                    for c in centroids])

# Distance to nearest mask boundary
boundary_dist = distance_transform_edt(mask)[centroids_y, centroids_x]

# N nearest neighbor distances (N configurable, default 5)
nndist, _ = tree.query(centroids, k=n_neighbors + 1)
nndist = nndist[:, 1:]  # exclude self
```

Columns: `microenv_cell_density`, `microenv_edge_distance`, `microenv_nndist_1` … `microenv_nndist_N`.

Supports multiple radii for density estimation via `density_radii` parameter.

### Level 3 — Tissue Compartments

Use Otsu thresholding and morphology to derive compartments without pre-existing annotations.

```python
from skimage.filters import threshold_otsu
from skimage.morphology import closing, disk, dilation

thresh = threshold_otsu(intensity_image)
epithelium = closing(intensity_image > thresh, disk(15))
stroma = closing(intensity_image <= thresh, disk(15))
border = dilation(epithelium, disk(5)) & stroma
necrosis = intensity_image < (thresh * 0.3)

for cy, cx in centroids_int:
    if epithelium[cy, cx]:
        label = "epithelium"
    elif border[cy, cx]:
        label = "border"
    elif stroma[cy, cx]:
        label = "stroma"
    else:
        label = "necrosis"
```

Fall back with warning if Otsu is unstable (e.g., uniform tissue).

Columns: `compartment_label`, `compartment_boundary_distance`, `compartment_intensity_gradient`.

## Multi-Image Merge

When multiple images exist (H&E + DAPI + FOVs), merge all features by centroid (x, y) in physical units.

```python
merged = None
for image_key, df in per_image_features.items():
    df = df.add_prefix(f"{image_key}_")
    if merged is None:
        merged = df
    else:
        merged = merged.merge(
            df, on=[f"{image_key}_centroid_x", f"{image_key}_centroid_y"],
            how="outer", suffixes=("", f"_{image_key}")
        )
```

- One-to-many: average features, log the count
- Validate coordinate systems are compatible before merging
- Prefix columns by source image key to prevent name collisions

## Output

```python
# Save features.csv (indexed by cell barcode)
features_df.to_csv("output/{session_id}/features.csv", index=True)

# Store in sdata["table"].obsm["image_features"] (dense numpy)
sdata["table"].obsm["image_features"] = features_df.values

# Key metrics in .obs
sdata["table"].obs["n_image_features"] = features_df.shape[1]
sdata["table"].obs["mean_cell_area"] = features_df.filter(like="morph_area").mean(axis=1)
sdata["table"].obs["mean_eccentricity"] = features_df.filter(like="morph_eccentricity").mean(axis=1)
if "compartment_label" in features_df.columns:
    sdata["table"].obs["tissue_compartment"] = features_df["compartment_label"]
```

```python
{
    "stage": "image-agent",
    "status": "completed",
    "n_images_processed": int,
    "n_features": int,
    "levels_completed": list[str],
    "features_csv": "output/{session_id}/features.csv",
    "summary": f"Extracted {n_features} features from {n_images} images",
    "plots": ["features_pca.png", "feature_distributions.png", "spatial_pc1.png"],
}
```

Preserve cell ordering: row i in `obsm["image_features"]` = `sdata["table"].obs_names[i]`.

## QC Plot Generation

```python
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# 1. Feature PCA scatter
pca = PCA(n_components=2)
pca_result = pca.fit_transform(features_df.select_dtypes(include=np.number))
plt.scatter(pca_result[:, 0], pca_result[:, 1], s=1, alpha=0.5)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("Image Feature PCA")
plt.savefig("output/{session_id}/plots/features_pca.png")

# 2. Feature distributions
features_df.hist(bins=50, figsize=(15, 10))
plt.savefig("output/{session_id}/plots/feature_distributions.png")

# 3. PC1 spatial overlay
plt.scatter(centroids[:, 0], centroids[:, 1], c=pca_result[:, 0],
            s=5, cmap="RdBu_r")
plt.colorbar(label="PC1")
plt.title("PC1 Spatial Overlay")
plt.savefig("output/{session_id}/plots/spatial_pc1.png")
```

## Hard Rules

- **NEVER modify** original images, masks, or `sdata["table"]` in place
- CPU-only by default: scikit-image is the only requirement
- Cellpose is **OPTIONAL** — for GPU segmentation only, never required
- All outputs go to a **NEW** output path
- For Visium data: one row per spot, coordinates from `sdata["table"].obsm["spatial"]`
- Skip cells smaller than `min_cell_area` (default 16 px²)
- Fail with descriptive error if referenced image key is missing or empty

## Validation

### Manual Validation Steps

1. **Xenium test data**: Run with known dataset with nuclei masks. Verify `morph_area` values (50–200 µm² for human cells). GLCM texture values in [0, 1). `features.csv` row count matches cell count.

2. **Visium test data**: Run on Visium (spot-level). Verify one row per spot, coordinates match `sdata["table"].obsm["spatial"]`.

3. **Reproducibility**: Run twice on same data → identical features (skimage is deterministic). Compare with `md5sum features.csv`.

4. **CPU-only test**: Run with cellpose not installed. Level 1–3 must complete without errors.

5. **Read-only verification**: Run `md5sum` on original image files before and after. Must be identical.

6. **Output shape check**: `sdata["table"].obsm["image_features"].shape[0] == sdata["table"].n_obs`. `.obs` columns populated: `n_image_features`, `mean_cell_area`, `mean_eccentricity`, `tissue_compartment` (if Level 3 ran).
