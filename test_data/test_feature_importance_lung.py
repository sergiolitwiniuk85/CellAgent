#!/usr/bin/env python3
"""
Feature Importance: Image Features + RNA Expression → Leiden Clusters
======================================================================
Train Random Forest classifiers to predict Leiden clusters from multimodal
features, extract feature importances (Gini + permutation), per-cluster one-vs-rest
signatures, spatial expression overlays, and generate a comprehensive PDF report.

Dataset: Lung Full FOV (Xenium 2 FOV)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, f1_score, classification_report
from sklearn.manifold import TSNE
from skimage import io
import json
import sys
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# ── Config ──────────────────────────────────────────────────────────────────

DATA_DIR = Path("lung_2fov/output_fullfov")
ALIGNED_CSV = DATA_DIR / "tables" / "aligned_cell_data.csv"
LEIDEN_CSV = DATA_DIR / "tables" / "leiden_clusters.csv"
OUTPUT_DIR = DATA_DIR / "feature_importance"
PLOT_DIR = OUTPUT_DIR / "plots"
TABLE_DIR = OUTPUT_DIR / "tables"
MORPHOLOGY_TIF = Path("lung_2fov/morphology.ome.tif")

LEIDEN_RES = "leiden_r0.5"
N_ESTIMATORS = 500
CV_FOLDS = 5
PERM_REPEATS = 5
RANDOM_STATE = 42
MIN_CLUSTER_SIZE = 30
TOP_N_GLOBAL = 20
TOP_N_PER_CLUSTER = 5
TOP_N_SPATIAL_GENES = 5

DATASET_NAME = "Lung Full FOV"

# ── Data Loading ────────────────────────────────────────────────────────────

def load_and_merge(aligned_path, leiden_path):
    """Load and merge aligned cell data with Leiden cluster assignments."""
    aligned = pd.read_csv(aligned_path)
    leiden = pd.read_csv(leiden_path)
    df = aligned.merge(leiden, on="cell_id")
    print(f"  Loaded {len(df)} cells, {df.shape[1]} columns")
    return df


def prepare_xy(df, leiden_col):
    """Split into feature matrix X and target y, filtering small clusters."""
    leiden_cols = [c for c in df.columns if c.startswith("leiden_")]
    exclude = {"cell_id", "cx_microns", "cy_microns"} | set(leiden_cols)

    feature_cols = [c for c in df.columns if c not in exclude]
    X = df[feature_cols].copy()
    y = df[leiden_col].copy()

    # Filter clusters below minimum size
    counts = y.value_counts()
    small = counts[counts < MIN_CLUSTER_SIZE].index.tolist()
    if small:
        print(f"  Warning: excluding clusters < {MIN_CLUSTER_SIZE} cells: {small}")
        mask = ~y.isin(small)
        X = X[mask]
        y = y[mask]

    return X.values, y.values, feature_cols, X.index.tolist()


# ── Training ────────────────────────────────────────────────────────────────

def train_rf(X, y):
    """Train Random Forest with balanced class weights."""
    clf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        class_weight="balanced",
        random_state=RANDOM_STATE,
        n_jobs=-1,
    )
    clf.fit(X, y)
    return clf


def evaluate_model(clf, X, y):
    """Cross-validation and baseline comparison."""
    cv_acc = cross_val_score(clf, X, y, cv=CV_FOLDS, scoring="accuracy")
    cv_f1 = cross_val_score(clf, X, y, cv=CV_FOLDS, scoring="f1_weighted")
    baseline = pd.Series(y).value_counts(normalize=True).max()

    return {
        "accuracy_mean": float(f"{cv_acc.mean():.4f}"),
        "accuracy_std": float(f"{cv_acc.std():.4f}"),
        "f1_mean": float(f"{cv_f1.mean():.4f}"),
        "f1_std": float(f"{cv_f1.std():.4f}"),
        "baseline_accuracy": float(f"{baseline:.4f}"),
        "beats_baseline": bool(cv_acc.mean() > baseline),
        "n_classes": int(pd.Series(y).nunique()),
        "n_samples": len(y),
        "n_features": X.shape[1],
    }


# ── Feature Importance ──────────────────────────────────────────────────────

def global_importance(clf, X_val, y_val, feature_names):
    """Gini + permutation importance."""
    gini = pd.DataFrame({
        "feature": feature_names,
        "importance_gini": clf.feature_importances_,
    }).sort_values("importance_gini", ascending=False)
    gini["rank_gini"] = range(1, len(gini) + 1)

    perm = permutation_importance(
        clf, X_val, y_val, n_repeats=PERM_REPEATS,
        random_state=RANDOM_STATE, n_jobs=1,
    )
    perm_df = pd.DataFrame({
        "feature": feature_names,
        "importance_perm_mean": perm.importances_mean,
        "importance_perm_std": perm.importances_std,
    }).sort_values("importance_perm_mean", ascending=False)
    perm_df["rank_perm"] = range(1, len(perm_df) + 1)

    merged = gini.merge(perm_df, on="feature")
    return merged


def split_modalidad(imp_df):
    """Split feature importance into image vs expression categories."""
    img = [c for c in imp_df["feature"] if not c.startswith("expr_")]
    expr = [c for c in imp_df["feature"] if c.startswith("expr_")]

    def _pct(col):
        img_s = imp_df.loc[imp_df["feature"].isin(img), col].sum()
        expr_s = imp_df.loc[imp_df["feature"].isin(expr), col].sum()
        total = img_s + expr_s
        return (float(img_s / total * 100) if total > 0 else 0,
                float(expr_s / total * 100) if total > 0 else 0)

    img_g, expr_g = _pct("importance_gini")
    img_p, expr_p = _pct("importance_perm_mean")

    return {
        "image_gini_pct": img_g, "expression_gini_pct": expr_g,
        "image_perm_pct": img_p, "expression_perm_pct": expr_p,
        "n_image_features": len(img), "n_expression_features": len(expr),
    }


def per_cluster_importance(X, y, feature_names, cluster_labels):
    """One-vs-rest RF per cluster. Returns DataFrame of top features per cluster."""
    rows = []
    for label in sorted(cluster_labels):
        y_bin = (y == label).astype(int)
        if y_bin.sum() < MIN_CLUSTER_SIZE:
            continue
        clf = RandomForestClassifier(
            n_estimators=200, class_weight="balanced",
            random_state=RANDOM_STATE, n_jobs=-1,
        )
        clf.fit(X, y_bin)
        imp = pd.DataFrame({
            "feature": feature_names,
            "importance": clf.feature_importances_,
        }).sort_values("importance", ascending=False).head(TOP_N_PER_CLUSTER)
        for _, r in imp.iterrows():
            rows.append({"cluster": int(label), "feature": r["feature"],
                         "importance": float(r["importance"])})
    return pd.DataFrame(rows)


# ── Spatial Helpers ─────────────────────────────────────────────────────────

def load_morphology_background(tif_path):
    """Load the first Z-plane of morphology.ome.tif for spatial backgrounds."""
    img = io.imread(str(tif_path))
    if img.ndim == 3:
        img = img[0]  # first Z-plane
    return img


def spatial_scatter(ax, df, coord_x, coord_y, c, cmap="tab20", s=5, alpha=0.7, **kw):
    """Scatter plot on morphology background."""
    # Background: DAPI
    ax.imshow(img, cmap="gray", alpha=0.15, extent=[x_min, x_max, y_max, y_min])
    sc = ax.scatter(df[coord_x], df[coord_y], c=c, cmap=cmap, s=s, alpha=alpha, **kw)
    ax.set_aspect("equal")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    return sc


def morphology_extent(tif_path):
    """Estimate spatial extent from morphology TIF dimensions."""
    from skimage import io
    img = io.imread(str(tif_path))
    if img.ndim == 3:
        h, w = img.shape[1:]
    else:
        h, w = img.shape
    mic_per_px = 0.2125
    return 0, w * mic_per_px, h * mic_per_px, 0  # x_min, x_max, y_max, y_min


# ── Plotting ────────────────────────────────────────────────────────────────

def plot_top20_features(imp_df, path):
    """Top-20 features by permutation importance."""
    top = imp_df.head(20)
    fig, ax = plt.subplots(figsize=(10, 8))
    y_pos = range(len(top))
    ax.barh(y_pos, top["importance_perm_mean"], xerr=top["importance_perm_std"],
            color="steelblue", alpha=0.8, label="Permutation importance")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top["feature"])
    ax.invert_yaxis()
    ax.set_xlabel("Importance (mean ± std, 5 reps)")
    ax.set_title("Top-20 Features — Permutation Importance")
    ax.legend(loc="lower right")
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ top20_features.png")


def plot_modalidad_split(split, path):
    """Pie chart split: image vs expression importance."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, key, title in [
        (axes[0], "gini", "Gini Importance"),
        (axes[1], "perm", "Permutation Importance"),
    ]:
        img_p = split[f"image_{key}_pct"]
        expr_p = split[f"expression_{key}_pct"]
        ax.pie([img_p, expr_p],
               labels=[f"Image ({img_p:.0f}%)", f"Expression ({expr_p:.0f}%)"],
               autopct="%1.1f%%", colors=["#66c2a5", "#fc8d62"],
               startangle=90)
        ax.set_title(title)
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ importance_by_modalidad.png")


def plot_per_cluster_heatmap(pc_imp, path):
    """Heatmap: clusters × top features."""
    pivot = pc_imp.pivot_table(index="cluster", columns="feature",
                               values="importance", fill_value=0)
    fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 0.8),
                                    max(6, pivot.shape[0] * 0.6)))
    sns.heatmap(pivot, ax=ax, cmap="YlOrRd", cbar_kws={"label": "Importance"})
    ax.set_title("Per-Cluster Feature Importance (one-vs-rest RF)")
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ per_cluster_importance.png")


def plot_spatial_clusters(df, coord_x, coord_y, cluster_col, path, title="Leiden Clusters — Spatial Map"):
    """Spatial scatter colored by cluster."""
    fig, ax = plt.subplots(figsize=(12, 10))
    n_clusters = df[cluster_col].nunique()
    sc = ax.scatter(df[coord_x], df[coord_y],
                    c=df[cluster_col], cmap="tab20", s=5, alpha=0.7)
    ax.set_aspect("equal")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title(title)
    cbar = plt.colorbar(sc, ax=ax, ticks=range(n_clusters))
    cbar.set_label("Leiden Cluster")
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ spatial_clusters.png")


def plot_spatial_top_feature(df, coord_x, coord_y, feature_cols, imp_df, path):
    """Spatial: each cell colored by its most important expressed feature."""
    top_n = min(10, len(imp_df))
    top_feats = imp_df.head(top_n)["feature"].tolist()
    # Only keep features present in the dataframe
    top_feats = [f for f in top_feats if f in df.columns]

    df["_top_feat"] = df[top_feats].idxmax(axis=1)

    fig, ax = plt.subplots(figsize=(12, 10))
    cmap = plt.colormaps["tab10"]
    for i, feat in enumerate(top_feats):
        mask = df["_top_feat"] == feat
        if mask.sum() == 0:
            continue
        ax.scatter(df.loc[mask, coord_x], df.loc[mask, coord_y],
                   c=[cmap(i % 10)], s=5, alpha=0.7, label=feat)
    ax.set_aspect("equal")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title("Dominant Feature per Cell (Top-10)")
    ax.legend(markerscale=5, fontsize=6, loc="upper right")
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    df.drop(columns=["_top_feat"], inplace=True)
    print(f"    ✓ spatial_top_feature.png")


def plot_spatial_expression_top5(df, coord_x, coord_y, imp_df, path):
    """Spatial expression of the top-5 most important genes."""
    # Filter to expression features
    expr_feats = [c for c in imp_df["feature"] if c.startswith("expr_")]
    if len(expr_feats) == 0:
        print("    ⚠ No expression features in top list, skipping top5 spatial")
        # Create empty placeholder
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No expression features in top list",
                ha="center", va="center", transform=ax.transAxes)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return

    top_expr = imp_df[imp_df["feature"].isin(expr_feats)].head(TOP_N_SPATIAL_GENES)
    n_genes = len(top_expr)

    fig, axes = plt.subplots(1, n_genes, figsize=(5 * n_genes, 5))
    if n_genes == 1:
        axes = [axes]
    for ax, (_, row) in zip(axes, top_expr.iterrows()):
        gene = row["feature"]
        sc = ax.scatter(df[coord_x], df[coord_y],
                        c=df[gene] if gene in df.columns else 0,
                        cmap="viridis", s=5, alpha=0.7)
        ax.set_aspect("equal")
        ax.set_xlabel("X (µm)")
        ax.set_ylabel("Y (µm)")
        ax.set_title(f"{gene} (rank {int(row['rank_perm'])})")
        plt.colorbar(sc, ax=ax, shrink=0.7)
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ spatial_expression_top5.png")


def plot_umap_features(X, y, feature_names, path, title="UMAP of Multimodal Features"):
    """UMAP (via TSNE as fallback) of feature space, colored by cluster."""
    print("    Computing UMAP (t-SNE fallback)...")
    # Use TSNE as we don't have UMAP installed; n_components=2
    from sklearn.manifold import TSNE
    # Subsample if too large for TSNE
    n_samples = min(5000, X.shape[0])
    if X.shape[0] > n_samples:
        rng = np.random.RandomState(RANDOM_STATE)
        idx = rng.choice(X.shape[0], n_samples, replace=False)
        X_sub = X[idx]
        y_sub = y[idx]
    else:
        X_sub = X
        y_sub = y

    tsne = TSNE(n_components=2, random_state=RANDOM_STATE, perplexity=30)
    emb = tsne.fit_transform(X_sub)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=y_sub, cmap="tab20", s=10, alpha=0.7)
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_title(title)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Leiden Cluster")
    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ feature_umap.png")


def plot_3panel_comparison(df, coord_x, coord_y, cluster_col, imp_df, path):
    """3-panel: top image feature vs top expression vs cluster."""
    # Pick top image feature and top expression feature
    img_feats = imp_df[~imp_df["feature"].str.startswith("expr_")]
    expr_feats = imp_df[imp_df["feature"].str.startswith("expr_")]

    top_img = img_feats.iloc[0]["feature"] if len(img_feats) > 0 else None
    top_expr = expr_feats.iloc[0]["feature"] if len(expr_feats) > 0 else None

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel 1: Top image feature
    ax = axes[0]
    if top_img and top_img in df.columns:
        sc = ax.scatter(df[coord_x], df[coord_y],
                        c=df[top_img], cmap="viridis", s=5, alpha=0.7)
        ax.set_title(f"Image: {top_img}")
        plt.colorbar(sc, ax=ax, shrink=0.7)
    else:
        ax.text(0.5, 0.5, "No image feature", ha="center", va="center", transform=ax.transAxes)
    ax.set_aspect("equal")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")

    # Panel 2: Top expression feature
    ax = axes[1]
    if top_expr and top_expr in df.columns:
        sc = ax.scatter(df[coord_x], df[coord_y],
                        c=df[top_expr], cmap="plasma", s=5, alpha=0.7)
        ax.set_title(f"Expression: {top_expr}")
        plt.colorbar(sc, ax=ax, shrink=0.7)
    else:
        ax.text(0.5, 0.5, "No expression feature", ha="center", va="center", transform=ax.transAxes)
    ax.set_aspect("equal")
    ax.set_xlabel("X (µm)")

    # Panel 3: Leiden clusters
    ax = axes[2]
    n_clusters = df[cluster_col].nunique()
    sc = ax.scatter(df[coord_x], df[coord_y],
                    c=df[cluster_col], cmap="tab20", s=5, alpha=0.7)
    ax.set_title("Leiden Clusters")
    cbar = plt.colorbar(sc, ax=ax, shrink=0.7, ticks=range(n_clusters))
    ax.set_aspect("equal")
    ax.set_xlabel("X (µm)")

    plt.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    ✓ 3panel_comparison.png")


# ── Report ──────────────────────────────────────────────────────────────────

def generate_report(metrics, imp_df, split, pc_imp, df, cluster_col, output_dir, dataset_name):
    """Generate Markdown report → convert to PDF."""
    report_md = f"""# Feature Importance Report — {dataset_name}

## Classification Performance

| Metric | Value |
|--------|-------|
| CV Accuracy (5-fold) | {metrics['accuracy_mean']:.3f} ± {metrics['accuracy_std']:.3f} |
| CV Weighted F1 | {metrics['f1_mean']:.3f} ± {metrics['f1_std']:.3f} |
| Baseline (majority class) | {metrics['baseline_accuracy']:.3f} |
| Beats baseline | {metrics['beats_baseline']} |
| Number of classes | {metrics['n_classes']} |
| Samples | {metrics['n_samples']} |
| Features | {metrics['n_features']} |

## Modality Split

| Modality | Gini Importance | Permutation Importance |
|----------|:--------------:|:---------------------:|
| Image features ({split['n_image_features']}) | {split['image_gini_pct']:.1f}% | {split['image_perm_pct']:.1f}% |
| Expression features ({split['n_expression_features']}) | {split['expression_gini_pct']:.1f}% | {split['expression_perm_pct']:.1f}% |

## Top-20 Features (by Permutation Importance)

| Rank | Feature | Gini Importance | Permutation (mean ± std) |
|:----:|---------|:--------------:|:------------------------:|
"""
    top20 = imp_df.head(20)
    for i, (_, row) in enumerate(top20.iterrows()):
        report_md += f"| {i+1} | {row['feature']} | {row['importance_gini']:.4f} | {row['importance_perm_mean']:.4f} ± {row['importance_perm_std']:.4f} |\n"

    report_md += f"""
## Cluster Distribution

| Cluster | Count | Percentage |
|:-------:|:-----:|:----------:|
"""
    cluster_counts = df[cluster_col].value_counts().sort_index()
    total = len(df)
    for label, count in cluster_counts.items():
        pct = count / total * 100
        report_md += f"| {label} | {count} | {pct:.1f}% |\n"

    report_md += """
## Generated Plots

![Top-20 Features](plots/top20_features.png)

![Importance by Modality](plots/importance_by_modalidad.png)

![Per-Cluster Importance](plots/per_cluster_importance.png)

![Spatial Clusters](plots/spatial_clusters.png)

![Dominant Feature per Cell](plots/spatial_top_feature.png)

![Top-5 Spatial Gene Expression](plots/spatial_expression_top5.png)

![3-Panel Comparison](plots/3panel_comparison.png)

![Feature Space UMAP](plots/feature_umap.png)
"""
    md_path = output_dir / "feature_importance_report.md"
    md_path.write_text(report_md)
    print(f"    ✓ feature_importance_report.md")

    # Convert to PDF via report-agent
    pdf_path = output_dir / "feature_importance_report.pdf"
    import subprocess
    # Compute absolute paths for PDF conversion
    # When running from test_data/, parent is the repo root
    repo_root = Path.cwd().parent
    # When running from repo root, cwd is already repo root
    if Path.cwd().name == "test_data":
        repo_root = Path.cwd().parent
    script_path = repo_root / "skills" / "report-agent" / "scripts" / "md_to_pdf.py"
    result = subprocess.run(
        ["python3", str(script_path),
         str(md_path), "-o", str(pdf_path)],
        capture_output=True, text=True, cwd=Path.cwd()
    )
    if result.returncode == 0:
        print(f"    ✓ feature_importance_report.pdf")
    else:
        print(f"    ⚠ PDF generation failed: {result.stderr[:200]}")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print(f"Feature Importance Analysis — {DATASET_NAME}")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load
    print("\n[1/6] Loading data...")
    df = load_and_merge(ALIGNED_CSV, LEIDEN_CSV)

    # 2. Prepare
    print(f"\n[2/6] Preparing X/y for {LEIDEN_RES}...")
    X, y, feature_cols, clean_idx = prepare_xy(df, LEIDEN_RES)
    # Subset df to match filtered indices
    df_clean = df.iloc[clean_idx].reset_index(drop=True)
    print(f"  X: {X.shape}, y: {len(np.unique(y))} clusters")
    print(f"  Cluster distribution:")
    for label, count in zip(*np.unique(y, return_counts=True)):
        print(f"    {label}: {count}")

    # 3. Train + Evaluate
    print(f"\n[3/6] Training RF (n_estimators={N_ESTIMATORS})...")
    clf = train_rf(X, y)
    metrics = evaluate_model(clf, X, y)
    print(f"  CV Accuracy: {metrics['accuracy_mean']:.3f} ± {metrics['accuracy_std']:.3f}")
    print(f"  CV F1: {metrics['f1_mean']:.3f} ± {metrics['f1_std']:.3f}")
    print(f"  Baseline: {metrics['baseline_accuracy']:.3f}")
    print(f"  Beats baseline: {metrics['beats_baseline']}")
    print(f"  Classes: {metrics['n_classes']}, Samples: {metrics['n_samples']}, Features: {metrics['n_features']}")

    # 4. Feature importance
    print(f"\n[4/6] Computing feature importance...")
    imp_df = global_importance(clf, X, y, feature_cols)
    split = split_modalidad(imp_df)
    print(f"  Image: {split['image_perm_pct']:.1f}%  |  Expression: {split['expression_perm_pct']:.1f}%")
    print(f"  Top-5 features:")
    for i, (_, row) in enumerate(imp_df.head(5).iterrows()):
        print(f"    {i+1}. {row['feature']} (Gini={row['importance_gini']:.4f}, Perm={row['importance_perm_mean']:.4f})")

    # 5. Per-cluster importance
    print(f"\n[5/6] Computing per-cluster importance...")
    pc_imp = per_cluster_importance(X, y, feature_cols, np.unique(y))
    print(f"  {len(pc_imp)} cluster-feature pairs across {pc_imp['cluster'].nunique()} clusters")

    # 6. Visualize
    print(f"\n[6/6] Generating plots...")

    print("  Feature importance plots...")
    plot_top20_features(imp_df, PLOT_DIR / "top20_features.png")
    plot_modalidad_split(split, PLOT_DIR / "importance_by_modalidad.png")
    plot_per_cluster_heatmap(pc_imp, PLOT_DIR / "per_cluster_importance.png")

    print("  Spatial plots...")
    plot_spatial_clusters(df_clean, "cx_microns", "cy_microns", LEIDEN_RES,
                          PLOT_DIR / "spatial_clusters.png")
    plot_spatial_top_feature(df_clean, "cx_microns", "cy_microns", feature_cols, imp_df,
                             PLOT_DIR / "spatial_top_feature.png")
    plot_spatial_expression_top5(df_clean, "cx_microns", "cy_microns", imp_df,
                                  PLOT_DIR / "spatial_expression_top5.png")
    plot_3panel_comparison(df_clean, "cx_microns", "cy_microns", LEIDEN_RES, imp_df,
                           PLOT_DIR / "3panel_comparison.png")

    print("  Embedding plot...")
    plot_umap_features(X, y, feature_cols, PLOT_DIR / "feature_umap.png")

    # Save tables
    print("  Saving tables...")
    imp_df.to_csv(TABLE_DIR / "feature_importance_global.csv", index=False)
    pc_imp.to_csv(TABLE_DIR / "feature_importance_per_cluster.csv", index=False)
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(TABLE_DIR / "classification_report.csv", index=False)

    # Summary JSON
    summary = {
        "dataset": DATASET_NAME,
        **metrics,
        "top5_features": imp_df.head(5)["feature"].tolist(),
        "modalidad_split": {k: v for k, v in split.items()},
    }
    with open(OUTPUT_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Outputs:")
    print(f"    Plots: {PLOT_DIR}")
    print(f"    Tables: {TABLE_DIR}")
    print(f"    Summary: {OUTPUT_DIR / 'summary.json'}")

    # Report
    print("\n  Generating report...")
    generate_report(metrics, imp_df, split, pc_imp, df_clean, LEIDEN_RES, OUTPUT_DIR, DATASET_NAME)

    print(f"\n{'=' * 60}")
    print(f"✓ Feature Importance Complete — {DATASET_NAME}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
