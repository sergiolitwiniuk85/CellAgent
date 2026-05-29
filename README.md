# CellAgent — AI-Powered Single-Cell Analysis Framework

CellAgent lets you run complete single-cell analyses (RNA-seq, multimodal, spatial) by describing what you want in plain English. No need to memorize scanpy, muon, or squidpy APIs — just tell CellAgent what you need, and it guides you step by step.

Behind the scenes, an orchestrator delegates each stage to a specialist agent: data loading, QC, normalization, clustering, multimodal integration, spatial analysis, image processing, and reporting. **Every analysis ALWAYS ends with two mandatory documents: an executive report (Markdown + PDF) and a full traceability document** with all commands, parameters, and results — guaranteed reproducibility.

## Features

- **Unimodal RNA-seq**: data → QC → normalize → cluster → markers → **report** → **trace**
- **Multimodal integration**: RNA + ATAC + protein via muon (WNN, MOFA+)
- **Spatial analysis**: Xenium, Visium, MERFISH via SpatialData + Squidpy
- **Image morphology pipeline**: extract cell morphology, GLCM texture, and tissue compartments from tissue images (Xenium morphology.ome.tif) — classical CV, no GPU needed
- **Multimodal image + expression integration**: hexagonal binning (50µm), shared kNN propagation, tissue compartment scaffold, Leiden clustering on the joint multimodal space
- **Feature importance**: interpret what drives multimodal clusters with Random Forest classifiers, permutation importance, per-cluster signatures, and spatial expression overlays
- **Cross-session memory**: remembers parameters and decisions across analyses via Engram
- **Report + Traceability**: every analysis produces an executive report (PDF) AND a forensic traceability document with all commands, parameters, library versions, and decisions

## Installation

### Prerequisites

You need **two things**: the scverse environment (scanpy, muon, etc.) and OpenCode (the AI agent runtime).

### Step 1: Install OpenCode (Linux & macOS)

```bash
# Option A — Install script (Linux & macOS)
curl -fsSL https://raw.githubusercontent.com/opencode-ai/opencode/refs/heads/main/install | bash

# Install a specific version
curl -fsSL https://raw.githubusercontent.com/opencode-ai/opencode/refs/heads/main/install | VERSION=0.1.0 bash

# Option B — Homebrew (macOS & Linux)
brew install opencode-ai/tap/opencode

# Verify
opencode --version
```

> **Note for macOS users**: OpenCode works natively on macOS. No Docker or extra setup needed.

If you prefer a specific version or manual install, see [opencode.ai/docs](https://opencode.ai/docs).

### Step 2: Install Engram (persistent memory)

Engram saves analysis context across sessions so CellAgent remembers past decisions.

```bash
# Linux
curl -fsSL https://engram.sh/install.sh | sh

# macOS
brew install engramhq/tap/engram
```

Start the Engram service:

```bash
engram service start
```

Verify it's running:

```bash
engram doctor
```

### Step 3: Clone CellAgent

```bash
git clone https://github.com/sergiolitwiniuk85/CellAgent.git
cd CellAgent
```

### Step 4: Install the scverse environment

**Option A — Conda (recommended for daily use)**

```bash
conda env create -f environment.yml
conda activate cellagent
```

This pins tested versions: scanpy 1.12.1, muon 0.1.7, squidpy 1.8.1, spatialdata 0.7.3.

**Option B — Docker (fully isolated)**

```bash
docker pull ghcr.io/sergiolitwiniuk85/cellagent:latest
docker run -it --rm \
  -v $(pwd):/data \
  -v /path/to/CellAgent:/cellagent \
  ghcr.io/sergiolitwiniuk85/cellagent:latest
# Inside: cd /cellagent && opencode --agent scai-orchestrator --project /data
```

**Option C — Singularity (HPC clusters)**

```bash
singularity build cellagent.sif Singularity.def
singularity run cellagent.sif
```

## Quick Start

```bash
conda activate cellagent
opencode --agent scai-orchestrator --project /path/to/your/data
```

Then tell CellAgent what you want. Examples:

> *"I have a PBMC h5ad, I want QC and clustering"*
>
> *"Load this CITE-seq data and integrate RNA with protein"*
>
> *"I have a Xenium dataset, run spatial analysis"*
> *
> *"I have a Xenium dataset, run the full multimodal image + expression pipeline"*
> *
> *"Which features drive my spatial clusters? Run feature importance"*

The orchestrator will guide you through every step, showing plots and results, and asking for confirmation before proceeding.

## What You Get

After an analysis completes, you get:

### 1. Executive Report (`output/{session}/report.md`)

A publication-ready markdown report with:
- Executive summary (cells, genes, clusters)
- QC metrics and filtering decisions
- Normalization parameters and HVG selection
- Dimensionality reduction (PCA, UMAP)
- Clustering results with resolution
- Marker genes per cluster with cell type annotation
- Complete parameter table for reproducibility

### 2. Traceability Document (`output/{session}/traceability_{project}_{date}.md`)

A forensic record of the entire analysis:
- Every stage executed in chronological order
- Exact Python commands and scanpy/muon APIs used
- All parameters with their values and justifications
- Numerical results (cells before/after, clusters, etc.)
- Paths to all generated plots
- Library versions (scanpy, muon, squidpy, etc.)
- Decisions made during the analysis

This document is **self-contained** — anyone reading it can reproduce the analysis.

### 3. Generated Files

```
output/{session_id}/
├── data_processed.h5ad       # Filtered, normalized, clustered data
├── report.md                 # Executive report
├── traceability_*.md         # Full traceability document
└── plots/
    ├── qc_violin_pre.png
    ├── qc_violin_post.png
    ├── pca_variance.png
    ├── hvg_plot.png
    ├── umap_clusters.png
    ├── marker_heatmap.png
    └── ...
```

## Example Analysis Sessions

### Example 1: Unimodal RNA-seq (PBMC)

```
User:  "I have a PBMC h5ad, I want QC and clustering"
Agent: "Got it. Let me load your data first."

  → data-agent loads pbmc.h5ad (2638 cells × 1838 genes)
  → qc-agent computes metrics, generates violin plots
      Filters: min_genes=200, max_genes=6000, max_pct_mito=20
      Retained: 2500/2638 cells (94.8%)
  → normalize-agent runs normalize_total + log1p + HVG
      2000 HVGs selected
  → cluster-agent runs PCA(50), neighbors(30), UMAP, Leiden(res=0.8)
      12 clusters identified
      Top markers: CD3D (cluster 0), CD14 (cluster 1), NKG7 (cluster 2)
  → report-agent generates executive report
  → trace-agent generates traceability document
  
Agent: "Analysis complete. Report and trace document are ready."
```

### Example 2: Multimodal Integration (CITE-seq)

```
User:  "I have CITE-seq data, I want to integrate RNA and protein"
Agent: "Loading your MuData file."

  → data-agent loads data.h5mu (RNA + protein modalities)
  → qc-agent processes both modalities
  → normalize-agent normalizes RNA; muon CLR for protein
  → integration-agent runs WNN integration
      Joint UMAP, modality weights per cell
      Leiden on WNN graph: 8 multimodal clusters
  → report-agent generates integrated report
  → trace-agent generates full traceability
  
Agent: "WNN shows good separation. 8 joint clusters found."
```

### Example 3: Spatial Analysis (Xenium)

```
User:  "I have a Xenium dataset, run spatial analysis"
Agent: "Let me load the SpatialData zarr."

  → data-agent loads sdata.zarr (5000 cells, H&E image, transcripts)
  → qc-agent filters low-quality cells
  → normalize-agent + cluster-agent runs standard pipeline
  → spatial-agent computes spatial neighbors, nhood enrichment
      Spatial scatter plots, co-occurrence analysis
  → report-agent generates spatial report
  → trace-agent generates traceability document
  
Agent: "Spatial analysis done. Clusters show clear spatial segregation."
```

### Example 4: Image + Expression Multimodal Integration (Xenium Full FOV)

Tested on **Xenium Human Lung 2 FOV** (11,898 cells) and **Human Breast Cancer 2 FOV** (7,020 cells) with full-field-of-view coverage:

```
User:  "Run the full multimodal image + expression pipeline on my Xenium data"
Agent: "Loading morphology and expression data."

  → Level 1 — Cell Morphology: regionprops (area, eccentricity, solidity, etc.)
  → Level 2 — Microenvironment: NND, cell density
  → Level 3 — Tissue Compartments: Otsu thresholding → epithelium/necrosis/stroma/border
  → Merge: 541 genes + 18 image features per cell
  → Strategy 1 — Hexagonal Grid (50µm): aggregate cells into spatial bins
  → Strategy 2 — Shared kNN (k=10): propagate features across spatial graph
  → Strategy 3 — Compartment Scaffold: stratify by tissue histology
  → Leiden Clustering (res=0.3/0.5/1.0) on multimodal PCA
  → Feature Importance: RF classifier, permutation importance, per-cluster signatures

Agent: "Full FOV complete. 20 multimodal clusters found. compartment_boundary_distance
       is the top predictor (45.6% permutation importance). Report ready."
```

## How It Works

```
                          ┌─────────────────────────────┐
                          │       scai-orchestrator       │
                          │  (decides, delegates, reviews) │
                          └──────────┬──────────────────┘
                   ┌──────────────────┼──────────────────┐
                   ▼                  ▼                  ▼
           ┌───────────┐    ┌──────────────┐    ┌──────────────┐
           │ data-agent │    │   qc-agent   │    │ normalize-   │
           │ load data  │───▶│ QC + filter  │───▶│ agent        │
           └───────────┘    └──────────────┘    │ norm + HVG   │
                                                └──────┬───────┘
                                                       ▼
           ┌───────────┐    ┌──────────────┐    ┌──────────────┐
           │ trace-    │◀───│ report-agent │◀───│ cluster-     │
           │ agent     │    │ generate     │    │ agent        │
           │ trace-    │    │ report       │    │ PCA, UMAP,   │
           │ ability   │    │              │    │ Leiden,      │
           │ document  │    │              │    │ markers      │
           └───────────┘    └──────────────┘    └──────────────┘
                                      ▲                  ▲
                                      │                  │
                              ┌───────┴────────┐  ┌──────┴───────┐
                              │ integration-   │  │ spatial-     │
                              │ agent          │  │ agent        │
                              │ WNN, MOFA+     │  │ Squidpy, SD  │
                              └────────────────┘  └──────────────┘
```

## Cross-Session Learning

CellAgent remembers past analyses using Engram. The more you use it, the smarter it gets:

- "Last time we used min_genes=500 for this tissue type"
- "Previous PBMC analysis had doublet issues, let's run scrublet"
- "The UMAP resolution that worked before was 0.8"

## Project Structure

```
CellAgent/
├── opencode.jsonc                # Agent configuration (OpenCode v3+)
├── agents/
│   └── scai-orchestrator/
│       └── INSTRUCTIONS.md       # Orchestrator persona and workflow
├── skills/
│   ├── _shared/pipeline-state.md # Shared state schema
│   ├── data-agent/SKILL.md       # Load h5ad/h5mu/zarr
│   ├── qc-agent/SKILL.md         # Quality control and filtering
│   ├── normalize-agent/SKILL.md  # Normalization and HVG
│   ├── cluster-agent/SKILL.md    # PCA, UMAP, Leiden, markers
│   ├── image-agent/SKILL.md      # Image morphology pipeline (classical CV)
│   ├── integration-agent/SKILL.md# Multimodal integration (muon)
│   ├── spatial-agent/SKILL.md    # Spatial analysis (Squidpy)
│   ├── report-agent/SKILL.md     # Executive report
│   └── trace-agent/SKILL.md      # Traceability document
├── environment.yml               # Conda environment
├── Dockerfile                    # Docker image
└── Singularity.def               # Singularity definition
```

## Roadmap

- **Phase 1** ✅ Unimodal pipeline (data → QC → normalize → cluster → report → trace)
- **Phase 2** ✅ Multimodal integration with muon (WNN, MOFA+)
- **Phase 3** ✅ Spatial analysis with SpatialData and Squidpy
- **Phase 4** ✅ Image morphology + expression integration (hex grid, kNN, compartments, Leiden)
- **Phase 5** 🔄 Feature importance & interpretability (RF, permutation importance, per-cluster signatures)

## License

MIT
