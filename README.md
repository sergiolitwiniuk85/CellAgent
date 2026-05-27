# SCAI — Single Cell AI Framework

## Installation

```bash
# 1. Clone the repo
git clone https://github.com/YOUR_USERNAME/scai.git
cd scai

# 2. Install Python dependencies
pip install scanpy muon squidpy spatialdata anndata mudata leidenalg

# 3. You're done. Run SCAI:
opencode --agent scai-orchestrator
```

### Prerequisites

- [OpenCode](https://opencode.ai) installed and configured
- Python 3.10+ with scverse packages
- A GitHub PAT for cloning (or public repo)

### Quick start

```bash
opencode --agent scai-orchestrator --project /path/to/your/data
# Then say: "I have a PBMC h5ad, I want QC and clustering"
```

The orchestrator will guide you through the entire analysis step by step.

## How Analyses Persist Across Sessions

Every time you run SCAI, it saves the analysis context, results, and decisions to **engram** (OpenCode's persistent memory).

When you start a new analysis, the orchestrator automatically searches past sessions for relevant context:

- "Last time we used min_genes=500 for this tissue type"
- "Previous PBMC analysis had doublet issues, let's run scrublet"
- "The UMAP resolution that worked before was 0.8"

This builds a **memory of decisions** — the framework gets smarter the more you use it.

### What gets saved

| What | Where | Why |
|------|-------|-----|
| Pipeline decisions | engram `scai/pipeline/{id}/*` | Current session state |
| Parameter choices | engram `scai/learned/{type}` | Reuse good defaults |
| Bugs and gotchas | engram `scai/discoveries/*` | Avoid repeating mistakes |
| Dataset characteristics | engram `scai/datasets/*` | Match patterns across datasets |

## Project Structure

```
scai/
├── opencode.jsonc                # Agent configuration
├── agents/
│   └── scai-orchestrator/
│       └── INSTRUCTIONS.md       # Orchestrator persona
└── skills/
    ├── _shared/pipeline-state.md # Shared state schema
    ├── data-agent/SKILL.md       # Load data
    ├── qc-agent/SKILL.md         # Quality control
    ├── normalize-agent/SKILL.md  # Normalization & HVG
    ├── cluster-agent/SKILL.md    # PCA, UMAP, Leiden, markers
    ├── integration-agent/SKILL.md # Multimodal (muon)
    ├── spatial-agent/SKILL.md    # Spatial (Squidpy)
    └── report-agent/SKILL.md     # Executive report
```

## Skills Reference

| Skill | Trigger | What it does |
|-------|---------|--------------|
| data-agent | "load data", "cargar datos" | Load h5ad/h5mu/zarr, describe dataset |
| qc-agent | "QC", "quality control" | Metrics, violins, filtering, doublets |
| normalize-agent | "normalize", "HVG" | normalize_total, log1p, HVG selection |
| cluster-agent | "cluster", "UMAP", "Leiden" | PCA, neighbors, UMAP, clustering, markers |
| integration-agent | "integrate", "WNN", "muon" | Multimodal integration (RNA+ATAC+protein) |
| spatial-agent | "spatial", "Squidpy" | Spatial analysis, neighbors, enrichment |
| report-agent | "report", "summary" | Executive report with plots and parameters |

## Roadmap

- **Phase 1** ✅ Unimodal pipeline (data → QC → normalize → cluster → report)
- **Phase 2** 🔄 Multimodal integration with muon (WNN, MOFA+)
- **Phase 3** ⏳ Spatial analysis with Squidpy (Xenium, Visium)
- **Phase 4** ⏳ Full multimodal + spatial pipeline with cross-session learning
