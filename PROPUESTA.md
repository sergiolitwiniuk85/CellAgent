# SCAI — Single Cell AI Framework

## Architecture Proposal

### What is SCAI?

SCAI is a single-cell analysis framework built on top of OpenCode agents. It enables **scientists without bioinformatics expertise** to run complex single-cell analyses (unimodal, multimodal, and spatial) by simply describing what they want in natural language.

Instead of memorizing scanpy, muon, squidpy and spatialdata APIs, the scientist describes what they want and SCAI executes it.

### Inspiration

SCAI is inspired by the **SDD (Spec-Driven Development)** pattern from Gentle AI:

```
SDD:   explore → propose → spec → design → tasks → apply → verify
SCAI:  load   → QC      → normalize → integrate → cluster → interpret → report
```

Just as SDD decomposes software development into phases with specialized agents, SCAI decomposes bioinformatics analysis into stages with specialist agents. Each stage transforms the data object and the next agent receives the result.

### Architecture

```
Scientist → "I have CITE-seq data, I want to integrate RNA and protein"
       ↓
[ SCAI Orchestrator ]  — understands the pipeline, decides what's next, speaks English
       ↓  delegates to:
├── data-agent      → Load data (h5ad, h5mu, zarr), describe dataset
├── qc-agent        → QC plots, filtering, diagnostics
├── normalize-agent → Normalization, HVG, scaling, log1p
├── cluster-agent   → PCA, neighbors, UMAP, leiden, marker genes
├── integration-agent → Multimodal integration with muon (WNN, MOFA+)
├── spatial-agent   → Spatial analysis with Squidpy (neighbors, enrichment)
└── report-agent    → Generate executive report with key results
```

### State Flow

Each agent receives and produces a serializable data state:

```
PipelineState {
  data: AnnData | MuData | SpatialData  // the main data object
  modality: "rna" | "atac" | "protein" | "multimodal" | "spatial"
  history: Stage[]                       // what was done so far
  qc_metrics: dict                       // computed metrics
  embeddings: dict                       // PCA / UMAP / LSI
  clusters: dict                         // clustering results
  markers: DataFrame                     // detected marker genes
  report_path: str                       // path to generated report
}
```

### MVP (Phase 1) — Unimodal Pipeline

The minimum viable product covers the standard scanpy pipeline for unimodal RNA-seq data:

1. **Load** h5ad → AnnData
2. **QC**: metrics, violins, filtering
3. **Normalize**: normalize_total, log1p, HVG
4. **Cluster**: PCA, neighbors, UMAP, leiden
5. **Markers**: rank_genes_groups
6. **Report**: executive summary with plots

### Phase 2 — Multimodal (Muon/MuData)

Add multimodal integration with muon:

- Load .h5mu
- Modality-specific preprocessing (ATAC, protein)
- Integration: WNN, MOFA+, concat
- Joint visualization

### Phase 3 — Spatial (SpatialData/Squidpy)

Full spatial analysis:

- Load SpatialData .zarr
- Images, labels, points, shapes
- Squidpy: spatial_neighbors, nhood_enrichment, spatial_scatter
- Rendering with spatialdata-plot

### What the Orchestrator DOES

The orchestrator does NOT execute analysis. The orchestrator:

1. **Listens** to the scientist in natural language
2. **Decides** what stage comes next based on pipeline state
3. **Delegates** to a specialist agent by loading its skill
4. **Reviews** the result and asks the scientist whether to continue or adjust
5. **Maintains** state between stages

It's the equivalent of the `gentle-orchestrator` from SDD, but for bioinformatics.

### Skills

Each skill is a `SKILL.md` file containing:

- **Frontmatter**: name, trigger, description
- **Domain knowledge**: what this stage does, why it matters
- **Libraries and functions**: which APIs to call (scanpy, muon, squidpy)
- **Recommended parameters**: sensible defaults for non-expert scientists
- **Validation**: how to verify the result is correct
- **Interpretation**: how to read the results

### How to Use

```bash
# The scientist opens OpenCode in their project with single-cell data
# They call the SCAI orchestrator
opencode --agent scai-orchestrator

# And says:
"I have a PBMC h5ad, I want QC and clustering"
```

The orchestrator guides the entire process step by step, showing results and asking before proceeding.
