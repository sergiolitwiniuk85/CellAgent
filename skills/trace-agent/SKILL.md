---
name: trace-agent
description: "Trigger: trazabilidad, trace, reproducibilidad, reproducibility, documento final. Generate a comprehensive traceability document after analysis completion."
license: MIT
metadata:
  author: sergiolitwiniuk
  version: "1.0"
---

# Trace Agent — Documento de Trazabilidad

## Propósito

Generar un documento de trazabilidad completo que registre **cada paso, comando, parámetro y resultado** del análisis single-cell. Esto garantiza reproducibilidad, auditoría y respaldo para publicaciones.

Este documento es el **registro forense** del análisis — no un resumen ejecutivo (eso lo hace report-agent), sino la bitácora completa.

## Cuándo se Ejecuta

**Siempre al final del pipeline**, después de report-agent. El orquestador debe invocarlo como etapa final obligatoria.

## Formato de Salida

```
output/{session_id}/traceability_{project}_{YYYY-MM-DD}.md
```

## Estructura del Documento

```markdown
# Trazabilidad de Análisis Single-Cell

## Metadatos del Análisis

| Campo | Valor |
|-------|-------|
| Proyecto | {project_name} |
| Fecha | {YYYY-MM-DD HH:mm} |
| Archivo original | {data_path} |
| Formato | {format} |
| Sesión SCAI | {session_id} |
| Usuario | {user} |

## Resumen del Dataset

- {n_obs} células × {n_vars} genes
- Modadalidad: {modality}
- Layers: {layers}
- Columnas en obs: {obs_columns}

## Pipeline Ejecutado

Las etapas se listan en orden cronológico, tal como fueron ejecutadas.

---

### Etapa 1: Carga de Datos

**Agente**: data-agent
**Timestamp**: {timestamp}
**Comandos ejecutados**:

\`\`\`python
import scanpy as sc
adata = sc.read_h5ad("path/to/data.h5ad")
\`\`\`

**Parámetros**:
- Ninguno (solo lectura)

**Resultados**:
- Dataset cargado: {n_obs} × {n_vars}
- Archivo guardado en: {current_data_path}

**Decisiones**:
- Formato detectado automáticamente
- No se modificó el archivo original

---

### Etapa N: {stage_name}

**Agente**: {agent_name}
**Timestamp**: {timestamp}

**Comandos ejecutados**:

\`\`\`python
# Comando real ejecutado
scanpy.pp.filter_cells(adata, min_genes=200)
scanpy.pp.filter_cells(adata, max_genes=6000)
\`\`\`

**Parámetros**:
| Parámetro | Valor | Justificación |
|-----------|-------|---------------|
| min_genes | 200 | Filtro estándar para eliminar debris |
| max_genes | 6000 | Eliminar posibles dobletes |

**Resultados**:
- Células antes: 10000
- Células después: 8500 (85% retenidas)
- Plots generados:
  - `output/{session_id}/plots/qc_violin_pre.png`
  - `output/{session_id}/plots/qc_violin_post.png`

**Decisiones**:
- Se usó min_genes=200 porque el dataset es PBMC y ese valor es estándar
- Se rechazó filter_genes porque no era necesario

---

### Tabla Completa de Parámetros

| Etapa | Parámetro | Valor | Rango recomendado |
|-------|-----------|-------|-------------------|
| QC | min_genes | 200 | 100-500 |
| QC | max_genes | 6000 | 5000-10000 |
| QC | max_pct_mito | 20 | 5-25 |
| Normalización | método | normalize_total + log1p | — |
| Normalización | target_sum | None | None o 1e4 |
| Normalización | n_top_genes (HVG) | 2000 | 1000-5000 |
| Clustering | n_pcs | 30 | 15-50 |
| Clustering | resolución Leiden | 0.8 | 0.1-2.0 |
| ... | ... | ... | ... |

### Versiones de Librerías

| Librería | Versión |
|----------|---------|
| scanpy | 1.12.1 |
| muon | 0.1.7 |
| squidpy | 1.8.1 |
| anndata | {version} |
| python | 3.12 |

### Archivos Generados

| Archivo | Path |
|---------|------|
| Datos procesados | output/{session_id}/data_processed.h5ad |
| Reporte ejecutivo | output/{session_id}/report.md |
| PCA variance plot | output/{session_id}/plots/pca_variance.png |
| UMAP clusters | output/{session_id}/plots/umap_clusters.png |
| ... | ... |

### Notas del Analista

{space para que el orquestador/científico agregue interpretaciones o comentarios}

---

*Documento generado automáticamente por SCAI Trace Agent*
*Framework: CellAgent (https://github.com/sergiolitwiniuk85/CellAgent)*
```

## Cómo Construir el Documento

El trace-agent recibe el `PipelineState` completo del orquestador con todo el historial.

Pasos:

1. **Leer el estado**: Recibir `PipelineState` del orquestador
2. **Iterar `history[]`**: Para cada etapa, extraer:
   - Nombre del agente
   - Timestamp
   - Comandos ejecutados (si están registrados)
   - Parámetros usados
   - Resultados (métricas numéricas, paths de plots)
   - Decisiones tomadas
3. **Recopilar metadatos**: Versiones de librerías (usar `pip list` o `scanpy.__version__`)
4. **Generar markdown**: Usar la estructura exacta de arriba
5. **Guardar archivo**: En `output/{session_id}/traceability_{project}_{date}.md`

## Recolección de Información

El trace-agent debe ejecutar estos comandos para obtener metadatos del entorno:

```python
import subprocess
import scanpy as sc
import muon as mu
import squidpy as sq
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

libraries = {
    "scanpy": sc.__version__,
    "muon": mu.__version__,
    "squidpy": sq.__version__,
    "anndata": ad.__version__,
    "pandas": pd.__version__,
    "numpy": np.__version__,
    "matplotlib": matplotlib.__version__,
    "seaborn": sns.__version__,
}
```

## Output

```python
{
    "stage": "trace",
    "status": "completed",
    "trace_path": "output/{session_id}/traceability_{project}_{date}.md",
    "n_stages": len(history),
    "summary": f"Documento de trazabilidad generado con {len(history)} etapas"
}
```

## Hard Rules

- **No inventar comandos**. Si el agente no registró qué comandos ejecutó, poner "Comandos: no registrados"
- **No inventar justificaciones**. Si no se registró por qué se eligió un parámetro, poner "No documentado"
- **No modificar los datos**. El trace-agent es read-only
- **Incluir SIEMPRE** la tabla de parámetros completa y las versiones de librerías
- Si el pipeline tiene menos de 2 etapas, documentar como "análisis parcial"
- El documento debe ser **autocontenido**: cualquiera que lo lea debe poder reproducir el análisis
