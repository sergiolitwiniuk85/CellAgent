# CellAgent — scverse base image
# Build:   docker build -t ghcr.io/sergiolitwiniuk85/cellagent:latest .
# Run:     docker run -it --rm -v $(pwd):/data -v /path/to/CellAgent:/cellagent cellagent
#
# Mount your data to /data and CellAgent to /cellagent.
# Inside the container:
#   cd /cellagent && opencode --agent scai-orchestrator --project /data

FROM python:3.12-slim

LABEL org.opencontainers.image.source="https://github.com/sergiolitwiniuk85/CellAgent"
LABEL org.opencontainers.image.description="CellAgent scverse environment - scanpy, muon, squidpy, spatialdata"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc gfortran libopenblas-dev liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Install pinned scverse stack (same versions as environment.yml)
RUN pip install --no-cache-dir \
    scanpy==1.12.1 \
    muon==0.1.7 \
    squidpy==1.8.1 \
    spatialdata==0.7.3 \
    spatialdata-plot==0.3.4 \
    mudata==0.3.8 \
    leidenalg==0.12.0

WORKDIR /data
CMD ["/bin/bash"]
