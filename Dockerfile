# CellAgent — reproducible environment
# Build:  docker build -t cellagent .
# Run:    docker run -it --rm -v $(pwd):/data cellagent

FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc gfortran libopenblas-dev liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Install scverse stack via pip (clean environment, no conflicts)
RUN pip install --no-cache-dir \
    scanpy==1.12.1 \
    muon==0.1.7 \
    squidpy==1.8.1 \
    spatialdata==0.7.3 \
    spatialdata-plot==0.4.0 \
    mudata==0.3.8 \
    leidenalg==0.12.0

# Copy CellAgent framework
WORKDIR /cellagent
COPY . .

# Default: launch bash for interactive use
CMD ["/bin/bash"]
