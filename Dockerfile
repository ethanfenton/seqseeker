# ── Stage 1: build ──────────────────────────────────────────────────────────
# Needs gcc + zlib headers to compile the C kmer search binary.
FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install into an isolated venv so it copies cleanly to the runtime stage.
COPY . /src
WORKDIR /src
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir .

# Pre-compile the C binary so the runtime image needs no compiler.
RUN python -c "from seqseeker._compile import get_binary; get_binary(force=True)"

# ── Stage 2: runtime ─────────────────────────────────────────────────────────
# Only zlib runtime (no dev headers, no gcc) — keeps the image small.
FROM python:3.12-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Users bind-mount their FASTQ directories under /data
WORKDIR /data
ENTRYPOINT ["seqseeker"]
