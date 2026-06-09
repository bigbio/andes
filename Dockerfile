# syntax=docker/dockerfile:1
#
# Andes — the data-driven peptide search engine of the quantms ecosystem.
#
# Multi-stage build: compile the `andes` Rust binary, then ship it on a slim
# Debian runtime with the bundled model store (`resources/ionstat/models.parquet`)
# placed next to the binary so it is self-contained.
#
# This image builds the pure-Rust path, which reads **mzML** and **MGF** spectra
# (no system runtime required). Native Thermo `.raw` (needs a bundled .NET 8
# runtime) and native Bruker timsTOF `.d` are shipped in the per-platform
# release archives instead — see the GitHub Releases page. For `.raw`/`.d` in a
# container, convert to mzML first (e.g. ThermoRawFileParser) or use a release
# binary.
#
# Build:  docker build -t andes .
# Run:    docker run --rm -v "$PWD":/data andes \
#             --spectrum /data/sample.mzML --database /data/proteins.fasta \
#             --output-pin /data/sample.pin

# ---- Stage 1: builder ------------------------------------------------------
FROM rust:1.87-bookworm AS builder

WORKDIR /build

# Copy the full workspace (Cargo.toml, crates/, resources/, etc.).
COPY . .

# Build the release binary for the `andes` package only. Default features =
# mzML + MGF input; no native-format C/.NET dependencies.
RUN --mount=type=cache,target=/build/target \
    --mount=type=cache,target=/usr/local/cargo/registry \
    cargo build --release -p andes \
    && cp target/release/andes /usr/local/bin/andes

# ---- Stage 2: runtime ------------------------------------------------------
FROM debian:bookworm-slim

LABEL org.opencontainers.image.title="andes"
LABEL org.opencontainers.image.description="The data-driven peptide search engine of the quantms ecosystem (clean-room Rust)."
LABEL org.opencontainers.image.source="https://github.com/bigbio/andes"
LABEL org.opencontainers.image.documentation="https://github.com/bigbio/andes"
LABEL org.opencontainers.image.licenses="see LICENSE + NOTICE"
LABEL about.tags="Proteomics"
LABEL maintainer="Yasset Perez-Riverol <ypriverol@gmail.com>"

# The binary opens files and writes PIN/TSV; no extra runtime packages are
# needed for the mzML/MGF path. ca-certificates is included for any TLS use.
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Lay out the install so `current_exe()` finds the model store next to the
# binary: <exe_dir>/resources/ionstat/models.parquet (see bundled_store_path()).
WORKDIR /opt/andes
COPY --from=builder /usr/local/bin/andes /opt/andes/andes
COPY resources/ /opt/andes/resources/
COPY LICENSE NOTICE README.md /opt/andes/

# Expose on PATH. current_exe() canonicalizes the symlink to /opt/andes/andes,
# so the resources lookup resolves correctly.
RUN ln -s /opt/andes/andes /usr/local/bin/andes

# Run as a non-root user (security hardening). A mounted /data must be writable
# by this uid; if the host directory is owned by another user, override with
# `docker run --user "$(id -u):$(id -g)"`.
RUN useradd --create-home --uid 10001 andes && mkdir -p /data && chown andes:andes /data
USER andes

WORKDIR /data
ENTRYPOINT ["andes"]
CMD ["--help"]
