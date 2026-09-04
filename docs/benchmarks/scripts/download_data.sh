#!/bin/bash
# Download benchmark datasets from PRIDE FTP and quantms-test-datasets
#
# Usage:
#   bash download_data.sh              # Download all datasets
#   bash download_data.sh yeast        # Download yeast (PXD001819) only
#   bash download_data.sh lfq          # Download LFQ mixed-species (PXD028735) only
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
mkdir -p "$DATA_DIR"

PRIDE_FTP="https://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/benchmarks"
QUANTMS_DATASETS="https://raw.githubusercontent.com/bigbio/quantms-test-datasets/quantms"

download_if_missing() {
    local URL="$1"
    local DEST="$2"
    if [ -f "$DEST" ]; then
        echo "  Already exists: $(basename "$DEST")"
        return 0
    fi
    echo "  Downloading: $(basename "$DEST")"
    wget -q --show-progress -O "$DEST" "$URL"
}

download_yeast() {
    echo ""
    echo "=== PXD001819: Yeast + UPS1 spike-in (LFQ, CID, LTQ Orbitrap Velos) ==="
    local DIR="$DATA_DIR/PXD001819"
    mkdir -p "$DIR"

    download_if_missing \
        "$PRIDE_FTP/lfq/LTQOrbitrapVelos/PXD001819/UPS1_5000amol_R1.mzML.gz" \
        "$DIR/UPS1_5000amol_R1.mzML.gz"

    if [ -f "$DIR/UPS1_5000amol_R1.mzML.gz" ] && [ ! -f "$DIR/UPS1_5000amol_R1.mzML" ]; then
        echo "  Decompressing..."
        gunzip -k "$DIR/UPS1_5000amol_R1.mzML.gz"
    fi

    download_if_missing \
        "$QUANTMS_DATASETS/benchmarks/lfq/PXD001819/yeast_ups1_reviewed.fasta" \
        "$DIR/yeast_ups1_reviewed.fasta" \
        2>/dev/null || echo "  FASTA not found in quantms-test-datasets, provide manually"

    # Create a simple LFQ mods file (oxidation M only)
    if [ ! -f "$DATA_DIR/mods_lfq.txt" ]; then
        cat > "$DATA_DIR/mods_lfq.txt" <<'EOF'
NumMods=3
C2H3N1O1,C,fix,any,Carbamidomethyl
O1,M,opt,any,Oxidation
EOF
    fi

    echo "  Done: $DIR"
}

download_lfq() {
    echo ""
    echo "=== PXD028735: Mixed-species ProteoBench (LFQ, HCD, Q Exactive HF) ==="
    local DIR="$DATA_DIR/PXD028735"
    mkdir -p "$DIR"

    echo "  NOTE: PXD028735 mzML files are large (~1GB+)."
    echo "  Download manually from PRIDE if needed:"
    echo "  $PRIDE_FTP/"
    echo ""

    download_if_missing \
        "$QUANTMS_DATASETS/benchmarks/lfq/PXD028735/ProteoBenchFASTA_MixedSpecies_HYE.fasta" \
        "$DIR/ProteoBenchFASTA_MixedSpecies_HYE.fasta" \
        2>/dev/null || echo "  FASTA not found in quantms-test-datasets, provide manually"

    echo "  Done: $DIR"
}

FILTER="${1:-all}"

echo "========== Downloading Benchmark Data =========="
echo "Target: $DATA_DIR"

if [ "$FILTER" = "all" ] || [ "$FILTER" = "yeast" ]; then
    download_yeast
fi

if [ "$FILTER" = "all" ] || [ "$FILTER" = "lfq" ]; then
    download_lfq
fi

echo ""
echo "Download complete. Data directory:"
du -sh "$DATA_DIR"/* 2>/dev/null || echo "(empty)"
