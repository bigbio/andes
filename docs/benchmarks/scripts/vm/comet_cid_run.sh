#!/usr/bin/env bash
# Run Comet in OpenMS docker; writes Percolator PIN to $4.
# Usage: comet_cid_run.sh <tmt|std> <mzML> <fasta> <out.pin>
set -euo pipefail
MODE="${1:?tmt|std}"
MZML="${2:?mzML}"
FASTA="${3:?fasta}"
OUTPIN="${4:?out.pin}"
OIMG="${OIMG:-ghcr.io/openms/openms-tools-thirdparty:latest}"

WORKDIR="$(mktemp -d /tmp/comet-cid.XXXXXX)"
trap 'rm -rf "$WORKDIR"' EXIT
mkdir -p "$(dirname "$OUTPIN")"
OUTNAME="$(basename "$OUTPIN")"

if [ "$MODE" = tmt ]; then
  PREC_UP=20.0; PREC_LO=-20.0; ISO=4
  TMT_EXTRA='
s#^add_K_lysine = .*#add_K_lysine = 229.162932#;
s#^add_Nterm_peptide = .*#add_Nterm_peptide = 229.162932#;'
else
  PREC_UP=5.0; PREC_LO=-5.0; ISO=2
  TMT_EXTRA='
s#^add_K_lysine = .*#add_K_lysine = 0.0#;
s#^add_Nterm_peptide = .*#add_Nterm_peptide = 0.0#;'
fi

docker run --rm \
  -v "$WORKDIR":/data \
  -v "$(dirname "$OUTPIN")":/out \
  "$OIMG" bash -c "
set -euo pipefail
cd /out
C=/opt/OpenMS/thirdparty/Comet/comet.exe
\$C -p >/dev/null 2>&1
P=comet.params.new
sed -i -E \"
s#^database_name = .*#database_name = /data/db.fasta#;
s#^decoy_search = .*#decoy_search = 1#;
s#^num_threads = .*#num_threads = 8#;
s#^peptide_mass_tolerance_upper = .*#peptide_mass_tolerance_upper = ${PREC_UP}#;
s#^peptide_mass_tolerance_lower = .*#peptide_mass_tolerance_lower = ${PREC_LO}#;
s#^peptide_mass_units = .*#peptide_mass_units = 2#;
s#^isotope_error = .*#isotope_error = ${ISO}#;
s#^search_enzyme_number = .*#search_enzyme_number = 1#;
s#^allowed_missed_cleavage = .*#allowed_missed_cleavage = 2#;
s#^num_enzyme_termini = .*#num_enzyme_termini = 2#;
s#^fragment_bin_tol = .*#fragment_bin_tol = 1.0005#;
s#^fragment_bin_offset = .*#fragment_bin_offset = 0.4#;
s#^theoretical_fragment_ions = .*#theoretical_fragment_ions = 1#;
s#^output_percolatorfile = .*#output_percolatorfile = 1#;
s#^output_txtfile = .*#output_txtfile = 0#;
s#^output_pepxmlfile = .*#output_pepxmlfile = 0#;
s#^add_C_cysteine = .*#add_C_cysteine = 57.021464#;
${TMT_EXTRA}
s#^precursor_charge = .*#precursor_charge = 2 4#;
s#^variable_mod01 = .*#variable_mod01 = 15.994915 M 0 3 -1 0 0 0.0#;
s#^variable_mod02 = .*#variable_mod02 = 0.0 X 0 3 -1 0 0 0.0#;
\" \"\$P\"
grep -q '^peptide_length_range' \"\$P\" && sed -i -E 's#^peptide_length_range = .*#peptide_length_range = 6 40#' \"\$P\"
cp '$MZML' /out/comet_in.mzML
cp '$FASTA' /data/db.fasta
\$C -P\"\$P\" /out/comet_in.mzML
rm -f /out/comet_in.mzML
[ -f /out/comet_in.pin ] && mv -f /out/comet_in.pin /out/${OUTNAME}
"

if [ -f "$OUTPIN" ]; then
  echo "comet OK rows=$(($(wc -l < "$OUTPIN") - 1)) -> $OUTPIN"
else
  echo "comet FAIL no PIN at $OUTPIN" >&2
  exit 1
fi
