#!/bin/bash
set -uo pipefail
AD=$BENCH/astral-data
MZML=LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
FASTA=ProteoBenchFASTA_MixedSpecies_HYE.fasta
RES=$BENCH/repo/bench-astral; mkdir -p $RES
OIMG=ghcr.io/openms/openms-tools-thirdparty:latest
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
echo "################ COMET ASTRAL $(date -Is) ################"
/usr/bin/time -v docker run --rm -v $AD:/data -v $RES:/out $OIMG bash -c '
cd /out && C=/opt/OpenMS/thirdparty/Comet/comet.exe
$C -p >/dev/null 2>&1
P=comet.params.new
sed -i -E "
s#^database_name = .*#database_name = /data/'"$FASTA"'#;
s#^decoy_search = .*#decoy_search = 1#;
s#^num_threads = .*#num_threads = 8#;
s#^peptide_mass_tolerance_upper = .*#peptide_mass_tolerance_upper = 10.0#;
s#^peptide_mass_tolerance_lower = .*#peptide_mass_tolerance_lower = -10.0#;
s#^peptide_mass_units = .*#peptide_mass_units = 2#;
s#^isotope_error = .*#isotope_error = 4#;
s#^search_enzyme_number = .*#search_enzyme_number = 1#;
s#^allowed_missed_cleavage = .*#allowed_missed_cleavage = 2#;
s#^num_enzyme_termini = .*#num_enzyme_termini = 2#;
s#^fragment_bin_tol = .*#fragment_bin_tol = 0.02#;
s#^fragment_bin_offset = .*#fragment_bin_offset = 0.0#;
s#^theoretical_fragment_ions = .*#theoretical_fragment_ions = 0#;
s#^output_percolatorfile = .*#output_percolatorfile = 1#;
s#^output_txtfile = .*#output_txtfile = 0#;
s#^output_pepxmlfile = .*#output_pepxmlfile = 0#;
s#^add_C_cysteine = .*#add_C_cysteine = 57.021464#;
s#^precursor_charge = .*#precursor_charge = 2 4#;
s#^variable_mod01 = .*#variable_mod01 = 15.994915 M 0 3 -1 0 0 0.0#;
s#^variable_mod02 = .*#variable_mod02 = 42.010565 n 0 1 0 0 0 0.0#;
" $P
grep -q "^peptide_length_range" $P && sed -i -E "s#^peptide_length_range = .*#peptide_length_range = 7 40#" $P
echo PARAMS:; grep -E "^database_name|^decoy_search|^peptide_mass_tol|^isotope_error|^variable_mod0[12]|^output_perc|^fragment_bin_tol|^num_enzyme_termini" $P
cp /data/'"$MZML"' /out/comet_in.mzML
$C -P$P /out/comet_in.mzML 2>&1 | tail -8
rm -f /out/comet_in.mzML
' > $RES/comet_astral.log 2>&1
COMET_EXIT=$?
echo "  comet exit=$COMET_EXIT"; grep -E "Elapsed \(wall|Maximum resident" $RES/comet_astral.log | sed "s/^/  /"
[ -f $RES/comet_in.pin ] && mv -f $RES/comet_in.pin $RES/comet_astral.pin
echo "  rows=$(($(wc -l < $RES/comet_astral.pin 2>/dev/null || echo 1)-1))"
docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y --results-psms /r/comet_astral.t.psms --decoy-results-psms /r/comet_astral.d.psms --only-psms=false /r/comet_astral.pin > $RES/comet_astral.perc.log 2>&1
tp=$RES/comet_astral.t.psms
q=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="q-value")print i}' $tp)
pc=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="peptide")print i}' $tp)
rcol=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++)if($i=="proteinIds")print i}' $tp)
ps=$(awk -F"\t" -v q=$q 'NR>1&&$q<=0.01{c++}END{print c+0}' $tp)
pep=$(awk -F"\t" -v q=$q -v p=$pc 'NR>1&&$q<=0.01{s=$p;gsub(/^[A-Z-]\./,"",s);gsub(/\.[A-Z-]$/,"",s);gsub(/\[[^]]*\]/,"",s);print s}' $tp|sort -u|wc -l)
pr=$(awk -F"\t" -v q=$q -v r=$rcol 'NR>1&&$q<=0.01{print $r}' $tp|tr "\t" "\n"|grep -v "^DECOY_\|^rev_\|^XXX_\|^$"|sort -u|wc -l)
echo "  RESULT comet_astral PSMs@1%=$ps peptides@1%=$pep proteins@1%=$pr"
echo "################ COMET_ASTRAL_DONE $(date -Is) ################"
