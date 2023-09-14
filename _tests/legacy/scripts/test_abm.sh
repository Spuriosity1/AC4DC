#!/bin/bash
# Run from the top-level directory

echo "============================================================"
echo "  Testing ODE method: SHO "
bin/tests/abm_verif 300 1e-6 5 > tmp/abm5_h6.txt
gnuplot -e "prefix='abm5_h6'" tests/scripts/abm_test.gnuplot
bin/tests/abm_verif 300 1e-6 3 > tmp/abm3_h6.txt
gnuplot -e "prefix='abm3_h6'" tests/scripts/abm_test.gnuplot
bin/tests/abm_verif 300 1e-3 4 > tmp/abm4_h3.txt
gnuplot -e "prefix='abm4_h3'" tests/scripts/abm_test.gnuplot
bin/tests/abm_verif 300 1e-4 4 > tmp/abm4_h4.txt
gnuplot -e "prefix='abm4_h4'" tests/scripts/abm_test.gnuplot
bin/tests/abm_verif 300 1e-6 4 > tmp/abm4_h6.txt
gnuplot -e "prefix='abm4_h6'" tests/scripts/abm_test.gnuplot
echo "Output saved to testOutput";
echo "============================================================"
rm tmp/abm*.txt
