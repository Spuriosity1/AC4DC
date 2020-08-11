#!/bin/bash
# Run from the top-level directory

echo "============================================================"
echo "  Testing ODE method: SHO "
bin/tests/abm_verif 300 1e-6 5 > tmp/abm5_h6.txt
gnuplot -e "prefix='abm5_h6'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-6 3 > tmp/abm3_h6.txt
gnuplot -e "prefix='abm3_h6'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-3 8 > tmp/abm8_h3.txt
gnuplot -e "prefix='abm8_h3'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-4 8 > tmp/abm8_h4.txt
gnuplot -e "prefix='abm8_h4'" tests/abm_test.gnuplot
bin/tests/abm_verif 300 1e-6 8 > tmp/abm8_h6.txt
gnuplot -e "prefix='abm8_h6'" tests/abm_test.gnuplot
echo "Output saved to testOutput";
echo "============================================================"
rm tmp/abm*.txt
