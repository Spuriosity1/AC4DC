bin/tests/splinecheck 20 10000 2> err.txt > tmp/spline.dat
gnuplot -e "ncols=20" tests/test_bspline.gnuplot
