set terminal png
set output "testOutput/bspline_visualisation.png"
plot for [i=2:(ncols+1)] "tmp/spline.dat" u 1:i w l title "k=".(i-2)
