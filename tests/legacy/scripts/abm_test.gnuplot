set terminal png size 640,480
filename="tmp/".prefix.".txt"
set output "testOutput/".prefix."SHO.png"
plot filename u 1:2 title "Numerical" w dots, filename u 1:3 title "Analytic" w lines
set output "testOutput/stdev_".prefix."SHO.png"
set logscale y 10
plot filename u 1:4 title "deviation +" w dots, filename u 1:(-$4) title "deviation -" w dots
