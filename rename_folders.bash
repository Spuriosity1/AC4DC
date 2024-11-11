TAG=("N" "S" "Fe" "Zn" "Se" "Gd")
for t in ${TAG[@]}; do
for i in output/Sprtn-out/*;
do
 mv -- "$i" "${i/SH2_${t}-/SH2_${t}-SPRT}"
done;
done