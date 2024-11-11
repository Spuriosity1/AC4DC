set -x
#!/usr/bin/env



HANDLE=("SH_all_light-9_1" "SH_Zn-8_1" "SH_Zn-5_3" "SH_Zn-9_1" "SH_Zn-6_1" "SH_Zn-7_1")
for f in ${HANDLE[@]}; do python3.9 "scripts/_generate_plots.py" "${f}"; done


