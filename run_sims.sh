set -x
#!/usr/bin/env

# Full batch 
#for f in input/batch_lys/lys-{1..150}.mol; do ./ac4dc "$f"; done

### Specific numbers ###
# NUM=(3 5)
# for n in ${NUM[@]}; do ./ac4dc "input/_batches/batch_lys/lys-${n}.mol"; done

## Specific strings ###
#TAG=("Ne" "Ar" "Zn" "Se" "Kr" "Zr" "Xe")
#for t in ${TAG[@]}; do ./ac4dc "input/superheavy/SH_${t}.mol"; done

### Loop all ###
# for f in input/_batches/batch_ES/*; do ./ac4dc "$f"; done

### Loop range ##
#for f in input/_batches/batch_SH_Zn/SH_Zn-{1..10}.mol; do ./ac4dc "$f"; done

#############################


#./ac4dc input/nass/lys_nass_water_solvent
#./ac4dc input/nass/lys_nass_water_solvent_9kev
#./ac4dc input/nass/lys_nass_gd_solvent
./ac4dc_continuum input/nass/lys_nass_gauss

#./ac4dc input/galli/lys_galli_LF_no_Gd
# ./ac4dc input/galli/lys_galli_LF
# ./ac4dc input/galli/lys_galli_HF_no_Gd
./ac4dc input/galli/lys_galli_HF
#./ac4dc_consta input/_batches/batch_ES_C/ES_C-0

# #for f in input/_batches/batch_SH2_N/*; do ./ac4dc "$f"; done
# #for f in input/_batches/batch_SH2_S/*; do ./ac4dc "$f"; done
# #for f in input/_batches/batch_SH2_Fe/*; do ./ac4dc "$f"; done
# #for f in input/_batches/batch_SH2_Zn/*; do ./ac4dc "$f"; done
# for f in input/_batches/batch_SH2_Zn/SH2_Zn-{2..9}.mol; do ./ac4dc "$f"; done
# for f in input/_batches/batch_SH2_Se/*; do ./ac4dc "$f"; done
# for f in input/_batches/batch_SH2_Ag/*; do ./ac4dc "$f"; done
# for f in input/_batches/batch_SH2_Zr/*; do ./ac4dc "$f"; done
# for f in input/_batches/batch_SH2_I/*; do ./ac4dc "$f"; done
# for f in input/_batches/batch_SH2_Gd/*; do ./ac4dc "$f"; done

#for f in input/_batches/batch_SH2_N/SH2_N-{14..53}.mol; do ./ac4dc "$f"; done
#for f in input/_batches/batch_SH2_S/SH2_S-{14..53}.mol; do ./ac4dc "$f"; done
#for f in input/_batches/batch_SH2_Fe/SH2_Fe-{14..53}.mol; do ./ac4dc "$f"; done
#for f in input/_batches/batch_SH2_Se/SH2_Se-{1..53}.mol; do ./ac4dc "$f"; done
#for f in input/_batches/batch_SH2_N/SH2_N-{54..93}.mol; do ./ac4dc "$f"; done
#for f in input/_batches/batch_SH2_S/SH2_S-{71..93}.mol; do ./ac4dc "$f"; done
#for f in input/_batches/batch_SH2_Fe/SH2_Fe-{54..93}.mol; do ./ac4dc "$f"; done


# vvv TO BE DONE LATER


# ###########
# ./ac4dc input/nass/lys_nass_gauss
# ./ac4dc input/nass/lys_nass_gauss_9kev
# ./ac4dc input/nass/lys_nass_light
# ./ac4dc input/nass/lys_nass_light_9kev
# for f in input/_batches/batch_SH_Ag/SH_Ag-{904..908}.mol; do ./ac4dc "$f"; done

# #for f in input/_batches/batch_SH2_Gd/SH2_Gd-{6..14}.mol; do ./ac4dc "$f"; done
# #for f in input/_batches/batch_SH_Xe/SH_Xe-{906..907}.mol; do ./ac4dc "$f"; done
# #for f in input/_batches/batch_SH_Fe/SH_Fe-{40..40}.mol; do ./ac4dc "$f"; done
# # for f in input/_batches/batch_ES_C/ES_C-{90..135}.mol; do ./ac4dc "$f"; done
# # for f in input/_batches/batch_ES_C/ES_C-{67..89}.mol; do ./ac4dc "$f"; done
# # for f in input/_batches/batch_SH_N/SH_N-{29..34}.mol; do ./ac4dc_12_thread "$f"; done
# # NUM=(4)
# # for n in ${NUM[@]}; do ./ac4dc_12_thread "input/_batches/batch_SH_excl_Zn_HF/SH_excl_Zn_HF-${n}.mol"; done


# # NUM=(110 137 136 114 116 118 121 125 128 135)
# # for n in ${NUM[@]}; do ./ac4dc_12_thread "input/_batches/batch_ES_C/ES_C-${n}.mol"; done
# ############

# #for f in input/_batches/batch_ES/*; do ./ac4dc "$f"; done
# #for f in input/_batches/batch_ES_L/*; do ./ac4dc "$f"; done
# # NEED TO DO THE BELOW v
# for f in input/_batches/batch_SH_Zr/SH_Zr-{11..11}.mol; do ./ac4dc "$f"; done
# ./ac4dc input/L-edge/L_control
# ./ac4dc input/L-edge/L_Xe-1
# ./ac4dc input/L-edge/L_Se-1