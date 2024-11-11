##### (C) Spencer Passmore 2023 
'''
File: generate_interactive.py
Purpose: For generation of interactive electron density figure, alongside misc plots.

1. generate interactive + figures for single output
'python3 generate_interactive.py lysozyme_3'

2. Compare interactive figures 
'python3 generate_interactive.py lysozyme_3 lysozyme_4'
If two handles provided dotted line is used for first.  

3. Use a folder containing a "batch" of simulation output folders to generate a folder of graphs (batch folder in /output/__Molecular/)
'python3 generate_interactive.py batch_lys'  
if made with generate_batch.py, batch_lys likely contains files like: lys-1_2, lys-2_2, lys-3_2,... but names do not matter.

TODO ac4dc doesnt auto put results into a batch folder currently.
TODO need to be able to specify folders containing outputs e.g. to compare a few outputs
TODO add -N flag for normalisation (cld also add normalisation option to plotly to replace sliders with normed trace vrsns.)
Notes:
Simulation outputs/batches should be in AC4DC/output/__Molecular/ while this file is in AC4DC/scripts


'''
####
from types import SimpleNamespace 
from interactive_handler import generate_graphs

P = SimpleNamespace()
P.NORMALISE = False
P.ELECTRON_DENSITY = False # if False, use energy density  
P.ALSO_MAKE_PLOTS = False # Generate static plots for each simulation using _generate_plots.py
P.SINGLE_FRAME = False # Save a png using plotly .
P.NAMING_MODE = 0  # For legend. 0: full details of sim parameters + sim name | 1: elements in sim| 
P.SCALE_DENSITY_BY_THOUSAND = False # Use cubic nm rather than cubic angstrom for measuring energy density  
P.END_T = 9999  # Put at value to cutoff times early. Note if multiple handles inputted will cutoff all to earliest end time.
P.POINTS = 70

generate_graphs(P)


############### Scratchpad
#python3.9 scripts/generate_interactive.py carbon_classic_static carbon_classic_dynamic
#python3.9 scripts/generate_interactive.py lys_nass_no_S_3 lys_nass_gauss lys_nass_Gd_full_1 