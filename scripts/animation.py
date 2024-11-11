##### (C) Spencer Passmore 2023 
'''
File: generate_interactive.py
Purpose: For generating a nice animation.

Usage:
'python3 animation.py lysozyme_fine_grid'

Notes:
By default, assumes that outputs/batches is in AC4DC/output/__Molecular/, and that this file is in AC4DC/scripts


'''
####
from types import SimpleNamespace 
from interactive_handler import generate_graphs

P = SimpleNamespace()
P.NORMALISE = False
P.ELECTRON_DENSITY = False # if False, use energy density  
P.ANIMATION = True # if true, generate animation rather than interactive figure (i.e. automatic slider movement) 
P.ALSO_MAKE_PLOTS = False # Generate static plots for each simulation using _generate_plots.py
P.SINGLE_FRAME = False # Save a png using plotly .
P.NAMING_MODE = 1  # For legend. 0: full details of sim parameters + sim name | 1: elements in sim| 
P.SCALE_DENSITY_BY_THOUSAND = False # Use cubic nm rather than cubic angstrom for measuring energy density  
P.END_T = 9999  # Put at value to cutoff times early.
P.POINTS = 200

generate_graphs(P)


############### Scratchpad
#python3.9 scripts/generate_interactive.py carbon_classic_static carbon_classic_dynamic
#python3.9 scripts/generate_interactive.py lys_nass_no_S_3 lys_nass_gauss lys_nass_Gd_full_1 