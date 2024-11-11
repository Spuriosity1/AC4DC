##### (C) Spencer Passmore 2023 
'''
File: generate_interactive.py
Purpose: For generation of snapshot of electron energy density function.

Usage: Generate snapshots at specified times:
'python3 generate_snapshots lysozyme_3 -15 0 15'

Notes:
By default, assumes that outputs/batches is in AC4DC/output/__Molecular/, and that this file is in AC4DC/scripts


'''
####
import sys
import os.path as path
from types import SimpleNamespace 
from interactive_handler import generate_graphs
from copy import deepcopy

P = SimpleNamespace()
P.NORMALISE = False
P.ELECTRON_DENSITY = False # if False, use energy density  
P.ANIMATION = False # if true, generate animation rather than interactive figure (i.e. automatic slider movement) 
P.ALSO_MAKE_PLOTS = False # Generate static plots for each simulation using _generate_plots.py
P.SINGLE_FRAME = True # Save a png using plotly .
P.NAMING_MODE = 1  # For legend. 0: full details of sim parameters + sim name | 1: elements in sim| 
P.SCALE_DENSITY_BY_THOUSAND = True # Use cubic nm rather than cubic angstrom for measuring energy density  
P.POINTS = 70
P.SINGLE_FRAME_DICT = dict(
    # In inches
    width = (7.24409436834*2/3),#3.49751*2*3/5,
    height = 7.24409436834/3, #3.49751*2/3
    line_width = 2,
    font_size = 10,
    superscript_font_size = 8,
    
    show_legend = False, 

    ylog = False, # < uses linear scale if False 
    xlog = False, # <

    #times = [-10],#[-10],#[0],#[-10,0], # fs.  Use multiple to plot traces from multiple times (colours will be same between frames however.)
    #y_range = [0,0.0014],#lin:[0,0.0014]/[0,0.009],#log: [-5,-1],      # < If axis is log scale, limits correspond to powers of 10.
    
    #times = [-10],
    #y_range = [0,0.0014], # -10 fs
    times = [0],
    y_range = [0,0.009], # 0 fs
    x_range = [None,7500], # (use [None,None] for default)
)
P.INSET = True #False
P.SUBPLOT_AX_FONT_COLOUR = "#495187"# "blue"
P.SUBPLOT_AX_GRID_COLOUR = "#c4c4eb"
P.INSET_DICT = dict(
    axes_kwargs = dict(
        xaxis2 =dict(
            domain=[0.1, 0.4],
            range = [0,120], #[0,50],
            anchor='y2',
            gridcolor=P.SUBPLOT_AX_GRID_COLOUR, 
            zeroline = False,
            #zerolinecolor = P.SUBPLOT_AX_GRID_COLOUR
        ),
        yaxis2=dict(
            domain=[0.15, 0.8],
            range = [0,None],
            anchor='x2',
            gridcolor=P.SUBPLOT_AX_GRID_COLOUR, 
            zeroline = False,
            tickformat = ".1f",
            tick0=P.SINGLE_FRAME_DICT['y_range'][0],
            dtick = (P.SINGLE_FRAME_DICT['y_range'][1]- P.SINGLE_FRAME_DICT['y_range'][0])/3*2*(1+999*P.SCALE_DENSITY_BY_THOUSAND),          
            #zerolinecolor = P.SUBPLOT_AX_GRID_COLOUR
        ),
    ),
)

if  len(sys.argv) < 3:
    print("Usage: Generate snapshots at some number of times, e.g. t = -10 and t = 0, with: 'python3 scripts/"+path.basename(__file__)+" lysozyme_3 tetrapeptide_1 -10 0'")
    exit()
n=0
for k in sys.argv:
   print(k)
   print(k.strip('-').isnumeric())
   if  k.strip('-').isnumeric():
       break 
   n+=1
   
    
if n >= len(sys.argv):
    print("Times not provided.")
    exit()
for snapshot_t in sys.argv[n:]:
    print("Taking snapshot at t = "+snapshot_t+" fs")
    #P.END_T = float(snapshot_t)  # Put at value to cutoff times early.
    P.SINGLE_FRAME_DICT["times"] = [float(snapshot_t)]  # Put at value to cutoff times early.
    generate_graphs(deepcopy(P),sys_argv = sys.argv[0:n]) # deepcopy because plotly modifies the dicts


############### Scratchpad
#python3.9 scripts/generate_interactive.py carbon_classic_static carbon_classic_dynamic
#python3.9 scripts/generate_interactive.py lys_nass_no_S_3 lys_nass_gauss lys_nass_Gd_full_1 
#python3.9 scripts/generate_snapshot.py lys_nass_no_S_3 lys_nass_gauss lys_nass_Gd_gauss_1 