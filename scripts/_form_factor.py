#%%
%matplotlib widget


import os
os.getcwd() 
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
print(sys.path)


import matplotlib
#matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': 10,
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plotter_core import Plotter
import sys, traceback
import os.path as path
from QoL import set_highlighted_excepthook
import numpy as np


set_highlighted_excepthook()

FIGWIDTH = 3.49751
FIGHEIGHT = FIGWIDTH*2/3
###
ANG_TO_BOHR =0.529177
HIGH_RES = 1
LOW_RES = 4
NUM_TRACES = 5
SHOW_AVG=False
PLOT_DIFF = True
YLIM = [-35,5]
###
handles = ["lys_nass_no_S_3","lys_nass_gauss","lys_nass_Gd_full_1"]
#handles=["lys_nass_no_S_3"]
##
#atoms = pl.atomdict
#atoms = ["N","Gd_fast"]
#atoms = ["C","N","O"]
atoms = ["C"]
#######
dname_Figures = "../../output/_Graphs/plots/"
dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
for handle in handles:
    label = handle +'_' 

    name = handle.replace('_',' ')

    pl = Plotter(handle)
    plt.close()
    photon_energy = 6000
    # q_max is in units of bohr^-1 (atomic units). 
    q_max = 2*np.pi/HIGH_RES*ANG_TO_BOHR
    pl.initialise_form_factor_params(-18,18,q_max,photon_energy,50,100,True)


    
    for atom in atoms:
        # Plot atomic form factor at different times as a function of q
        #pl.plot_ffactor_get_R_AC4DC(atom)   # More accurate form factors from AC4DC

        # Slater form factors. Should be fine for light atoms.
        pl.plot_form_factor(NUM_TRACES,[atom],resolution=True,q_min=2*np.pi/LOW_RES*ANG_TO_BOHR,fig_width=FIGWIDTH,fig_height=FIGHEIGHT,show_average=SHOW_AVG,percentage_change=PLOT_DIFF)
        #plt.title(handle+" - "+atom) 
        plt.tight_layout()
        plt.ylim(YLIM)
        plt.savefig(dname_Figures+atom+"-ff-"+handle+".pdf")
    if len(atoms) > 1:
        pl.plot_form_factor(NUM_TRACES,atoms,resolution=True,q_min=2*np.pi/LOW_RES*ANG_TO_BOHR,fig_width=FIGWIDTH,fig_height=FIGHEIGHT,show_average=SHOW_AVG,percentage_change=PLOT_DIFF)
        plt.title(handle+" - combined") 

    #print(pl.get_A_bar(-10,-7.5,12,12,"C_fast","C_fast",100))

    vmin = 0
    vmax = 1

    #pl.plot_A_map("C_fast","C_fast",vmin=vmin,vmax=vmax)
    #pl.plot_A_map("C_fast","C_fast",vmin=vmin,vmax=vmax,title = r"Naive Foreground coherence $\bar{A}$")
    #pl.plot_R_factor()
    #pl.print_bound_slice(-10)
    #print(" GAP ")
    #pl.print_bound_slice(-7.5)

    print("Done! Remember to stretch!")

# %%
