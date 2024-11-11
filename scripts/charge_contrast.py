import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    #thaumatin thing
    'axes.titlesize':10,     # fontsize of the axes title
    'axes.labelsize':10,    # fontsize of the x and y labels   
    'ytick.labelsize':10,
    'xtick.labelsize':10,
    'legend.fontsize':10,    
    'lines.linewidth':2,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook


####
ELECTRON_DENSITY = False # Whether to use electron density for free distribution plots. Energy density if False
###
PLOT_ELEMENT_CHARGE= False #
PLOT_FREE_CONTINUUM = False
PLOT_FREE_SLICES=False
PLOT_ION_RATIOS=False
PLOT_ION_RATIOS_BARS= False
PLOT_ORBITAL_DENSITIES = True #
PLOT_PHOTO_RATES = False
###
COLUMNWIDTH = 6


# FIGWIDTH = COLUMNWIDTH#/2
# FIGHEIGHT = FIGWIDTH*1/2#*9/16

FIGWIDTH = COLUMNWIDTH
FIGHEIGHT = FIGWIDTH*9/16

# FIGWIDTH = COLUMNWIDTH/2
# FIGHEIGHT = FIGWIDTH*9/16

DPI = 800
##
END_T = None
##
def main():
    set_highlighted_excepthook()


    # Basic num arguments check
    if  len(sys.argv) < 2:
        print("Usage: python3 generate_plots.py Carbon_1")
        print({"Pass extra arguments to plot for each"})
        exit()
        
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    valid_folder_names= True
    for i, data_folder in enumerate(sys.argv[1:]):
        if not path.isdir(molecular_path+data_folder):
            valid_folder_names = False
            print("\033[91mInput error\033[0m (argument \033[91m"+str(i)+ "\033[0m): folder name not found.")
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    for data_folder in sys.argv[1:]:
        label = data_folder +'_Plt'
        make_plot(data_folder,molecular_path,label,dname_Figures)

def make_plot(mol_name,sim_output_parent_dir, label,figure_output_dir):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = ".png" #.png
    fname_tot_charge = "tot_charge"
    fname_free = "free"
    fname_HR_style = "HR_style"
    fname_bound_dynamics = "bound_dynamics"
    load_specific_atoms = None#["Fe_singleShell","C"] #None #["C","N","O"] #None  # If plotting free dsitribution, will combine the contributions from those specified here (e.g. if "C","N" then dist_C.csv and dist_N.csv ). If none is specified, will just use the full continuum freeDist.csv.
    label += "_chargeContrast"

    plot_ratio = True
    if plot_ratio:
        pl = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY,end_t = END_T,load_specific_atoms=load_specific_atoms)
        num_atoms = len(pl.statedict)
        pl.setup_axes(1)

        empirical_data_folder = path.abspath(path.join(__file__ ,"../nass_charge_contrast_data")) + "/"
        _, extra_artists = pl.plot_charge_contrast_custom_thing("Gd_fast",empirical_data_paths = [empirical_data_folder+"Gd2.csv",empirical_data_folder+"Gd1.csv",])

        plt.gcf().set_figheight(FIGHEIGHT)
        plt.gcf().set_figwidth(FIGWIDTH)
        plt.savefig(figure_output_dir + label + figures_ext,dpi=DPI,bbox_extra_artists = extra_artists, bbox_inches='tight')
    plot_ratio_continuous = False
    if plot_ratio_continuous:
        COLUMNWIDTHTMP = 3.4975


        # FIGWIDTH = COLUMNWIDTH#/2
        # FIGHEIGHT = FIGWIDTH*1/2#*9/16

        FIGWIDTHTMP = COLUMNWIDTHTMP
        FIGHEIGHTTMP = FIGWIDTHTMP*9/16

        pl = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY,end_t = END_T,load_specific_atoms=load_specific_atoms)
        num_atoms = len(pl.statedict)
        pl.setup_axes(1)

        empirical_data_folder = path.abspath(path.join(__file__ ,"../nass_charge_contrast_data")) + "/"
        pl.plot_charge_contrast("Gd_fast",light_element="C",ylim=[0,64/(6*20+7*10+8*10)],ylim_heavy=[0,64])

        plt.gcf().set_figheight(FIGHEIGHTTMP)
        plt.gcf().set_figwidth(FIGWIDTHTMP)
        
        plt.savefig(figure_output_dir + label +"_c" + figures_ext,dpi=DPI,bbox_inches='tight')
    do_plot_pulse_energy = False
    if do_plot_pulse_energy:
        plot_pulse_energy()

        plt.gcf().set_figheight(FIGHEIGHT/2.6)
        plt.gcf().set_figwidth(FIGWIDTH*0.8937)
        plt.tight_layout()
        plt.savefig(figure_output_dir + "pulse_energy" + figures_ext,dpi=DPI)

    plt.close()

def plot_pulse_energy():
    fig,ax = plt.subplots(1,1)
    custom_t = [0,35,37,62,102,112] 
    pulse_energy = [0.33,0.95/2,0.99/2,0.95/2,0.79/2,0.93/2]
    col = "red"
    ax.scatter(custom_t,pulse_energy,marker="x",color=col,label = "Probe",linewidths=1,s=15)
    ax.scatter(custom_t[1:],pulse_energy[1:],marker="+",color=col,label = "Pump",linewidths=1)
    ax.set_ylabel("Pulse energy (mJ)")
    ax.set_xlabel("Probe pulse delay (fs)")
    ax.set_ylim([0.3,0.55])
    ax.legend(loc="best")
    plt.tight_layout()

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option