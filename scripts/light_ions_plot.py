import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': 10,
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook

CHARGE_DIFFERENCE = False # Set True if want ionised atoms' trace to start from origin
PLOT_DERIVATIVE = False # Plot the rate of avg charge gain
PLOT_MODE = 1  # 0: plot all charges, 1: plot element total charges # 2: plot orbital charges
YLIM = [0,None]
FIGWIDTH = 3.49751*2/3
FIGHEIGHT = 3.49751/2
#XLIM = [None,None]
XLIM = [None,None]
#YLIM=[0,6]

CUSTOM_LEGEND=("1s: CNO","2s: CNO","2p: CNO","1s: CNO,S,Gd","2s: CNO,S,Gd","2p: CNO,S,Gd")

def main():
    set_highlighted_excepthook()


    # Basic num arguments check
    assert len(sys.argv[1:]) == 2, "Usage: python compare_ion.py <sim_handle_1> <sim_handle_2>"
        
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    valid_folder_names= True
    for i, data_folder in enumerate(sys.argv[1:]):
        if not path.isdir(molecular_path+data_folder):
            valid_folder_names = False
            print("\033[91mInput error\033[0m (argument \033[91m"+str(i)+ "\033[0m): folder name not found.")
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    data_folders = sys.argv[1:]
    label = data_folders[0] + "_"+data_folders[1]
    make_some_plots(data_folders,molecular_path,label,dname_Figures,plot_derivative=PLOT_DERIVATIVE)

def make_some_plots(mol_names,sim_output_parent_dir, label,figure_output_dir,plot_derivative = False):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = "" #.png
    for plot_mode in (PLOT_MODE,):#(0,1):     # 0: plot all charges, 1: plot element total charges

        fig, axs = plt.subplots(3, 3, sharey=True, facecolor='w')

        dashes = ["dashed","solid"]
        atoms = ("C","N","O")#,"S")
        cmap = plt.get_cmap("Dark2")
        assert len(mol_names) <= 2
        # Hacky way to get overlaid plots TODO
        pl = Plotter(mol_names[0],sim_output_parent_dir)
        if plot_mode == 0:
            pl.setup_axes(4)
        if plot_mode == 1 or plot_mode == 2:
            pl.setup_axes(1)
        for m, mol_name in enumerate(mol_names):    
            pl.__init__(mol_name,sim_output_parent_dir)       
            pl.num_plotted = 0 # ƪ（˘へ˘ ƪ）
            if plot_mode==0:
                pl.plot_all_charges(plot_legend=(m==0),linestyle=dashes[m])
            if plot_mode == 1:
                colours = [cmap(i) for i in range(len(atoms))]
                pl.plot_tot_charge(every=1,linestyle=dashes[m],colours = colours,atoms = atoms,plot_legend=(m==0),xlim=XLIM,ylim=YLIM,charge_difference=CHARGE_DIFFERENCE,plot_derivative=PLOT_DERIVATIVE)
            if plot_mode == 2:
                ax = pl.plot_orbitals_charge(every=1,linestyle=dashes[m],atom = "C",plot_legend=False,xlim=XLIM,ylim=YLIM,plot_derivative=PLOT_DERIVATIVE)       
        if plot_mode == 2:
            custom_legend = CUSTOM_LEGEND
            if custom_legend is None: 
                ax.legend(bbox_to_anchor=(1.02, 1),loc='upper left', ncol=1,handlelength=1)  # Top right legend.
            else:
                handles,_ = ax.get_legend_handles_labels()
                handles = list(handles)
                assert len(custom_legend)==len(handles)
                ax.legend(handles,custom_legend,bbox_to_anchor=(1.02, 1),loc='upper left', ncol=1,handlelength=1)             
        plt.gcf().set_figwidth(FIGWIDTH)
        plt.gcf().set_figheight(FIGHEIGHT)
        #plt.gcf().tight_layout()
        #plt.tight_layout()
        if plot_mode == 0:
            qualifier = "_"+ "BoundComp"
        if plot_mode == 1:
            qualifier = "_"+ "ElementComp"
        if plot_mode == 2:
            qualifier = "_"+ "C-OrbitalComp"
        if PLOT_DERIVATIVE:
            qualifier+="-deriv"
        plt.savefig(figure_output_dir + label +qualifier + figures_ext,bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option