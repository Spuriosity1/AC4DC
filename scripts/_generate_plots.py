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
PLOT_ELEMENT_CHARGE= True #
PLOT_FREE_CONTINUUM = True
PLOT_FREE_SLICES=False
PLOT_ION_RATIOS=False
PLOT_ION_RATIOS_BARS= False
PLOT_ORBITAL_DENSITIES = False #
PLOT_PHOTO_RATES = False
###
COLUMNWIDTH = 3.4975
#COLUMNWIDTH = 2
FIGWIDTH = COLUMNWIDTH#/2
FIGHEIGHT = FIGWIDTH*1/2#*9/16

# FIGWIDTH = COLUMNWIDTH/2
# FIGHEIGHT = FIGWIDTH*9/16

DPI = 800
##
END_T = None#None
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
        make_some_plots(data_folder,molecular_path,label,dname_Figures,PLOT_ELEMENT_CHARGE,PLOT_ION_RATIOS,PLOT_FREE_CONTINUUM,PLOT_FREE_SLICES,PLOT_ION_RATIOS_BARS,PLOT_ORBITAL_DENSITIES,PLOT_PHOTO_RATES)

def make_some_plots(mol_name,sim_output_parent_dir, label,figure_output_dir, tot_charge=False,bound_ionisation=False,free=False,free_slices=False,bound_ionisation_bar=False,orbital_densities_bar=False,photo_rates = False):
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
    load_specific_atoms = None #["C","N","O"] #None  # If plotting free dsitribution, will combine the contributions from those specified here (e.g. if "C","N" then dist_C.csv and dist_N.csv ). If none is specified, will just use the full continuum freeDist.csv.
    if load_specific_atoms is not None:
        label+="_"
        for elem in load_specific_atoms:
            label+=elem

    pl = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY,end_t = END_T,load_specific_atoms=load_specific_atoms)
    num_atoms = len(pl.statedict)
    num_subplots = tot_charge + free + free_slices + bound_ionisation_bar + (bound_ionisation+ orbital_densities_bar+ photo_rates)*num_atoms 
    pl.setup_axes(num_subplots)
    if num_subplots > 1:
        pl.fig.tight_layout()
        pl.fig.subplots_adjust(left=0.12/pl.axs.shape[0], bottom=None, right=None, top=None, wspace=0.2, hspace=None)

    if tot_charge: 
        #NOTE ensure load_specific_atoms is None or does not exclude `atoms` if `atoms` is passed.
        pl.plot_tot_charge(ylim=[None,None],every=1,charge_difference=True,legend_loc="best")  
        #pl.plot_tot_charge(ylim=[None,None],every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O"])  
        #pl.plot_tot_charge(ylim=[0,6],every=1,charge_difference=False,legend_loc="best",atoms=["C","N","O"])  
        #pl.plot_tot_charge(ylim=[0,6],plot_legend=False,every=1,charge_difference=True,legend_loc="best")  
        #pl.plot_tot_charge(ylim=[0,6],plot_legend=False,every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O","S","Gd_fast"])  
        #pl.plot_tot_charge(ylim=[0,1.1],plot_legend=False,every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O","S","Gd_fast"])  
        #pl.plot_tot_charge(ylim=[0,1.1],plot_legend=False,every=1,charge_difference=True,legend_loc="best",atoms=["C","N","O","S"])  
        #pl.plot_tot_charge(ylim=[0,6],plot_legend=False,every=1,charge_difference=False,legend_loc="best",atoms=["C","N","O","S"])  
        #pl.plot_tot_charge(ylim=[0,2.5],plot_legend=False,every=1,charge_difference=True,scale_intensity=0.939)  #TODO automatically set to charge_difference to True if starting with ions...
        #pl.plot_tot_charge(ylim=[0,2.5],xlim=[-15.5,0.5],plot_legend=False,every=1,charge_difference=True)  #TODO automatically set to charge_difference to True if starting with ions...
 

    if bound_ionisation_bar:
        pl.plot_charges_bar("C",show_pulse_profile=False)
        #plt.gcf().set_figwidth(15)        
    if orbital_densities_bar:
        pl.plot_orbitals_bar(atoms=None,atoms_excluded=["N","O"],show_pulse_profile=False,normalise = True)
        #pl.plot_orbitals_bar("Gd_fast",show_pulse_profile=True,orbitals=["3p","4p","5p"])
    if photo_rates:
        pl.plot_photoionisation(atoms=None,show_pulse_profile=True)
    if bound_ionisation:
        pl.plot_all_charges(show_pulse_profile=False,ylim=[0,1])

        #Abdallah
        #pl.plot_all_charges(show_pulse_profile=False,xlim=[-40,0],ylim=[0,1])
        #pl.fig.set_figwidth(6)
        #pl.fig.set_figheight(4)
        # #Royle B
        # pl.plot_charges_royle_style("Al", True)
        # pl.fig.subplots_adjust(left=0.2,bottom=0.18,top=0.95)
        # pl.fig.set_figheight(2.5) # total dimensions, for when other plots turned off...
        # pl.fig.set_figwidth(6)
        # # leonov
        # pl.plot_charges_leonov_style("Si",show_pulse_profile=True,xlim=[-20,20],ylim=[0,1.099])
        # pl.fig.set_figwidth(6.662*0.7)  
        # pl.fig.set_figheight(6*0.7)  
        
    if free:
        #pl.plot_free(log=True,cmin=10**(-7.609),cmax=1e-3,ylim=[10,8000])
        pl.plot_free(log=True,ylog=False,cmin=10**(-8),cmax=10**(-3.609),ylim=[0,8000],keV=True)
        # #Leonov
        # ymax = 9e3
        # pl.plot_free(log=True, cmin=10**(-6.609),cmax = 10**(-2), every=5,mask_below_min=True,cmap='turbo',ymax=ymax,leonov_style=True)
        # pl.fig.set_figwidth(6.662*0.7*1.16548042705)  
        # pl.fig.set_figheight(6*0.7)      
    if free_slices:
        pl.initialise_step_slices_ax()
        from plotter_core import fit_maxwell, maxwell

        cmap = plt.get_cmap("tab10")

        plt.rcParams["font.size"] = 10 # 8

        ### defaults
        #Cutoff energies for fitting MB curves. TODO get from AC4DC
        thermal_cutoff_energies = [2000]*10     
        xmin,xmax = 1,1e4
        
        ###

        # slices = [-15,0,15]
        #TODO get slices from AC4DC or user input           
        # # Hau-Riege (Sanders results)
        # thermal_cutoff_energies = [200, 500, 500, 1000]
        # slices = [-7.5,-5,-2.5,0]         
        # # -7.5 fs Hau-Riege
        # thermal_cutoff_energies = [200]
        # slices = [-7.5] 

        #abdallah
        #slices = [-39,-38,-36,-34,-32,-30,-0.01]         
        # Royle Sect. B
        thermal_cutoff_energies = [500,1000,1250,2000]
        slices = [-90, -50, 0, 50]  
        # # Royle Sect. C
        # thermal_cutoff_energies = [500,500,600,2000]
        # slices = [-15, 0, 15, 30]

        colrs = [cmap(i) for i in range(len(slices))]

        plot_legend = True
        plot_fits = False # Whether to fit MB curves to distribution below thermal cutoff energies.
        plot_those_darn_knots = False
        ####### 
        # Here we can plot MB curves e.g. for fit comparison
        plot_custom_fits = False
        ####
        # # example
        # custom_T = [30.8,75.5,131.8,207.5]
        # custom_n = [0.06*3/2,0.06*3/2,0.06*3/2,0.06*3/2]
        # #Sanders/H-R for -7.5 fs  
        # custom_T = [44.1,31]  
        # custom_n = [0.06*3/2,0.06*3/2] # - (not anything meaningful, just density of MB approximated as 50% of total). 
        # custom_colrs = ['r','b']
        #H-R
        # custom_T = [31,70,125,195]   # [44.1,84.9,135.6,205.8]
        # custom_n = [0.06*3/2]*len(custom_T) # - (not anything meaningful, just density of MB approximated as 50% of total). 
        # custom_colrs = colrs   
        # xmin,xmax = 1,1e4 
        # plot_legend = True 
        #Royle B
        xmin, xmax = 0,2000
        colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]  
        custom_colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]  
        plot_fits = True
        plot_custom_fits = False
        custom_T = [8,33,105,120]  
        custom_n = [0.12,0.05,0.12,0.12] 
        # #Royle C
        # custom_T = [100]  
        # custom_n = [0.155]     
        # colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]          
        # custom_colrs = ['black']
        # plot_fits = False
        # plot_custom_fits = False
        # # Fit our last point in time 
        # if plot_custom_fits == False:
        #     T = pl.plot_fit(slices[-1], thermal_cutoff_energies[-1], normed=True, color=cmap(3), lw=1.5,alpha=0.7)
        # xmin,xmax = 0,1e4
        #######

        #v_anchors = [0.16,0.12,0.08,0.04]
        #v_anchors = [0.2,0.15,0.1,0.05]
        v_anchors = [0.19,0.12,0.05]

        lw = 1.5
        #pl.ax_steps.set_ylim([0.4e-4, 0.4])
        pl.ax_steps.set_ylim([1e-4, 1])
        #pl.ax_steps.set_ylim([1e-4, 1])
        #TODO get from AC4DC
        pl.ax_steps.set_xlim([xmin,xmax]) #Hau-Riege        
        lines = []
        for (t, e, col ) in zip(slices, thermal_cutoff_energies, colrs):
                lines.extend(pl.plot_step(t, normed=True, color = col, lw=lw))
                if plot_fits:
                    T = pl.plot_fit(t, e, normed=True, color=col, lw=lw,alpha=0.7)
        if plot_custom_fits:
            for (T, n, col ) in zip(custom_T, custom_n, custom_colrs):
                pl.ax_steps.plot([0],[0],alpha=0,label=None)
                pl.plot_maxwell(T,n,color = col, lw=lw,alpha=0.7)
        if plot_those_darn_knots:
            ymin,_ = pl.ax_steps.get_ylim()
            pl.ax_steps.set_ylim([ymin*0.67,None])
            pl.plot_the_knots(slices,v_anchors,colrs,padding=0.14)
            #pl.fig_steps.set_figheight(4.8)
            # pl.fig_steps.set_figwidth(7)  
            # for l in lines:
            #     l.set_linewidth(3)          
        if ELECTRON_DENSITY:
            pl.ax_steps.set_ylim([2e-7, 1e-2]) #royle sect. B
            # pl.ax_steps.set_ylim([1e-8, 1e-3]) #royle sect. C
            pl.ax_steps.set_xscale("linear")

        #pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)
        pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
        pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False

        if plot_legend:
            handles, labels = pl.ax_steps.get_legend_handles_labels()
            ncols = 2
            # e.g. for 8 labels (4 time steps), order = [0,2,4,6,1,3,5,7]  , if ncols = 2.
            order = []
            #2 cols
            order = list(range(0,len(labels) - 1,ncols)) + list(range(1,len(labels),ncols)) 
            # 2 rows 4 cols
            #order = list(range(0,len(labels) - 1,ncols)) + list(range(1,len(labels),ncols)) +  list(range(2,len(labels) - 1,ncols)) + list(range(ncols-1,len(labels),ncols)) 
            if len(labels)%2 != 0:
                order.append(len(order))  # shouldnt happen though.
            if plot_fits == False and plot_custom_fits == False:
                ncols = 1
                order = list(range(0,len(labels)))
            pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper center',ncol=ncols)

        name = label.replace('_',' ')
        pl.ax_steps.set_title(name + " - Free-electron distribution")

        #plt.savefig(figure_output_dir + label + fname_HR_style + figures_ext)
    #Abdallah
    # pl.ax_steps.set_xscale("linear")
    # pl.ax_steps.set_xlim([0,1800])
    # pl.ax_steps.set_ylim([0.5e-4,0.5])
    pl.delete_remaining_axes()
    if num_subplots > 2:
        plt.gcf().set_figwidth(FIGWIDTH*pl.axs.shape[0])
        plt.gcf().set_figheight(FIGHEIGHT*pl.axs.shape[1])
    elif num_subplots == 2:
        plt.gcf().set_figheight(FIGWIDTH*pl.axs.shape[0])
    else:
        plt.gcf().set_figwidth(FIGWIDTH)
        plt.gcf().set_figheight(FIGHEIGHT)          
    #plt.tight_layout()
    plt.savefig(figure_output_dir + label + figures_ext,dpi=DPI)
    plt.close()

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option