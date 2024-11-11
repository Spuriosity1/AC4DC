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
    'lines.linewidth':1,
    'lines.markersize':4,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook
import copy
from scipy.interpolate import InterpolatedUnivariateSpline,SmoothBivariateSpline
import math
import json
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/scattering')
from scatter import XFEL,Crystal,stylin
from damage_landscape import multi_damage
from core_functions import get_sim_params, get_pdb_path
from core_variables import *
from mergedeep import merge
import imaging_params


FIGWIDTH_FACTOR = 1
FIGWIDTH = 3.49751 *FIGWIDTH_FACTOR
FIGHEIGHT = FIGWIDTH*3/4
FIT = True
IGNORE_MISSING_SIMULATIONS = True  # If False, throws an error if not all specified simulations in the batch(es) are found.

matplotlib.rcParams.update({
    'figure.subplot.left':0.16/FIGWIDTH_FACTOR,
    'figure.subplot.right':1-0.02/FIGWIDTH_FACTOR-0.02,
    'figure.subplot.bottom': 0.16/FIGWIDTH_FACTOR,
    'figure.subplot.top':1-0.07/(FIGWIDTH_FACTOR*FIGHEIGHT/FIGWIDTH*4/3),
})



####### User parameters #########
INDEP_VARIABLE = 1 # | 0: energy of photons |1: Energy separation from edge given in edge_dict |2: Artificial electron source energy
DEP_VARIABLE = 0 # | 0: mean carbon charge, either integrated over intensity or at end of pulse ) | 1: R factors (performs scattering simulations for each) |
BATCH = 0 #| 0: Real elements | 1: Electron source. 

## Subcategories of dependent variables
INTENSITY_AVERAGED = True # Used if DEP_VARIABLE = 0 (average charge). | False: Average carbon charge at the simulation's termination | True: Average charge throughout the pulse, weighted by the pulse profile -- NOT the same as the mean "observed" charge by elastic scattering, as does not account for loss of scattering power. 
SCATTERING_TARGET = 0 # Used if DEP_VARIABLE = 1.  | 0: unit lysozyme (light atoms) | 1: 3x3x3 lysozyme (light atoms no solvent)|

## Graphical
LEGEND = True
LABEL_TIMES = False # Set True to graphically check simulation end times are all aligned
DPI = 800
MAX_TICKS = 5

## Numerical

NORMALISING_STEM = None#"SH2_N_1e12" #"SH_N" # None  #Stem (str) to normalise traces by. (s.t. normalised trace becomes horizontal line). If None, does not normalise.
SUBTRACT_NORMALISING_STEM = True # If false, divide
#NORMALISING_STEM = "SH_N"
ylim=[None,None]
CUSTOM_YLABEL = None#"$R_{dmg}^{CNO,Zn}-R_{dmg}^{CNO}$"
#ylim=[None,2.3-0.7]
#CUSTOM_YLABEL = None
#ylim=[None,2.3]
L_TYPE = "three" # Which L edge to use. "one","two","three"; anything else plots all and uses separation from L1 edge.

# #### TEMPORARY
# #mode_A
# A = False
# if A:
#     ylim=[0,0.025]
#     INDEP_VARIABLE= 0
#     BATCH = 0
#     NORMALISING_STEM = None
# #mode_B
# else:
#     ylim=[None,None]
#     INDEP_VARIABLE = 2
#     BATCH = 1
#     NORMALISING_STEM = None
##########

xlim = [6,18]
#xlim = [9.7,19.7]
if INDEP_VARIABLE is EDGE_SEPARATION:
    xlim[0] = -1.5
    xlim[1] = 18#8.5
if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
    #xlim = [None,None]
    xlim = [0,15]




if BATCH is REAL_ELEMENTS:
    assert(INDEP_VARIABLE != ELECTRON_SOURCE_ENERGY)
    # Lists of ranges of simulations to plot
    # Batch stems and index range (inclusive). currently assuming form of key"-"+n+"_1", where n is a number in range of stem[key]
    # e.g. item "Carbon":[2,4] means to plot simulations corresponding to outputs Carbon-2, Carbon-3, Carbon-4, mol files. Last simulation output is used in case of duplicates (highest run number).
    HF_15fs = { # 10^12 ph um^-2 fluence, 10 fs
            "SH_N":[[1,10],[901,903]], #[700,704],
            "SH_S":[[0,10],[901,903]],  
            "SH_Fe":[[0,10],[901,904]],#,[20,40]],
            "SH_Zn":[[0,10],[901,903]],
            "SH_Se":[[0,13],],
            "SH_Zr":[[0,10],[901,902]],
            #-----L------
            "SH_Ag":[[0,8],[901,904]],
            "SH_Xe":[[0,9],[901,907]],
            #"SH2_Gd":[[1,14]],
            #------other------
            #"SH_excl_Zn":[1,9],
            #"L_Gd":[1,9],
    }
    LF_15fs = { 
            "SH2_N":[[1,13]],
            "SH2_S":[[1,13]],  
            "SH2_Fe":[[1,13]],
            "SH2_Zn":[[1,13]],
            "SH2_Se":[[1,13]],
            "SH2_Zr":[[1,13]],
            #-----L------
            "SH2_Ag":[[1,13]],
            "SH2_Xe":[[1,13]],
            "SH2_Gd":[[15,29]],
            #------other------
            #"SH_excl_Zn":[1,9],
            #"L_Gd":[1,9],
    }
    HF_15fs_SH2 = { 
            #"SH2_N":[[24,33]],
            #"SH2_S":[[24,33]], 
            "SH2_Fe":[[24,33]],
            "SH_Fe":[[0,10],[901,904]],
            #"SH2_Zn":[[24,33]],
            "SH2_Se":[[24,33]],
            "SH_Se":[[0,13],],
            #"SH2_Zr":[[24,33]],
            #-----L------
            #"SH2_Ag":[[24,33]],
            #"SH2_Xe":[[24,33]],
    }
    Zn_N_comp = {
        #"SH2_N_1e12":[[24,33],[9919,9927]],
        "SH_N":[[1,10],[901,903]], #[700,704],
        "SH_Zn":[[0,10],[901,903]],

        #"SH2_Zn_1e12":[[9939,9957]],
    }
    HF_10fs = {}
    LF_10fs = {} 
    stem_dict = HF_15fs

same_targets_dict = dict(
    SH2_N=["SH2_N_1e12",], 
    SH2_Zn=["SH2_Zn_1e12",]
    ) 


if BATCH is ELECTRON_SOURCE:
    assert(INDEP_VARIABLE == ELECTRON_SOURCE_ENERGY)
    ##
    #photon_energy = 10
    fluence = 1
    width = 15
    source_fraction = 3
    source_duration = 1
    ##
    # Manual mapping... TODO just pop outputs in folder or something.
    def batch_mapping(fluence,width,source_fraction,source_duration):
        params = fluence,width,source_fraction,source_duration
        return_dict = None
        if params == (1,15,1,1):
            return_dict = dict(ES_C = ([1,22],[0,0]))
        elif params == (1,15,3,1):
            return_dict =  dict(ES_C = ([23,44],[0,0],[902,902]))
        elif params == (1,30,1,1):
            return_dict = dict(ES_C =([45,66],))
        elif params == (100,15,3,1):
            return_dict = dict(ES_C=([90,113],))
        elif return_dict is None:
            raise KeyError("Invalid/missing batch params")
        return return_dict
    
    stem_dict = batch_mapping(fluence,width,source_fraction,source_duration)
        
# Ionisation edges of dopants. Highest absorption edge specified in first element in tuple is plotted if in energy range. Second element in the tuple is label for said edge.
edge_dict = {
    "SH_N": (0.3,"K"),
    "SH2_N": (0.3,"K"),
    "SH_S": (1,"K"),
    "SH2_S": (1,"K"),
    "SH_Fe":(7.1,"K"),  
    "SH2_Fe":(7.1,"K"),  
    "SH_Zn":(9.7,"K"),
    "SH2_Zn":(9.7,"K"),
    "SH_excl_Zn":(9.7,"K"),
    "SH_Se":(12.7,"K"),
    "SH2_Se":(12.7,"K"),
    "SH_Zr":(17.98,"K"),  # Slight offset for graphical purposes.
    "SH2_Zr":(17.98,"K"),  # Slight offset for graphical purposes.
    "SH_Ag":([3.3,3.5,3.8],["L_{3}","L_{2}","L_{1}"]),
    "SH2_Ag":([3.3,3.5,3.8],["L_{3}","L_{2}","L_{1}"]),
    "SH_Xe":([4.8,5.1,5.5],["L_{3}","L_{2}","L_{1}"]),
    "SH2_Xe":([4.8,5.1,5.5],["L_{3}","L_{2}","L_{1}"]),
    #"SH2_Gd":([7.2,7.9,8.4],["L_{3}","L_{2}","L_{1}"]),  
    "SH2_Gd":(7.4,["L_{2}",]),  
    "L_Gd":(7.4,["L_{2}",]),
    "ES":(0.3,"K"),
    "ES_L":(0.3,"K"),
    "ES_C":(0.3,"K"),
}
for elem in edge_dict.keys():
    if type(edge_dict[elem][0]) in [float,int]:
        edge_dict[elem] = ([edge_dict[elem][0],],edge_dict[elem][1]) 
         
# Initial charges of dopant
ground_charge_dict = {
    "SH_N": 0,
    "SH2_N": 0,
    "SH_S": 0,
    "SH2_S": 0,
    "SH_Fe":2,  
    "SH2_Fe":2,  
    "SH_Zn": 2,
    "SH2_Zn": 2,
    "SH_excl_Zn": 2,
    "SH_Se":0,
    "SH2_Se":0,
    "SH_Zr":2,
    "SH2_Zr":2,
    "SH_Ag":1,
    "SH2_Ag":1,
    "SH_Xe":8,
    "SH2_Xe":8,
    "SH2_Gd":3,
    
    "L_Gd":11,

    "ES":0,
    "ES_L":0,
    "ES_C":0,
}
col_dict = {
    "SH_N": 0,
    "SH2_N": 0,
    "SH_S": 1,
    "SH2_S": 1,
    "SH_Fe":7,  
    "SH2_Fe":7,  
    "SH_Zn": 3,
    "SH2_Zn": 3,
    "SH_excl_Zn": 3,
    "SH_Se":4,
    "SH2_Se":4,
    "SH_Zr":5,
    "SH2_Zr":5,
    "SH_Ag":8,
    "SH2_Ag":8,
    "SH_Xe":6,
    "SH2_Xe":6,
    "SH2_Gd":2,    
    
    "L_Gd":2,    

    "ES":0,
    "ES_L":1,
    "ES_C":2,
}

for k, v in same_targets_dict.items():
    for stem in v:
        edge_dict[stem] = edge_dict[k]
        col_dict[stem] = col_dict[k]
        ground_charge_dict[stem] = ground_charge_dict[k]



################
## Constants

SCATTERING_TARGET_DICT = {0: (imaging_params.goldilocks_dict_unit,"-unit"),1: (imaging_params.goldilocks_dict_3x3x3,"3x3x3"),99:(imaging_params.goldilocks_dict_unit,"-unit")}
#
SCATTER_DIR = path.abspath(path.join(__file__ ,"../")) +"/scattering/"
PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"lys") 
MOLECULAR_PATH = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/" # directory of damage sim output folders
GOLDI_RESULTS_PATH = path.abspath(path.join(__file__ ,"../.goldilocks_saves")) + "/" 



set_highlighted_excepthook()
def main():       


    assert all([k in edge_dict.keys() and k in ground_charge_dict.keys() and k in col_dict for k in stem_dict.keys()])
    # Plot tag
    plot_type = {AVERAGE_CHARGE:'avg_charge_plot',R_FACTOR:'R_dmg'}      

    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    
    batches = {}
    # Get folders names in molecular_path
    all_outputs = os.listdir(MOLECULAR_PATH)
    valid_folder_names= True
    for key,val_range in stem_dict.items():
        data_folders = []
        sim_tags = []
        for val in val_range:
        # Get folder name, excluding run tag
            sim_tags += [*range(val[0],val[1]+1)]
        for n in sim_tags:
            data_folders.append(key+"-"+str(n)) 
            # Add run tag "_"+n corresponding to latest run (n = highest "R" for "stem-n_R")
        for i, handle in enumerate(data_folders):
            matches = [match for match in all_outputs if match.startswith(handle+"_")]
            run_nums = [int(R.split("_")[-1]) for R in matches if R.split("_")[-1].isdigit()]            
            # no folders - > not valid 
            if len(run_nums) == 0: 
                valid_folder_names = False
                error_string_prefix = "\033[91mInput error\033[0m a"
                if IGNORE_MISSING_SIMULATIONS:
                    error_string_prefix = "A"
                print(error_string_prefix+" run corresponding to \033[91m"+handle+"\033[0m was not found in"+MOLECULAR_PATH)
                if IGNORE_MISSING_SIMULATIONS:
                    data_folders[i] = None
                    continue
                break
            data_folders[i] = handle + "_"+str(max(run_nums))
        batches[key] = [f for f in data_folders if f is not None]
    if not IGNORE_MISSING_SIMULATIONS:
        assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder. See input error message(s)"
    fig_title = "".join([stem.split("-")[0].split("_")[-1] for stem in batches.keys()])
    label = fig_title+"_"+plot_type[DEP_VARIABLE]
    plot(batches,label,dname_Figures,mode=DEP_VARIABLE)

def plot(batches,label,figure_output_dir,mode = 0):
    ylabel = {AVERAGE_CHARGE:"Averaged carbon charge",R_FACTOR:"$R_{dmg}$"}
    '''
    Arguments:
    '''    
    ############
    # File/directory names
    #######  
    figures_ext = ".png" #.png
    fname_charge_conservation = "charge_conservation"
    fname_free = "free"
    fname_HR_style = "HR_style"
    fname_bound_dynamics = "bound_dynamics"

    fig, ax = plt.subplots(figsize=(FIGWIDTH,FIGHEIGHT))
    cmap = plt.get_cmap('tab10')

    output_tag = ""
    if DEP_VARIABLE is R_FACTOR:
        _,output_tag = SCATTERING_TARGET_DICT[SCATTERING_TARGET]
    if NORMALISING_STEM not in [None,False]:
        output_tag += "-normed"
    if INDEP_VARIABLE is EDGE_SEPARATION:
        output_tag +="-PhEl"  # X-ray photoelectron = XPhEl...
    if BATCH is ELECTRON_SOURCE:
        output_tag += str(fluence)+"-"+str(width)+"-"+str(source_fraction)+"-"+str(source_duration)
                 
    if NORMALISING_STEM not in [None,False]:
        norm = []
        mol_names = batches[NORMALISING_STEM]
        for mol_name in mol_names:
            norm.append(get_data_point(ax,NORMALISING_STEM,mol_name,mode))

    for stem, mol_names in batches.items():
        print("Plotting "+stem+" trace.")
        c = col_dict[stem]
        dopant = stem.split("-")[0].split("_")[-1]
        X = []
        Y = []
        times = []  # last times of sims, for checking everything lines up.
        for mol_name in mol_names:
            if mol_name == NORMALISING_STEM:
                for x,y,t in norm:
                    if x is not None:
                        X.append(x)
                        Y.append(y)  
                        times.append(t)                  
          
        for mol_name in mol_names:
            if mol_name==NORMALISING_STEM:
                continue
            x,y,t = get_data_point(ax,stem,mol_name,mode)
            if x is not None: # Ignore any sourceless data points.
                X.append(x)
                Y.append(y)  
                times.append(t)   
            elif INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:        
                    # Plot horizontal dashed line corresponding to damage with no artificial source.
                    ax.axhline(y=y,ls=(5,(10,3)),color=cmap(0),label="\,No source")          
                    
        print("Energies:",X)
        if mode is AVERAGE_CHARGE:
            print("Charges:",Y)
        if mode is R_FACTOR:
            print("R factors:",Y)
        if NORMALISING_STEM not in [None,False]:
            #assert(INDEP_VARIABLE is not EDGE_SEPARATION) #TODO get working with edge separation
            for i, energy in enumerate(X): 
                found_norm = False
                for elem in norm:
                    if elem[0] == energy:
                        found_norm = True
                        if SUBTRACT_NORMALISING_STEM:
                            Y[i]-=elem[1]
                        else:
                            Y[i]/=elem[1]
                            print(stem)
                            print("ASDASD")
                            print(Y[i])
                        break
                assert(found_norm)
                #     del(X[i])
                #     del(Y[i])
                #     del(energies[i])
                #     del(times[i])

        if BATCH is REAL_ELEMENTS:
            _label = dopant 
            if ground_charge_dict[stem]>1: 
                _label+="$^{"+str(ground_charge_dict[stem])+"+}$"
            elif ground_charge_dict[stem]==1:
                _label+="$^{+}$"
        if BATCH is ELECTRON_SOURCE:
            _label = "\,Source"
        print(mol_name)
        if stem is NORMALISING_STEM:
            continue

        ax.scatter(X,Y,label=_label,color = cmap(c))
        
        k =2 # Spline order
        if FIT and len(X)>=k+1:
            ordered_dat = sorted(zip(X,Y))
            # Split dataset with ionisation edge of dopant
            split_idxes = [0]
            if INDEP_VARIABLE is PHOTON_ENERGY:
                split_idxes = [np.searchsorted(np.array([x[0] for x in ordered_dat]),edge) for edge in edge_dict[stem][0]]
            if INDEP_VARIABLE is EDGE_SEPARATION:
                done = False
                if type (edge_dict[stem][0]) not in [int,float]:
                    # Assume L-edge
                    if L_TYPE in ["one","two","three"]:
                        split_idxes = [np.searchsorted(np.array([x[0] for x in ordered_dat]),- 0.000000001),]
                        if stem == "SH_Xe": split_idxes = [3,4] # Xenon temporary thing.
                        done = True
                # default case
                if not done:
                    split_idxes = np.searchsorted(np.array([x[0] for x in ordered_dat]),[edge-max(edge_dict[stem][0]) - 0.000000001 for edge in edge_dict[stem][0]])
            splines = []
            finer_x = []
            if not [split_idx in (0, len(X)) for split_idx in split_idxes]:
                # No split in fit
                splines.append(InterpolatedUnivariateSpline(*zip(*ordered_dat),k=k))
                finer_x.append(np.linspace(np.min(X),np.max(X),endpoint=True))
            else:
                # Split in fit                   START POINTS             END POINTS
                for start_idx,end_idx  in zip([0,max(split_idxes)],[min(split_idxes)-1,len(X)-1]):
                    if end_idx+1 - start_idx>=k+1:
                        splines.append(InterpolatedUnivariateSpline(*zip(*ordered_dat[start_idx:end_idx+1]),k=k))
                        finer_x.append(np.linspace(ordered_dat[start_idx][0],ordered_dat[end_idx][0],endpoint=True))
            for sp, lil_x in zip(splines,finer_x):
                ax.plot(lil_x, sp(lil_x),color = cmap(c)) 

        if LABEL_TIMES:
            for i, txt in enumerate(times):
                ax.annotate(txt, (X[i],Y[i]))
        _ylab = ylabel[mode]
        if NORMALISING_STEM not in [None,False]:
            _ylab += " relative difference" # TODO better clarity
        if CUSTOM_YLABEL is not None:
            _ylab = CUSTOM_YLABEL
        ax.set_ylabel(_ylab)
        if INDEP_VARIABLE is PHOTON_ENERGY:
            ax.set_xlabel("Photon energy (keV)")
        if INDEP_VARIABLE is EDGE_SEPARATION:
            ax.set_xlabel("Lowest photoelectron emission energy (keV)")  
            #ax.set_xlabel("LEIS photoelectron energy (keV)")  
        if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
            ax.set_xlabel("Source electron energy (keV)")  # TODO define better
    
    if xlim[0] is None:
        xlim[0] = np.min(X)
    if xlim[1] is None:
        xlim[1] = np.max(X)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)           
    ax.yaxis.set_major_locator(plt.MaxNLocator(MAX_TICKS))


    if INDEP_VARIABLE is PHOTON_ENERGY: 
        for stem, mol_names in batches.items():
            c = col_dict[stem]
            dopant = stem.split("-")[0].split("_")[-1]
            i = 0

            for lin, text in zip(edge_dict[stem][0],edge_dict[stem][1]):
                i+=1
                if L_TYPE == "one":
                    if type(edge_dict[stem][1]) in (list,tuple) and  i!=len(edge_dict[stem][1]):
                        continue  
                if L_TYPE == "two":
                    if type(edge_dict[stem][1]) in (list,tuple) and  i!=len(edge_dict[stem][1]) - 1:
                        continue  
                if L_TYPE == "three":
                    if type(edge_dict[stem][1]) in (list,tuple) and  i!=len(edge_dict[stem][1]) - 2:
                        continue 
                if ax.get_xlim()[0] < lin < ax.get_xlim()[1]:
                    ax.axvline(x=lin,ls=(5,(10,3)),color=cmap(c))
                    ax.text(lin-0.1, 0.96*ax.get_ylim()[1] +0.04*ax.get_ylim()[0],dopant+"--$"+text+"$",verticalalignment='top',horizontalalignment='right',rotation=-90,color=cmap(c))

            #if ax.get_xlim()[0] < max(edge_dict[stem][0]) < ax.get_xlim()[1]:
                #ax.axvline(x=max(edge_dict[stem][0]),ls=(5,(10,3)),color=cmap(c))
                #ax.text(max(edge_dict[stem][0])-0.1, 0.96*ax.get_ylim()[1] +0.04*ax.get_ylim()[0],dopant+"--$"+edge_dict[stem][1]+"$",verticalalignment='top',horizontalalignment='right',rotation=-90,color=cmap(c))
    if INDEP_VARIABLE is EDGE_SEPARATION:
        ax.axvline(x=0,ls=(5,(10,3)),color="black")

    if DEP_VARIABLE is R_FACTOR:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(-1,1),)
    #yticks = ax.get_yticks()
    #OOM = math.floor(math.log10(np.max(yticks)))
    #ax.set_yticks((yticks//(10**(OOM-1))*(10**(OOM-1))))
    #ax.set_yticklabels(([str(y/OOM) for y in yticks ]))

    if LEGEND:
        if BATCH is REAL_ELEMENTS:
            ax.legend(title="Dopant",fancybox=True,ncol=2,loc='upper center',bbox_to_anchor=(0.758, 1.02),handletextpad=0.01,columnspacing=0,handlelength=1.35,borderpad = 0.18,frameon=True)
        if BATCH is ELECTRON_SOURCE:
            handles, labels = ax.get_legend_handles_labels()
            handles.reverse()
            labels.reverse()
            #ax.legend(fancybox=True,ncol=1,loc='upper center',bbox_to_anchor=(0.8, 1.02),handletextpad=0.01,columnspacing=0.2,borderpad = 0.18)
            ax.legend(labels=labels,handles=handles,fancybox=True,ncol=1,loc='upper center',bbox_to_anchor=(0.8, 0.37),handletextpad=0.01,columnspacing=0.2,borderpad = 0.18)

    #ax.set_title(fig_title)
    #fig.tight_layout()
    fig.savefig(figure_output_dir + label + output_tag + figures_ext,dpi=DPI)#bbox_inches="tight", pad_inches=0)

def grow_crystals(run_params,pdb_path,allowed_atoms,identical_deviations=True,plot_crystal=True,CNO_to_N=False,S_to_N=True):
    ''' 
    Grows a crystal. The first is damaged, corresponding to I_real, the second is undamaged, corresponding to I_ideal
    Ignores atoms not in allowed_atoms
    '''
    # The undamaged crystal form factor is constant (the initial state's) but still performs the same integration step with the pulse profile weighting.
    crystal_dmged = Crystal(pdb_path,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=True, **run_params["crystal"])
    if identical_deviations:
        # we copy the other crystal so that it has the same deviations in coords
        crystal_undmged = copy.deepcopy(crystal_dmged)
        crystal_undmged.is_damaged = False
    else:
        crystal_undmged = Crystal(pdb_path,allowed_atoms,is_damaged=False,convert_excluded_elements_to_N=True, **run_params["crystal"])
    if plot_crystal:
        #crystal.plot_me(250000)
        crystal_undmged.plot_me(25000,template="plotly_dark")
    return crystal_dmged, crystal_undmged

def get_R(sim_handle,sim_handle_parent_folder,scattering_sim_parameters,allowed_atoms=["C","N","O"]):        
        run_params = copy.deepcopy(scattering_sim_parameters)
        print("Performing scattering simulation using data from damage simulation '" + sim_handle + "'.")
        exp1_qualifier = "_real"
        exp2_qualifier = "_ideal"          
        exp_name1 = sim_handle + "_" + exp1_qualifier
        exp_name2 = sim_handle + "_" + exp2_qualifier

        param_dict,_,_ = get_sim_params(sim_handle)
        start_time = param_dict["start_t"]
        end_time = param_dict["end_t"]
        energy = param_dict["energy"]
        # 
        print("Preparing XFEL...")
        experiment1 = XFEL(exp_name1,energy,**run_params["experiment"])
        experiment2 = XFEL(exp_name2,energy,**run_params["experiment"])

        print("Growing crystals at breakneck speeds...")     
        structure = PDB_STRUCTURE
        if SCATTERING_TARGET == 99:
            # Debug, use tetrapeptide (small structure, fast)
            structure = get_pdb_path(SCATTER_DIR,"tetra")
            
        crystal_real, crystal_ideal = grow_crystals(run_params,structure,allowed_atoms,plot_crystal=False)
    
        if run_params["laser"]["SPI"]:
            SPI_result1 = experiment1.spooky_laser(start_time,end_time,sim_handle,sim_handle_parent_folder,crystal_real, **run_params["laser"])
            SPI_result2 = experiment2.spooky_laser(start_time,end_time,sim_handle,sim_handle_parent_folder,crystal_ideal, **run_params["laser"])
            #TODO apply_background([SPI_result1,SPI_result2])
            damage_dict = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=True,SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
        #pulse_params = [energy,fwhm,photon_count]
        return damage_dict["R"],damage_dict["resolutions"]  # resolution index 0 corresponds to the max resolution

def get_saved_data(handle):
    # Look for save of data
    file_path = GOLDI_RESULTS_PATH+handle+".json"
    if os.path.exists(file_path):
        with open(file_path,"r") as f:
            return json.load(f)
    return None
    
    
def save_data(handle,output_data_dict,delete_old = False):
    os.makedirs(GOLDI_RESULTS_PATH,exist_ok=True)
    if output_data_dict == {}:
        return
    if delete_old:
        # Delete saves for older simulations (smaller run number) with this mol file.
        all_outputs = os.listdir(MOLECULAR_PATH)
        for match in all_outputs:
            if match.startswith(handle+"_") and match != handle:
                os.remove(GOLDI_RESULTS_PATH+match+".json")

    file_path = GOLDI_RESULTS_PATH+handle+".json"
    # Get any data from existing json file 
    data = output_data_dict
    if os.path.exists(file_path):
        with open(GOLDI_RESULTS_PATH+handle+".json","r") as f:
            saved_data_dict = json.load(f)
            data = merge(data,saved_data_dict)

    # Save data
    with open(GOLDI_RESULTS_PATH+handle+".json", 'w') as f: 
        json.dump(data, f)    


def get_data_point(ax,stem,mol_name,mode):
    im_params,_ = SCATTERING_TARGET_DICT[SCATTERING_TARGET]
    indep_variable_key = str(INDEP_VARIABLE)
    dep_variable_key = str(DEP_VARIABLE)
    if INDEP_VARIABLE == EDGE_SEPARATION:
        indep_variable_key += "L" + L_TYPE
    if DEP_VARIABLE is R_FACTOR:
        dep_variable_key = str(DEP_VARIABLE)+"-"+str(SCATTERING_TARGET)
    if DEP_VARIABLE is AVERAGE_CHARGE:
        if INTENSITY_AVERAGED: 
            dep_variable_key += '_avged'
        else:
            dep_variable_key += '_EoP'
    saved_x = saved_y = saved_end_time = None
    dat= get_saved_data(mol_name) 
    if dat is not None:
            saved_x,                              saved_y,                            saved_end_time = (
            dat.get("x").get(indep_variable_key) if dat.get("x") is not None else None,dat.get("y").get(dep_variable_key),dat.get("end_time")
        )                          
    ## Measure damage 
    if saved_y is None:              
        pl = Plotter(mol_name)
        if mode is AVERAGE_CHARGE:
            if INTENSITY_AVERAGED: 
                y = np.average(pl.get_total_charge(atoms="C")*pl.intensityData)/np.average(pl.intensityData)
            else:
                step_index = -1  # Average charge at End-of-pulse 
                y = pl.get_total_charge(atoms="C")[step_index]
        elif mode is R_FACTOR:            
            y = get_R(mol_name,MOLECULAR_PATH,im_params)[0][0]
    else: 
        y = saved_y
    # Get end times (for purposes of checking all simulations ran to completion)
    if saved_end_time is None:
        if saved_y is not None:
            pl = Plotter(mol_name)
        t = pl.timeData[-1] 
    else:
        t = saved_end_time        
    ## Get independent variable (must be last because of sourceless case)
    if INDEP_VARIABLE is PHOTON_ENERGY: 
        energy = get_sim_params(mol_name)[0]["energy"]
    if INDEP_VARIABLE is EDGE_SEPARATION:
        if type(edge_dict[stem][1]) in (list,tuple) and L_TYPE == "one":
            energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0][2]*1000
        elif type(edge_dict[stem][1]) in (list,tuple) and  L_TYPE == "two":
            energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0][1]*1000
        elif type(edge_dict[stem][1]) in (list,tuple) and  L_TYPE == "three":
            energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0][0]*1000
        else:
            energy = get_sim_params(mol_name)[0]["energy"] - max(edge_dict[stem][0])*1000
    if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
        s = get_sim_params(mol_name)[0]
        energy = s["source_energy"]
        if energy is None:
            # No source
            if (s["fluence"],s["width"]) != (fluence,width):
                print(s)
                print((s["fluence"],s["width"]),"!=",(fluence,width))
                raise ValueError("Params do not match expected values for batch")                        
        elif (s["fluence"],s["width"],s["source_fraction"],s["source_duration"]) != (fluence,width,source_fraction,source_duration):
            print(s)
            print((s["fluence"],s["source_fraction"],s["source_duration"]),"!=",(fluence,width,source_fraction,source_duration))
            raise ValueError("Params do not match expected values for batch")                    
    x = energy
    if x!= None:
        # Standard case
        x = energy/1000    
    if saved_x is not None:
        assert saved_x == x                                 
    save_dict = {}
    if saved_x is None and energy != None:
        save_dict["x"] = {indep_variable_key:x}
    if saved_y is None:    
        save_dict["y"] = {dep_variable_key:y}
    if saved_end_time is None:
        save_dict["end_time"] = t
    save_data(mol_name,save_dict,delete_old=True)
    return x,y,t

if __name__ == "__main__":
    main()

