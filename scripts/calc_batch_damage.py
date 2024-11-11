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
from core_functions import get_sim_params,get_pdb_path
from core_variables import *
import imaging_params
from mergedeep import merge


FIGWIDTH_FACTOR = 1
FIGWIDTH = 3.49751 *FIGWIDTH_FACTOR
FIGHEIGHT = FIGWIDTH*3/4
FIT = False
IGNORE_MISSING_SIMULATIONS = True  # If False, throws an error if not all specified simulations in the batch(es) are found.

matplotlib.rcParams.update({
    'figure.subplot.left':0.16/FIGWIDTH_FACTOR,
    'figure.subplot.right':1-0.02/FIGWIDTH_FACTOR,
    'figure.subplot.bottom': 0.16/FIGWIDTH_FACTOR,
    'figure.subplot.top':1-0.07/(FIGWIDTH_FACTOR*FIGHEIGHT/FIGWIDTH*4/3),
})


###
INDEP_VARIABLE = PHOTON_ENERGY # photon energy.

####### User parameters #########
INDEP_VARIABLE = PHOTON_ENERGY
DEP_VARIABLE = AVERAGE_CHARGE # | 0: mean carbon charge (at end of pulse) | 1: R factors (performs scattering simulations for each) |
INTENSITY_AVERAGED = True # Used if DEP_VARIABLE = 0 (average charge).
SCATTERING_TARGET = 0 # Used if DEP_VARIABLE = 1 (R factor).  | 0: unit lysozyme (light atoms) | 1: 3x3x3 lysozyme (light atoms no solvent)|

## Graphical
LEGEND = True
LABEL_TIMES = False # Set True to graphically check simulation end times are all aligned.

## Numerical
NORMALISING_STEM = "lys_all_light"  #Stem (str) to normalise traces by. (s.t. normalised trace becomes horizontal line). If None, does not normalise.

ylim=[None,None]

# Batch stems and index range (inclusive). currently assuming form of key"-"+n+"_1", where n is a number in range of stem[key]
# e.g. item "Carbon":[2,4] means to plot simulations corresponding to outputs Carbon-2, Carbon-3, Carbon-4, mol files. Last simulation output is used in case of duplicates (highest run number).
stem_dict = {
        "lys_full":[1,150],
        "lys_all_light":[1,150],
}
# Directory relative to AC4DC/output/__Molecular
folder_dict = {
    "lys_full":"lys_full/",
    "lys_all_light":"lys_all_light/",
}
# Each item corresponds to a list of additional simulations to add on e.g. "Carbon":[8,12,5] adds those 3 mol files.
additional_points_dict = {
    #"example":[8,12,5]
}    

#Integers corresponding to colours in discrete cmap.
col_dict = {
    "lys_all_light":0,
    "lys_full":1,
}

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
    assert all([k in folder_dict.keys() and k in col_dict for k in stem_dict.keys()])
    for val in folder_dict.values():
        assert(val[-1] == "/")

    for key in additional_points_dict.keys():
        if key not in stem_dict.keys():
            stem_dict[key] = None
    # Plot tag
    plot_type = {AVERAGE_CHARGE:'avg_charge_plot',R_FACTOR:'R_dmg'}      

    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    
    batches = {}
    sim_tags = {}
    # Get folders names in molecular_path
    valid_folder_names= True
    for key,val in stem_dict.items():
        all_outputs = os.listdir(MOLECULAR_PATH+folder_dict[key])
        data_folders = []
        # Get folder name, excluding run tag
        sim_tags[key] = []
        if val is not None:
            sim_tags[key] = [*range(val[0],val[1]+1)] 
        if key in additional_points_dict:
            for sim_tag in additional_points_dict[key]:
                sim_tags[key].append(sim_tag)
        for n in sim_tags[key]:
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
                    sim_tags[key][i] = None
                    continue
                break
            data_folders[i] = handle + "_"+str(max(run_nums))
        batches[key] = [f for f in data_folders if f is not None]

    if not IGNORE_MISSING_SIMULATIONS:
        assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder. See input error message(s)"
    fig_title = "".join([stem.split("-")[0].split("_")[-1] for stem in batches.keys()])
    label = fig_title+"_"+plot_type[DEP_VARIABLE]
    plot(batches,sim_tags,label,dname_Figures,mode=DEP_VARIABLE)

def plot(batches,sim_tags,label,figure_output_dir,mode = 0):
    ylabel = {AVERAGE_CHARGE:"Average carbon's charge",R_FACTOR:"$R_{dmg}$"}
    '''
    Arguments:
    '''    
    ############
    # File/directory names
    #######  
    figures_ext = "" #.png
    fname_charge_conservation = "charge_conservation"
    fname_free = "free"
    fname_HR_style = "HR_style"
    fname_bound_dynamics = "bound_dynamics"


    fig, ax = plt.subplots(figsize=(FIGWIDTH,FIGHEIGHT))
    cmap = plt.get_cmap('tab10')

    output_tag = ""
    if INTENSITY_AVERAGED:
        output_tag+="-IAvg"    
    if DEP_VARIABLE is R_FACTOR:
        _,output_tag = SCATTERING_TARGET_DICT[SCATTERING_TARGET]
    if NORMALISING_STEM not in [None,False]:
        output_tag += "-normed"
    
    if NORMALISING_STEM not in [None,False]:
        norm = {}
        mol_names = batches[NORMALISING_STEM]
        for mol_name,sim_tag in zip(mol_names,sim_tags[NORMALISING_STEM]):
            norm[sim_tag]=get_data_point(ax,NORMALISING_STEM,mol_name,mode)


    for stem, mol_names in batches.items():
        print("Plotting "+stem+" trace.")
        c = col_dict[stem]
        X = [int(t) for t in sim_tags[stem]]
        Y = []
        times = []  # last times of sims, for checking everything lines up.
        if stem == NORMALISING_STEM:
            for energy,y,t in norm.values():
                Y.append(y)
                times.append(t)
        else:
            for mol_name in mol_names:
                energy,y,t = get_data_point(ax,stem,mol_name,mode)
                Y.append(y)  
                times.append(t)   
                    
        if mode is AVERAGE_CHARGE:
            print("Charges:",Y)
        if mode is R_FACTOR:
            print("R factors:",Y)
        if NORMALISING_STEM not in [None,False]:
            for i, x in enumerate(X):
                    Y[i]/=norm[x][1]        
        _label = stem
        ax.scatter(X,Y,label=_label,color = cmap(c))

        if LABEL_TIMES:
            for i, txt in enumerate(times):
                ax.annotate(txt, (X[i],Y[i]))
        _ylab = ylabel[mode]
        if NORMALISING_STEM:
            _ylab += " (relative)"
        ax.set_ylabel(_ylab)
        ax.set_xlabel("sim num")

    #ax.set_title(fig_title)
    #fig.tight_layout()

    ax.ticklabel_format(style='sci', axis='y', scilimits=(-1,1),)

    if LEGEND:
            #ax.legend(title="Dopant",fancybox=True,ncol=2,loc='upper center',bbox_to_anchor=(0.75, 1.02),handletextpad=0.01,columnspacing=0.2,borderpad = 0.18)
            ax.legend(fancybox=True,ncol=1,loc='upper center',bbox_to_anchor=(0.8, 1.02),handletextpad=0.01,columnspacing=0.2,borderpad = 0.18)

    fig.savefig(figure_output_dir + label + output_tag + figures_ext)#bbox_inches="tight", pad_inches=0)


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
            R, binned_R_res, cc, resolutions,binned_R_q = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=True,SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
        #pulse_params = [energy,fwhm,photon_count]
        return R,resolutions  # resolution index 0 corresponds to the max resolution

def get_saved_data(handle,subdir):
    # Look for save of data
    results_folder = GOLDI_RESULTS_PATH+subdir
    file_path = results_folder+handle+".json"
    if os.path.exists(file_path):
        with open(file_path,"r") as f:
            return json.load(f)
    return None
    
    
def save_data(handle,output_data_dict,subdir,delete_old = False):
    results_folder = GOLDI_RESULTS_PATH+subdir
    os.makedirs(results_folder,exist_ok=True)
    if output_data_dict == {}:
        return    
    # Delete saves for older simulations (smaller run number) with this mol file.
    all_outputs = os.listdir(MOLECULAR_PATH)
    for match in all_outputs:
        if match.startswith(handle+"_") and match != handle:
            os.remove(results_folder+match+".json")

    file_path = results_folder+handle+".json"
    # Get any data from existing json file 
    data = output_data_dict
    if os.path.exists(file_path):
        with open(results_folder+handle+".json","r") as f:
            saved_data_dict = json.load(f)
            data = merge(data,saved_data_dict)

    # Save data
    with open(results_folder+handle+".json", 'w') as f: 
        json.dump(data, f)    



def get_data_point(ax,stem,mol_name,mode):
    batch_path = MOLECULAR_PATH+folder_dict[stem]
    im_params,_ = SCATTERING_TARGET_DICT[SCATTERING_TARGET]
    indep_variable_key = str(INDEP_VARIABLE)
    dep_variable_key = str(DEP_VARIABLE)
    if INTENSITY_AVERAGED:
        dep_variable_key += '_avged'
    else:
        dep_variable_key += '_EoP'  
    if DEP_VARIABLE is R_FACTOR:
        dep_variable_key = str(DEP_VARIABLE)+"-"+str(SCATTERING_TARGET)
    saved_x = saved_y = saved_end_time = None
    dat= get_saved_data(mol_name,subdir=folder_dict[stem]) 
    if dat is not None:
            saved_x,                              saved_y,                            saved_end_time = (
            dat.get("x").get(indep_variable_key),dat.get("y").get(dep_variable_key),dat.get("end_time")
        )                          
    print(dat)
    ## Measure damage 
    if saved_y is None:              
        pl = Plotter(mol_name,abs_molecular_path=batch_path)
        if mode is AVERAGE_CHARGE:
            if INTENSITY_AVERAGED: 
                y = np.average(pl.get_total_charge(atoms="C")*pl.intensityData)/np.average(pl.intensityData)
            else:
                step_index = -1  # Average charge at End-of-pulse 
                y = pl.get_total_charge(atoms="C")[step_index]
        elif mode is R_FACTOR:            
            y = get_R(mol_name,batch_path,im_params)[0][0]
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
        energy = get_sim_params(mol_name,molecular_path=batch_path)[0]["energy"]            
    if energy!= None:
        # Standard case
        x = energy/1000    
    if saved_x is not None:
        assert saved_x == x                                 
    save_dict = {}
    if saved_x is None:
        save_dict["x"] = {indep_variable_key:x}
    if saved_y is None:    
        save_dict["y"] = {dep_variable_key:y}
    if saved_end_time is None:
        save_dict["end_time"] = t
    save_data(mol_name,save_dict,delete_old=True,subdir=folder_dict[stem])
    return x,y,t
        
if __name__ == "__main__":
    main()

