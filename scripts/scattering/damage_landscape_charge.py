#%%
#Generates damage landscapes by performing scattering simulation across multiple simulations, either specifically chosen or through an entire folder.

import numpy as np
import os.path as path
import os
from shutil import rmtree
import inspect
import copy
import imaging_params
import time
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import colorcet as cc; import cmasher as cmr
import matplotlib.pyplot as plt
import pickle
import inspect
src_file_path = inspect.getfile(lambda: None)
my_dir = path.abspath(path.join(src_file_path ,"../")) + "/"
from core_functions import get_sim_params, get_pdb_paths_dict
from scatter import XFEL,Crystal,stylin
from plotter_core import Plotter
#############
COLUMNWIDTH = 523.5307
FIGWIDTH = COLUMNWIDTH/2
FIGHEIGHT = FIGWIDTH*10/16
FIGWIDTH -=20 
margin=go.layout.Margin(
    l=0, 
    r=0, 
    b=0, 
    t=0, 
)

#Stupid thing have to do with plotly to get the plot content to have the same scale.
BUFFER1 = 0 
BUFFER2 = -15
BUFFER3 = 0
show_colorbar = False

SCALE = 8 

###############
# full sim that does point sampling, thus assumes a "smooth" distribution. Thus good for low cell count input. Num cells is entire whole target.
# Compared with non-SPI (crystal) which assumes crystal and does bragg points. Num cells is supercell.
#TODO rename the diff sims, 'non-SPI' unclear.
#TODO currently all cells are used for calculation of R, but we may wish to ignore cells outside the considered radius? Though I haven't see anyone do this, it would depend on what image reconstruction algos do I guess. But we get empty corners with high res. limits so


MODE_DICT = {0:'crystal',1:'spi',2:'both'}

PDB_PATHS = get_pdb_paths_dict(my_dir) # <value:> the target that should be used for <key:> the name of simulation output batch folder.

DATA_FOLDER = "dmg_data/charge_data/"
PLOT_FOLDER = "../../output/_Graphs/R_plots/"

def multi_damage(params,pdb_path,allowed_atoms_1,CNO_to_N,S_to_N,same_deviations,plasma_batch_handle = "", plasma_handles = None, sctr_results_batch_dir = None, get_R_only=True, realistic_crystal_growing_mode = False,specific_energy = None,eop_charge_only=False,element_considered="C"):
    '''
    NB: "Plasma" is used as a shortened way to refer to the output of the damage simulation  
    parameters:

    If plasma_batch_handle is not provided, handle_list must be provided to prevent.
    plasma_batch_handle<optional>:  name of folder containing individual simulation output folders. If not provided, with use handle_list in the default directory.
    handle_list <optional>: list of folder names to iterate through within the batch folder. If not provided, all are used. 
    '''
    if plasma_batch_handle == "":
        assert plasma_handles is not None

    if params["laser"]["SPI"] == False:  
        # Folder structure:
        #.../'root_results_dir'/'sctr_results_batch_dir'/<experiment_results>/<orientation>
        # Folder containing the folders corresponding to each batch of handles.
        root_results_dir = path.join(my_dir,"results_R_plotting/")+ "/"
        if path.isdir(root_results_dir):
            rmtree(root_results_dir)
        if path.isfile(root_results_dir):
            raise Exception("target directory is existing file's name.")
       
    
        # Determine the scatter results' folders' parent directory. (used for infinite crystal only)
        if sctr_results_batch_dir is None:
            version_number = 1
            count = 0
            while True:
                if count > 99:
                    raise Exception("could not find valid parent folder name in " + str(count) + " loops")
                sctr_results_batch_dir = root_results_dir+plasma_batch_handle +"_v"+str(version_number) +"/" # needs to be synced with other functions
                if path.exists(path.dirname(sctr_results_batch_dir)):
                    version_number+=1
                    count+=1
                    continue 
                break
        os.makedirs(sctr_results_batch_dir)
    # Run imagings for each simulation 
    exp1_qualifier = "_real"
    exp2_qualifier = "_ideal"        

    sim_data_batch_dir = path.abspath(path.join(my_dir ,"../../output/__Molecular/"+plasma_batch_handle)) + "/"
    sim_input_dir = path.abspath(path.join(my_dir ,"../../input/")) + "/"
    if plasma_handles == None:
        plasma_handles = []
        print("Iterating through entire batch of plasma simulation results.")
        for i, sim_handle in enumerate(os.listdir(sim_data_batch_dir)):
            plasma_handles.append(sim_handle)
    else:
        print("Iterating through",plasma_handles)

    start_time, end_time, energy, fwhm, photon_count = ([None]*len(plasma_handles) for _ in range(5))
    necessary_files = ["freeDist.csv","intensity.csv"]
    for i, sim_handle in enumerate(plasma_handles):
        assert (sim_handle != "" and sim_handle != None)
        # Check some of the expected files are present (not checking for existence of bound state files at present.) 
        plasma_output = sim_data_batch_dir + sim_handle
        assert path.isdir(plasma_output), "Directory not found for "+plasma_output
        for dat_file in necessary_files: 
            assert path.isfile(plasma_output + "/" +dat_file), "Missing "+dat_file+" in "+plasma_output
        # Index parameters of simulation
        param_dict, param_name_list,unit_list = get_sim_params(sim_handle,sim_input_dir,sim_data_batch_dir)
        start_time[i],end_time[i],energy[i],fwhm[i],photon_count[i] = param_dict["start_t"],param_dict["end_t"],param_dict["energy"],param_dict["width"],param_dict["fluence"]
    dmg_data = []
    pulse_params=[]
    names = []
    resolutions_vect = [] # Sometimes resolution we provided in input params is cutoff due to a low energy being used, so we need to save them all to be safe.
    for i, sim_handle in enumerate(plasma_handles):        
        # TODO
        # if params["laser"]["SPI"]:
        #     def apply_background():

        if specific_energy is not None and energy[i] != specific_energy:
            continue
        names.append(sim_handle)
        run_params = copy.deepcopy(params)
        
        #TODO print properly
        print(sim_handle)
        print(energy[i],fwhm[i],photon_count[i])
        print(param_name_list)
        print(unit_list)

        pulse_params.append([energy[i],fwhm[i],photon_count[i]])
        pl = Plotter(sim_handle,abs_molecular_path=sim_data_batch_dir,load_specific_atoms=[element_considered],sample_end_points_only=eop_charge_only)
        c = pl.get_total_charge(atoms=element_considered)[-1]  # Average charge at end of pulse 
        if eop_charge_only:
            c_IA = None
        else:
            c_IA = np.average(pl.get_total_charge(atoms=element_considered)*pl.intensityData)/np.average(pl.intensityData)  # Intensity scaled Average charge 
            
        dmg_data.append([c,c_IA])   
        print("Average final carbon charge:",c)
        print("Intensity-averaged carbon charge:",c_IA)
        print("")
    
    return np.array(pulse_params,dtype=float),np.array(dmg_data,dtype=object), np.array(resolutions_vect,dtype=object), names, param_name_list

def save_data(fname, data):
    os.makedirs(DATA_FOLDER,exist_ok=True)

    with open(DATA_FOLDER+fname+".pickle", 'wb') as file: 
        pickle.dump(data, file)

    
def load_df(fname,check_batch_nums=True):
    with open(DATA_FOLDER+fname +".pickle", "rb") as file:
        pulse_params, dmg_data,resolutions,names,param_name_list = pickle.load(file) 
    
    dmg_data = np.array(dmg_data,dtype=float)
    
    
    # resolutions_slice = resolutions[:,resolution_idx]
    
    # assert np.all(resolutions_slice==resolutions_slice[0]),   
    # resolution = resolutions_slice[0]
    # print("Resolution: ", resolution, "Angstrom")
    log_photons = np.log(pulse_params[:,2])-np.min(np.log(pulse_params[:,2]))
    log_photons += np.max(log_photons)/2
    # photon_size_thing = dmg_data[:,2] + np.max(dmg_data[:,2])/10
    photon_size_thing = 1+np.square(np.log(pulse_params[:,2])-np.min(np.log((pulse_params[:,2]))))

    assert param_name_list[3] == "R"
    photon_measure = param_name_list[2]
    photon_data = [1e12*p for p in pulse_params[:,2]]
    df = pd.DataFrame({
        "name": names, 
        "energy": pulse_params[:,0]/1000,
        "fwhm": pulse_params[:,1],
        photon_measure: photon_data,
        "eop_charge": dmg_data[:,0], # End of pulse average carbon charge
        "IA_charge": dmg_data[:,1], # Intensity-averaged carbon charge
        #"log_photons": log_photons,
        "photon_size_thing": photon_size_thing, 
        "_": pulse_params[:,0]*0+1, # dummy column for size.
    })
    # Check if missing any files.
    if check_batch_nums:
        nums = []
        for elem in df["name"]:
            elem = elem.split('-')[-1]
            elem = elem.split('_')[0]
            # extract only numbers.
            elem = ''.join([l for l in elem if l.isnumeric()])
            nums.append(int(elem))
        nums.sort()
        print("R values found for",len(nums),"simulations.")
        for i, elem in enumerate(nums):
            if i == 0: 
                if elem != 1: 
                    print("mising file number 1?")
            elif nums[i-1] != elem -1:
                print("Missing file? File number",elem,"was found after",nums[i-1])  
    df = df.sort_values(by=["energy","fwhm",photon_measure])
    df = df.reset_index(drop=True)
    df.resolution = resolution
    return df   

from math import log10,floor
def plot_that_funky_thing(df,cmin=0.1,cmax=0.3,clr_scale="temps",use_neutze_units=False,name_of_set="",energy_key=9000,photon_key=1e14,fwhm_key=50,cmin_contour=0,cmax_contour=None,contour_colour = "amp",contour_interval=0.05,dmg_measure="eop_charge",photon_min=None,max_fwhm=None,connect_contour_gaps=False,round_fluence=False,**kwargs):
    out_folder = PLOT_FOLDER
    os.makedirs(out_folder,exist_ok=True)
    ext = ".png"
    photon_measure = df.columns[3]
    # TODO temporary patch - rounds to two sig figs
    if round_fluence:
        df[photon_measure] = [round(p, 1 -int(floor(log10(abs(p))))) for p in df[photon_measure]] 
    if photon_measure[:5] == "Count":
        photon_measure_default = "$\\text{Fluence (ph um}^{-2}\\text{)}$"
        if use_neutze_units:
            df[photon_measure] = [7.854e-3*p for p in df[photon_measure]] # square micrometre to 100 nm diameter spot
            photon_measure = "photons per 100 nm spot"
        else:
            photon_measure = photon_measure_default
        df.rename(columns=dict(Count=photon_measure),inplace=True)
        if max_fwhm is None:
            max_fwhm = np.inf
        df= df[df['fwhm'] <= max_fwhm]
        print(df)
    else:
        if use_neutze_units:
            print("Neutze units only implemented for photon count measure. Disabled.")
            use_neutze_units = False    
    #3D - no idea why x axis is bugged.
    fig = px.scatter_3d(df, x=photon_measure, y='fwhm', z='energy',
              color=dmg_measure,  color_continuous_scale=clr_scale, range_color=[cmin,cmax],**kwargs)
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.25, y=-1.25, z=1.25)
    )
    fig.update_layout(scene_camera=camera)      
    fig.update_layout(title=name_of_set+"-3D")      
    fig.show()    
    #print(inspect.getmembers(fig["layout"]["title"]),dir(fig["layout"]["title"]),type(fig["layout"]["title"]))
    fig.write_html(out_folder+fig["layout"]["title"]["text"]+".html")

    # Plot designed for small number of simulations to compare.
    if len(df['name']) <= 10:
        label_col = 'name'   
        fig = px.scatter(df, text=label_col,x='fwhm', y=dmg_measure, symbol='energy', size='photon_size_thing',#'log_photons',
                color=dmg_measure, color_continuous_scale=clr_scale, range_color=[cmin,cmax],size_max=15,opacity=1,**kwargs) 
        fig.update_layout(
            showlegend=True,
            legend=dict(
            yanchor="top",
            y=0.98,
            xanchor="left",
            x=1.13,
        ))      
        fig.update_traces(textposition='top center')
        fig.update_layout(title=name_of_set+df['name'][0]+''.join(df['name'][:min(2,len(df['name']))]) + "-comparison")  
        fig.show()
        fig.write_image(out_folder+fig["layout"]["title"]["text"]+ext,scale=SCALE)

    # 'intuitive' plot, where:
    # Dot size represents num photons
    # colour is R (damage)
    # horiz axis is fwhm
    # vert axis is energy of photons 
    df = df.sort_values('photon_size_thing',ascending=False)
    fig = px.scatter(df, x='fwhm', y='energy', size='photon_size_thing',
              color=dmg_measure, color_continuous_scale=clr_scale, range_color=[cmin,cmax],size_max=30,opacity=1,
              **kwargs)
    fig.update_layout(title=name_of_set+"-Flat")     
    fig.show()    
    fig.write_image(out_folder+fig["layout"]["title"]["text"]+ext,scale=SCALE)


    ranges = {}
    for col in df:
        if col == "name":
            continue
        #rng = df[col].agg(['min','max'])
        rng = [df[col].min(),df[col].max()]
        rng[0]-=rng[1]/10 
        if photon_min != None and col == photon_measure:
            rng[0] = photon_min
        rng[0] = max(0,rng[0])
        rng[1]+= rng[1]/10
        ranges[col] = rng

    if cmax_contour == None:
        cmax_contour = ranges[dmg_measure][1]

    contour_args = dict(
        colorscale = contour_colour,#clr_scale, #'electric',
        line_smoothing=0,
        connectgaps = connect_contour_gaps,
        zmin = 0, zmax = 0.4,
        #contours_coloring='heatmap',
        contours = go.contour.Contours(start = contour_interval, end= cmax_contour, size = contour_interval),
        colorbar = dict(title = dmg_measure, tickvals = np.arange(contour_interval,cmax_contour+0.00001,contour_interval)),
    )
    original_df = df
    #get data frames for each energy
    unique_E = np.sort(df.energy.unique())
    unique_fwhm = np.sort(df.fwhm.unique())
    unique_photons = np.sort(df[photon_measure].unique())
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_E}    
    for key in data_frame_dict.keys():
        data_frame_dict[key] = df[:][df.energy == key]    

    #plot separate energies
    # Get min and max to use for axes here, so we can be consistent.
    """
    for energy in unique_E:
        df = data_frame_dict[energy] 
        fig = px.scatter(df,x='fwhm',y=photon_measure,
            color='R',  color_continuous_scale=clr_scale, range_color=[cmin,cmax], size='_',size_max=15,opacity=1,**kwargs)
        fig.update_xaxes(range=ranges['fwhm'])
        fig.update_yaxes(type="log",range=np.log10(ranges[photon_measure]))
        fig.update_layout(title = name_of_set + ", E = "+str(energy)+" eV")
        fig.show()
        fig.write_image(out_folder+fig["layout"]["title"]["text"]+"-Sctr"+ext,scale=SCALE)
    """
    # Contour: Fluence-FWHM
    if energy_key in data_frame_dict.keys(): 
        df = data_frame_dict[energy_key]
        fig = go.Figure(data =  
            go.Contour(
                #z=R_mesh,
                #x=unique_fwhm, 
                #y=unique_photons,
                z = df[dmg_measure],
                x = df["fwhm"],
                y = [p for p in df[photon_measure]],

                **contour_args,
            ))
        fig.update_xaxes(title="FWHM (fs)",type="log")
        fig.update_yaxes(title=photon_measure)
        fig.update_yaxes(type="log",range=np.log10(ranges[photon_measure]),tickvals = [1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15],tickformat = '.0e')
        fig.update_layout(width = FIGWIDTH + BUFFER1, height = FIGHEIGHT,font_family="Serif", font=dict(size=10),margin=margin,)
        fig.update_traces(showscale=show_colorbar)
        #fig.update_layout(title= name_of_set + ", E = "+str(energy_key)+" eV")
        title = name_of_set + ", E = "+str(energy_key)+" eV"
        print(title)
        fig.show()
        fig.write_image(out_folder+title+ext,scale=SCALE)
        df = original_df
        
    #Contour Energy-FWHM
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_photons}   
    default_unit_photon_key = '{:.2e}'.format(photon_key)
    if use_neutze_units:
        photon_key *=7.854e-3
    photon_key_for_title = '{:.2e}'.format(photon_key)
    if photon_key in data_frame_dict.keys():  
        for key in data_frame_dict.keys():
            data_frame_dict[key] = df[:][df[photon_measure] == key]    
        df = data_frame_dict[photon_key]
        print(df)

        fig = go.Figure(data =  
            go.Contour(
                #z=R_mesh,
                #x=unique_fwhm, 
                #y=unique_photons,
                z = df[dmg_measure],
                y = df["fwhm"],
                x = df["energy"],
                **contour_args,
            ))
        fig.update_layout(width = FIGWIDTH+BUFFER2, height = FIGHEIGHT,font_family="Serif", font=dict(size=10),margin=margin,)
        fig.update_traces(showscale=show_colorbar)
        #fig.update_layout(title = name_of_set + ", " +str(photon_key_for_title)+" " +photon_measure)
        title = name_of_set + ", FWHM = "+str(fwhm_key)+" fs"
        print(title)
        fig.update_yaxes(title="FWHM (fs)",type="log")
        fig.update_xaxes(title="Energy (keV)")

        fig.show()   
        #fig.write_image(out_folder+name_of_set + ", " +str(default_unit_photon_key)+" " +photon_measure_default+ext,scale=SCALE)
        fig.write_image(out_folder+name_of_set + ", " +str(default_unit_photon_key)+" " +"ph"+ext,scale=SCALE)
        df = original_df

    #Contour: Energy-Fluence
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_fwhm}    
    if fwhm_key in data_frame_dict.keys(): 
        for key in data_frame_dict.keys():
            data_frame_dict[key] = df[:][df.fwhm == key]    
        df = data_frame_dict[fwhm_key]

        fig = go.Figure(data =  
            go.Contour(
                #z=R_mesh,
                #x=unique_fwhm, 
                #y=unique_photons,
                z = df[dmg_measure],
                x = df["energy"],
                y = [p for p in df[photon_measure]],
               **contour_args,
            ))
        fig.update_layout(width = FIGWIDTH+ BUFFER3, height = FIGHEIGHT,font_family="Serif", font=dict(size=10),margin=margin,)
        fig.update_traces(showscale=show_colorbar)
        #fig.update_layout(title=name_of_set + ", FWHM = "+str(fwhm_key)+" fs")
        title = name_of_set + ", FWHM = "+str(fwhm_key)+" fs"
        print(title)
        fig.update_xaxes(title="Energy (keV)")
        fig.update_yaxes(title=photon_measure,type="log",range=np.log10(ranges[photon_measure]),tickvals = [1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15],tickformat = '.0e')
        fig.show()    
        fig.write_image(out_folder+title+ext,scale=SCALE)    
        df = original_df
#%%
#------------Generate data--------------------

if __name__ == "__main__":
    # Initialise kwargs
    kwargs = dict(
        plasma_batch_handle = "", 
        plasma_handles = None, 
        sctr_results_batch_dir = None, 
        get_R_only=True,
        eop_charge_only=False, # Only calculate the end of pulse charge (much faster as doesn't load whole data file)
        element_considered = "Gd_fast"#"Gd_fast"#"C"
    )


    batch_mode = False # Just doing this as I want to quickly switch between doing batches and comparing specific runs.
    multi_batch_mode = False

    mode = 1  #0 -> infinite crystal, 1 -> finite crystal/SPI, 2-> both  
    same_deviations = True # whether same position deviations between damaged and undamaged crystal (SPI only)
    imaging_params_preset = imaging_params.default_dict

    # Change options here... #TODO make clearer this is for user options
    assert (batch_mode == False or multi_batch_mode == False)
    if batch_mode or multi_batch_mode:
        allowed_atoms = ["C","N","O"] 
        CNO_to_N = False
        S_to_N = False
        batch_dir = None # Optional: Specify existing parent folder for batch of results, to add these orientation results to.
        pdb_path = PDB_PATHS["lys"]
        # Params set to batch_handle
    if batch_mode:
        batch_handle = "SH2_Fe" 
        kwargs["plasma_batch_handle"] = batch_handle
        fname = batch_handle
    # Equivalent to performing `batch_mode` for each element of `batch_handles` as `batch_handle`
    elif multi_batch_mode:
        #batch_handles = ["SH2_Gd","SH2_Fe","SH2_Se","SH2_S","SH2_N","SH2_Zn"] 
        batch_handles = ["SH2_S","SH2_N"] 
    else: # Compare specific simulations
        fname = "comparison"
        allowed_atoms = ["C","N","O"]; S_to_N = True
        #allowed_atoms = ["N"]; S_to_N = True
        #allowed_atoms = ["C","N","O","S"]; S_to_N = False
        CNO_to_N = False
        kwargs["plasma_batch_handle"] = ""
        #kwargs["plasma_handles"] = ["lys_nass_3","lys_nass_no_S_1"]    #Note that S_to_N must be true to compare effect on nitrogen R factor. Comparing S_to_N true with nitrogen only sim, then S_to_N false with nitrogen+sulfur sim, let's us compare the true effect of sulfur on R facator.
        #kwargs["plasma_handles"] = ["lys_nass_no_S_2","lys_nass_6","lys_nass_Gd_16"]  
        #kwargs["plasma_handles"] = ["lys_nass_no_S_2","lys_nass_HF","lys_nass_Gd_HF"]  
        #kwargs["plasma_handles"] = ["lys_nass_gauss","lys_nass_square"]  
        #kwargs["plasma_handles"] = ["lys_galli_LF_no_Gd_6","lys_galli_LF_17"]
        #kwargs["plasma_handles"] = ["lys_galli_HF_15","lys_galli_LF_17","lys_galli_HF_no_Gd_2","lys_galli_LF_no_Gd_6",]
        kwargs["plasma_handles"] = ["lys_galli_HF_15","lys_galli_LF_17","lys_galli_LF_19"]
        #kwargs["plasma_handles"] = ["lys_nass_HF","lys_nass_Gd_HF"]  
        #kwargs["plasma_handles"] = ["lys_full-typical","lys_all_light-typical"]  
        #kwargs["plasma_handles"] = ["glycine_abdullah_4"]
        #pdb_path = PDB_PATHS["fcc"]
        pdb_path = PDB_PATHS["lys"]
        #pdb_path = PDB_PATHS["lys_solvated"]
        #pdb_path = PDB_PATHS["glycine"]
        #pdb_path = PDB_PATHS["tetra"]
    

    num_loops = 1 if not multi_batch_mode else len(batch_handles)
    for i in range(num_loops):
        if multi_batch_mode: 
            kwargs["plasma_batch_handle"] = fname = batch_handles[i] 

        if MODE_DICT[mode] == "crystal" or MODE_DICT[mode] == "both":
            imaging_params_preset["laser"]["SPI"] = False         
            scatter_data = multi_damage(imaging_params_preset,pdb_path,allowed_atoms,CNO_to_N,S_to_N,same_deviations,**kwargs)
        if MODE_DICT[mode] == "spi" or MODE_DICT[mode] == "both":
            imaging_params_preset["laser"]["SPI"] = True         
            scatter_data = multi_damage(imaging_params_preset,pdb_path,allowed_atoms,CNO_to_N,S_to_N,same_deviations,**kwargs)
    
        save_data(fname,scatter_data)
    

#--------------------------------------
#%%
#------------Plot----------------------
if __name__ == "__main__":
    # data_names corresponds to fname passed to save_data in previous cell. "comparison" if not batch mode, otherwise
    #data_names = ["SH2_Gd","SH2_Fe","SH2_Se","SH2_S","SH2_N","SH2_Zn"] 
    #data_names = ["SH2_N"] 
    data_names = ["comparison"]   
    #data_names = batch_handles
    #data_names = ["SH2_Se"]
    batch_mode = False; mode = 1  #TODO store batch_mode and mode in saved object.
    damage_measure = "eop_charge"
    cmax_contour = 6; contour_interval = 0.5
    #damage_measure = "IA_charge"
    #cmax_contour = 3; contour_interval = 0.25
    connect_contour_gaps=False; round_fluence = True # Warning: enabling these options will create misleading graphs!
    for data_name in data_names:
        name_of_set = data_name
        resolution = 1.94 #1.9 2.8
        df = load_df(data_name,check_batch_nums=batch_mode) # resolution list is orders from best to worst resolutions.
        plot_2D_constants = dict(energy_key = 15000, photon_key = 1e12,fwhm_key = 15)
        neutze = False
        if MODE_DICT[mode] != "spi":
            print("-----------------Crystal----------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
            plot_that_funky_thing(df,0,0.20,"temps",name_of_set=name_of_set,**plot_2D_constants,template="plotly_dark",connect_contour_gaps=connect_contour_gaps,round_fluence=round_fluence,use_neutze_units = neutze,cmax_contour=cmax_contour,contour_interval=contour_interval,dmg_measure = damage_measure) # 'electric' #"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
        if MODE_DICT[mode] != "crystal":
            print("-------------------SPI------------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
            plot_that_funky_thing(df,0,0.20,"temps",name_of_set=name_of_set,**plot_2D_constants,template="plotly_dark",connect_contour_gaps=connect_contour_gaps,round_fluence=round_fluence,use_neutze_units=neutze,cmax_contour=cmax_contour,contour_interval=contour_interval,dmg_measure = damage_measure) # 'electric'#"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
    print("Done")


#%%
#------------Generate damage landscape difference maps----------------------
contour_interval_delta_difference = 0.05 # 0.01
contour_interval_percentage_difference = 0.05
if __name__ == "__main__":
    show_colorbar = False
    num_remaining_electrons_mode = True
    connect_contour_gaps=False; round_fluence = True  # Warning: enabling these options will create misleading graphs!
    batch_mode = True; mode = 1
    resolution = 1.9#2.2
    photon_min=None # seems slightly bugged
    fname1 = "SH2_S"; fname2 = "SH2_N"; 
    outfig_base_name = "X_impact" #Output base name
    #fname1 = "lys_full"; fname2 = "lys_all_light"
    plot_2D_constants = dict(energy_key = 15000, photon_key = 1e12,fwhm_key = 15)
    compared_damage_measure = "eop_charge" #"eop_charge" or "IA_charge" (end of pulse average carbon charge or intensity-averaged carbon charge)
    charge_cutoff = 1 # Just for printing alternate statistics, not important
    max_fwhm = 100

    #for percentage_difference in [False,True]:
    for percentage_difference in [True]:
        if percentage_difference:
            contour_interval = contour_interval_percentage_difference
            #contour_colour = "viridis"
            contour_colour = "plasma"
            vmax =  1
        else:
            contour_interval = contour_interval_delta_difference
            contour_colour = "plasma"
            vmax = 1.5

        #name_of_set = "Difference" + "_"+fname1 + "_" + fname2
        name_of_set = outfig_base_name
        df1 = load_df(fname1, check_batch_nums=batch_mode)
        df2 = load_df(fname2,  check_batch_nums=batch_mode)
        if num_remaining_electrons_mode:
            df_old = df1
            df1 = df2 
            df2 = df_old
            df1[compared_damage_measure] = 6 - df1[compared_damage_measure]
            df2[compared_damage_measure] = 6 - df2[compared_damage_measure]
        neutze = False
        df_diff = copy.deepcopy(df1)
        with pd.option_context('display.max_rows', None, 'display.max_columns', None):  #
            print(df1)
            print(df2)

        for col in df1.columns:
            if col != "eop_charge" and col != "IA_charge" and col!="name":
                print("Checking all datapoints align...")
                assert np.array_equal(df_diff[col], df2[col]) and np.array_equal(df1[col], df_diff[col])
                
            elif col != "name":
                if percentage_difference:
                    #df_diff[col] = df1[col]/df2[col] - 1
                    df_diff[col] = (df1[col]-df2[col])/df1[col]
                else:
                    df_diff[col] -= df2[col]
                    if col == compared_damage_measure:
                        num_non_pos_count = np.sum(df_diff[col] <= 0)
                        print("NUMBER OF damage_diff <= 0:", num_non_pos_count)  
                        print("AVERAGE "+compared_damage_measure,np.sum(df1[col]/df2[col])/len(df1[col]))
                        print("MIN/MAXIMUM",np.min(df1[col]/df2[col]),"|",np.max(df1[col]/df2[col]))
                        print("sample STD",np.std(df1[col]/df2[col],ddof=1))
                        print("====Cutoff v======")
                        mask = np.array(df1[col] > charge_cutoff,dtype=bool)
                        print(np.sum(mask))
                        print("AVERAGE "+compared_damage_measure,np.nansum(np.array(df1[col]/df2[col]/np.sum(mask))[mask==True]))
                        print("MIN/MAXIMUM",np.nanmin((df1[col]/df2[col])[mask]),"|",np.max((df1[col]/df2[col])[mask]))
                        print("sample STD",np.std(np.array(df1[col]/df2[col])[mask],ddof=1))                        
                        
        damage_variable = compared_damage_measure
        dmg_measure = compared_damage_measure + " "
        if percentage_difference:
            dmg_measure += "fract diff"
        else:
            dmg_measure += "difference"
        if num_remaining_electrons_mode:
            dmg_measure+= " remaining"
        df_diff = df_diff.rename(columns={damage_variable:dmg_measure})
        df_diff.resolution = df1.resolution
        assert df_diff.resolution == df2.resolution
        name_of_set += "_"+compared_damage_measure 
        if percentage_difference:
            name_of_set+="_PERC"
        else:
            name_of_set+="_DIFF"
        print(df_diff)                                                 #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
        plot_that_funky_thing(df_diff,max_fwhm=max_fwhm,cmin=0,cmax=vmax,cmin_contour = 0, cmax_contour = vmax, contour_colour=contour_colour, contour_interval = contour_interval,name_of_set=name_of_set,**plot_2D_constants,template="plotly_dark",connect_contour_gaps=connect_contour_gaps,round_fluence=round_fluence,use_neutze_units = neutze,dmg_measure=dmg_measure,photon_min=photon_min) # 'electric' #"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
    print("Done")
# %%