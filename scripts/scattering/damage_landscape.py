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
from matplotlib.ticker import FormatStrFormatter
import pickle
import inspect
src_file_path = inspect.getfile(lambda: None)
my_dir = path.abspath(path.join(src_file_path ,"../")) + "/"
from core_functions import get_sim_params, get_pdb_paths_dict
from scatter import XFEL,Crystal,stylin,res_to_q,q_to_res



# full sim that does point sampling, thus assumes a "smooth" distribution. Thus good for low cell count input. Num cells is entire whole target.
# Compared with non-SPI (crystal) which assumes crystal and does bragg points. Num cells is supercell.
#TODO rename the diff sims, 'non-SPI' unclear.
#TODO currently all cells are used for calculation of R, but we may wish to ignore cells outside the considered radius? Though I haven't see anyone do this, it would depend on what image reconstruction algos do I guess. But we get empty corners with high res. limits so


MODE_DICT = {0:'crystal',1:'spi',2:'both'}

PDB_PATHS = get_pdb_paths_dict(my_dir) # <value:> the target that should be used for <key:> the name of simulation output batch folder.

DATA_FOLDER = "dmg_data/"
PLOT_FOLDER = "R_plots/"

def multi_damage(params,pdb_path,allowed_atoms_1,CNO_to_N,S_to_N,same_deviations,plasma_batch_handle = "", plasma_handles = None, sctr_results_batch_dir = None, get_R_only=True, realistic_crystal_growing_mode = False,specific_energy = None):
    '''
    NB: "Plasma" is used a shortened way to refer to the output of the damage simulation  
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
        exp_name1 = sim_handle + "_" + exp1_qualifier
        exp_name2 = sim_handle + "_" + exp2_qualifier

        # Spend a few billion on XFELs
        experiment1 = XFEL(exp_name1,energy[i],**run_params["experiment"])
        experiment2 = XFEL(exp_name2,energy[i],**run_params["experiment"])
        # Wait a few months for crystals to grow
        if realistic_crystal_growing_mode:
            # For added realism
            for day in range(90):
                print("Crystals growing,",90-day,"days remaining...")
                time.sleep(216000)
        else:
            print("Crystals growing at breakneck speeds...")
            time.sleep(0.3)            

        # The undamaged crystal form factor is constant (the initial state's) but still performs the same integration step with the pulse profile weighting.
        crystal = Crystal(pdb_path,allowed_atoms_1,is_damaged=True,CNO_to_N = CNO_to_N, S_to_N = S_to_N, **run_params["crystal"])
        if same_deviations:
            # we copy the other crystal so that it has the same deviations in coords
            crystal_undmged = copy.deepcopy(crystal)
            crystal_undmged.is_damaged = False
        else:
            crystal_undmged = Crystal(pdb_path,allowed_atoms_1,is_damaged=False,CNO_to_N = CNO_to_N,S_to_N = S_to_N, **run_params["crystal"])
        if i == 0:
            #crystal.plot_me(250000)
            crystal.plot_me(25000,template="plotly_dark")
    
        if run_params["laser"]["SPI"]:
            SPI_result1 = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal, **run_params["laser"])
            SPI_result2 = experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_undmged, **run_params["laser"])
            #TODO apply_background([SPI_result1,SPI_result2])
            damage_dict = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=get_R_only,SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
        else:
            exp1_orientations = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal, results_parent_dir=sctr_results_batch_dir, **run_params["laser"])
            if exp_name2 != None:
                # sync orientation with damage simulation
                experiment2.set_orientation_set(exp1_orientations)  
                run_params["laser"]["random_orientation"] = False 
                experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_undmged, results_parent_dir=sctr_results_batch_dir, **run_params["laser"])
            damage_dict = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=get_R_only,SPI=False,results_parent_dir = sctr_results_batch_dir)
        dmg_data.append(damage_dict) 
        pulse_params.append([energy[i],fwhm[i],photon_count[i]])
    
    return np.array(pulse_params,dtype=float),np.array(dmg_data,dtype=object), names, param_name_list

def save_data(fname, data):
    os.makedirs(DATA_FOLDER,exist_ok=True)

    with open(DATA_FOLDER+fname+".pickle", 'wb') as file: 
        pickle.dump(data, file)

    
def load_df(fname, damage_measure,bin_limit,check_batch_nums=True):
    with open(DATA_FOLDER+fname +".pickle", "rb") as file:
        pulse_params, dmg_data,names,param_name_list = pickle.load(file) 
    if damage_measure == "bin_q_R" or damage_measure == "bin_q_I":
        unique_bins = []
        dmg_vect = []
        for i, dmg_dict in enumerate(dmg_data):
            bins = []
            dmg = []
            for k,v in dmg_dict[damage_measure].items():
                bins.append(k)
                dmg.append(v)
            bins,dmg = zip(*(sorted(zip(bins,dmg), key=lambda x: x[0])))
            bin_idx = np.searchsorted(bins,bin_limit)
            if bins[bin_idx] not in unique_bins:
                unique_bins.append(bins[bin_idx])
            dmg_vect.append(dmg[bin_idx]) 
        assert(np.max(unique_bins) - np.min(unique_bins)) < np.max(unique_bins)*0.025 , "nearest q  above the q specified for each simulation differ by a value above the tolerance - unique values: " + str(unique_bins) + " Try a different value (sorry for this bad coding)." 
        average_bin_limit = np.average(unique_bins)
    else:
        resolutions = np.array([d["resolutions"] for d in dmg_data])
        resolutions_idxes = np.empty(resolutions.shape)
        unique_resolutions = []
        dmg_vect = []
        for i, vect in enumerate(resolutions):
            res_idx = np.searchsorted(vect,bin_limit)
            if vect[res_idx] not in unique_resolutions:
                unique_resolutions.append(vect[res_idx])
            dmg_vect.append(dmg_data[i][damage_measure][res_idx]) 
        print("RESOLUTIONS USED",unique_resolutions)
        # Sorry for this
        assert(np.max(unique_resolutions) - np.min(unique_resolutions)) < 0.05 , "nearest resolution above the resolution specified for each simulation differ by a value above the tolerance - unique values: " + str(unique_resolutions) + " Try a different value (sorry for this bad coding)." 
        average_bin_limit = np.average(unique_resolutions)
    dmg_data = np.array(dmg_vect,dtype=float)
    
    
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
        "energy": pulse_params[:,0],
        "fwhm": pulse_params[:,1],
        photon_measure: photon_data,
        damage_measure: dmg_data,
        "damage_measure": damage_measure,
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
    df.average_bin_limit = average_bin_limit
    return df   
    
def plot_that_funky_thing(df,cmin=0.1,cmax=0.3,clr_scale="temps",use_neutze_units=False,name_of_set="",energy_key=9000,photon_key=1e14,fwhm_key=50,cmin_contour=0,cmax_contour=None,contour_colour = "amp",contour_interval=0.05,dmg_measure="R",photon_min=None,**kwargs,):
    out_folder = PLOT_FOLDER
    os.makedirs(out_folder,exist_ok=True)
    ext = ".svg"
    photon_measure = df.columns[3]

    if photon_measure[:5] == "Count":
        photon_measure_default = "ph per um^2"
        if use_neutze_units:
            df[photon_measure] = [7.854e-3*p for p in df[photon_measure]] # square micrometre to 100 nm diameter spot
            photon_measure = "photons per 100 nm spot"
        else:
            photon_measure = photon_measure_default
        df.rename(columns=dict(Count=photon_measure),inplace=True)
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
        fig.write_image(out_folder+fig["layout"]["title"]["text"]+ext)

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
    fig.write_image(out_folder+fig["layout"]["title"]["text"]+ext)


    ranges = {}
    for col in df:
        if col in ["name","damage_measure"]:
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
        connectgaps = False,
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
        fig.write_image(out_folder+fig["layout"]["title"]["text"]+"-Sctr"+ext)
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
        fig.update_layout(width = 750, height = 600,)
        fig.update_layout(title= name_of_set + ", E = "+str(energy_key)+" eV")
        fig.show()
        fig.write_image(out_folder+fig["layout"]["title"]["text"]+ext)
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
        #print(df)

        fig = go.Figure(data =  
            go.Contour(
                #z=R_mesh,
                #x=unique_fwhm, 
                #y=unique_photons,
                z = df[dmg_measure],
                x = df["fwhm"],
                y = df["energy"],
                **contour_args,
            ))
        fig.update_layout(width = 750, height = 600)
        fig.update_layout(title = name_of_set + ", " +str(photon_key_for_title)+" " +photon_measure)
        fig.update_xaxes(title="FWHM (fs)",type="log")
        fig.update_yaxes(title="Energy (eV)")

        fig.show()   
        fig.write_image(out_folder+name_of_set + ", " +str(default_unit_photon_key)+" " +photon_measure_default+ext)
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
                y = df["energy"],
                x = [p for p in df[photon_measure]],
               **contour_args,
            ))
        fig.update_layout(width = 750, height = 600,)
        fig.update_layout(title=name_of_set + ", FWHM = "+str(fwhm_key)+" fs")
        fig.update_yaxes(title="Energy (eV)")
        fig.update_xaxes(title=photon_measure,type="log",range=np.log10(ranges[photon_measure]),tickvals = [1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15],tickformat = '.0e')
        fig.show()    
        fig.write_image(out_folder+fig["layout"]["title"]["text"]+ext)    
        df = original_df
#%%
#------------Generate data (get R)--------------------

if __name__ == "__main__":
    # Initialise kwargs
    kwargs = dict(
        plasma_batch_handle = "", 
        plasma_handles = None, 
        sctr_results_batch_dir = None, 
        get_R_only=True,
    )


    batch_mode = False # Just doing this as I want to quickly switch between doing batches and comparing specific runs.

    mode = 1  #0 -> infinite crystal, 1 -> finite crystal/SPI, 2-> both  
    same_deviations = True # whether same position deviations between damaged and undamaged crystal (SPI only)
    
    imaging_params_preset = copy.deepcopy(imaging_params.default_dict)

    # Change options here... #TODO make clearer this is for user options
    if batch_mode:
        allowed_atoms = ["C","N","O","S"] 
        CNO_to_N = False
        S_to_N = True
        batch_handle = "lys_all_light" 
        batch_dir = None # Optional: Specify existing parent folder for batch of results, to add these orientation results to.
        pdb_path = PDB_PATHS["lys"]
        # Params set to batch_handle
        kwargs["plasma_batch_handle"] = batch_handle
        fname = batch_handle
    else: # Compare specific simulations.  Resolution map code intended to be used with this.
        #imaging_params_preset = copy.deepcopy(imaging_params.asymmetric_unit_dict)
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
        #kwargs["plasma_handles"] = ["lys_nass_no_S_3","lys_nass_gauss","lys_nass_Gd_full_1"]
        #kwargs["plasma_handles"] = ["lys_nass_gauss","lys_nass_Gd_full_1","lys_nass_gauss_solvated_1","lys_nass_Gd_gauss_solvated_2"]
        #kwargs["plasma_handles"] = ["lys_nass_gauss_dry_9kev_1","lys_nass_gauss_dry_1","lys_nass_light_dry_9kev_1","lys_nass_light_dry_1","lys_nass_Gd_gauss_dry_1","lys_nass_Gd_gauss_dry_9kev_1"]
        kwargs["plasma_handles"] = ["lys_nass_water_solvent_2","lys_nass_water_solvent_9kev_2","lys_nass_gd_solvent_5","lys_nass_gd_solvent_9kev_5"]
        #kwargs["plasma_handles"] = ["lys_nass_water_solvent_1","lys_nass_water_solvent_9kev_1","lys_nass_gd_solvent_2","lys_nass_gd_solvent_9kev_1"]
        #kwargs["plasma_handles"] = ["lys_nass_no_S_3","lys_nass_Gd_full_1"]
        #kwargs["plasma_handles"] = ["lys_nass_HF","lys_nass_Gd_HF"]  
        #kwargs["plasma_handles"] = ["lys_full-typical","lys_all_light-typical"]  
        #kwargs["plasma_handles"] = ["glycine_abdullah_4"]
        #pdb_path = PDB_PATHS["fcc"]
        pdb_path = PDB_PATHS["lys"]
        #pdb_path = PDB_PATHS["lys_solvated"]
        #pdb_path = PDB_PATHS["glycine"]
        #pdb_path = PDB_PATHS["tetra"]

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
    #data_name = "lys_full"; batch_mode = True; mode = 1  #TODO store batch_mode and mode in saved object.
    data_name = "comparison"; batch_mode = False; mode = 1  #TODO store batch_mode and mode in saved object.
    #####
    name_of_set = data_name
    _damage_measure = "R"
    resolution = 2 #1.9 2.8
    df = load_df(data_name, damage_measure=_damage_measure, bin_limit=resolution, check_batch_nums=batch_mode) # resolution list is orders from best to worst resolutions.
    plot_2D_constants = dict(energy_key = 15000, photon_key = 1e12,fwhm_key = 25)
    neutze = False
    name_of_set += _damage_measure+"lim"+str("{:.1f}".format(df.average_bin_limit)) # Currently there's an uncertainty in the 2nd decimal because I coded bad, so 2 Sig figs it shall be.
    if MODE_DICT[mode] != "spi":
        print("-----------------Crystal----------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
        plot_that_funky_thing(df,0,0.20,"temps",name_of_set=name_of_set,**plot_2D_constants,template="plotly_dark",use_neutze_units = neutze,cmax_contour=0.4) # 'electric' #"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
    if MODE_DICT[mode] != "crystal":
        print("-------------------SPI------------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
        plot_that_funky_thing(df,0,0.20,"temps",name_of_set=name_of_set,**plot_2D_constants,template="plotly_dark",use_neutze_units=neutze,cmax_contour=0.4) # 'electric'#"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
    print("Done")

#%%
# ------------- RESOLUTION MAP ------------ #
if __name__ == "__main__":
    data_name = "comparison"; batch_mode = False; mode = 1
    #####
    with open(DATA_FOLDER+data_name +".pickle", "rb") as file:
        _, dmg_data,_,_ = pickle.load(file) 
    
    #_damage_measure = "R"
    #_damage_measure = "bin_res_R"
    #_damage_measure = "bin_q_R"
    #_damage_measure = "cc"
    _damage_measure = "bin_q_I"
    dmg = []
    bins = []
    if _damage_measure in ["bin_q_R","bin_q_I"]:
        for elem in dmg_data:
            tmp_X = []
            for k, v in elem[_damage_measure].items():
                tmp_X.append(k)
            bins.append(tmp_X)
            assert(tmp_X == bins[0])       
        
    else:
        bins = np.array([d["resolutions"] for d in dmg_data])   # TODO Use dict for this same as above
    names = []
    # Iterate through different bins
    for bin_value in bins[0]:
        df = load_df(data_name, _damage_measure,bin_value, check_batch_nums=batch_mode) # resolution index 0 corresponds to the max resolution
        dmg.append(df[_damage_measure])
    dmg = np.array(dmg)
    names = np.array(df["name"])
    for i,name in enumerate(names):
        assert np.array_equal(bins[i],bins[0])      # Assert that the sampled X points are identical between sims.
    X = np.array(bins[0])

    print(dmg)
    print( np.array(["nan"]*len(dmg[-1])))
    # while np.isnan(dmg[-1,:]).all():
    #     X.pop()
    #     dmg = dmg[:-1,:]

    cmap = "viridis"
    y_label = ""
    if _damage_measure == "bin_res_R":
        y_label = "res. ring $R_{dmg}$"
        vmax = 0.1
    if _damage_measure == "bin_q_R":
        y_label = "$q$ ring $R_{dmg}$"
        vmax = 0.1      
    if _damage_measure == "R":
        y_label = "$R_{dmg}$"                  
        vmax = 0.1
    if _damage_measure == "cc":
        y_label = "cc"
        vmax = 0.1  
    if _damage_measure == "bin_q_I":
        y_label = "log(mean ring ratio)"
        vmin = -0.4
        vmax = 0.4
        cmap = "seismic"

    # Linear
    fig,ax = plt.subplots()
    for i,name in enumerate(names):
        x = copy.deepcopy(X)
        if _damage_measure in ["bin_q_R","bin_q_I"]:
            x = q_to_res(x)
        ax.scatter(x,dmg[:,i],label=name)

    labels=["Lysozyme, no S","Lysozyme","Lysozyme.Gd"]
    #labels=["Gaussian pulse","Square pulse",]
    ax.legend(labels=labels,reverse=True)
    ax.set_ylabel(y_label)
    ax.set_ylim([None,0.105])
    ax.set_xlabel("Resolution ($\AA$)")
    ax.set_xscale("log")
    if _damage_measure in ["bin_q_R","bin_q_I"]:
        ax.set_xlim([np.min(x),10])
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
    else:
        ax.set_xlim([np.min(x),np.max(x)])
    plt.show()
    #Radial note/TODO uses small angle approximation (plots against q rather than actual spatial distance.)
    for i,name in enumerate(names):
        #plt.plot(sim_resolutions[0],R_factors[:,i],label=name)
        fig, ax = plt.subplots(subplot_kw={'projection':'polar'})
        dmg_binned = dmg
        if _damage_measure in ["R","cc","bin_res_R"]:
            q = np.array([res_to_q(r) for r in X])
        if _damage_measure in ["bin_q_R","bin_q_I"]:
            q = np.array([r for r in X])
            
            # Make the bloody centre point actually get filled. Pain for forgetting numpy intricacies
            q = np.array(list(q) + [-np.min(q)])
            ###dmg = [list(d) for d in dmg] + [[0]*len(names)]
            ###dmg = np.array(dmg)
            dmg_binned = np.append(dmg,np.array([np.array([0]*dmg.shape[-1])]),axis=0)
            q, dmg_binned = zip(*(sorted(zip(q,dmg_binned), key=lambda x: x[0])))
            q = np.array(q)
            dmg_binned = np.array(dmg_binned)
        azm = np.linspace(0, 2 * np.pi, 100)
        r, th = np.meshgrid(q, azm)
        dmg_meshed = (r*0+1)*dmg_binned[None,:,i]
        cbar = ax.pcolormesh(th,r,dmg_meshed,vmin=vmin,vmax=vmax,cmap=cmap)          
        cbar = plt.colorbar(cbar,shrink=0.8)
        cbar.set_label(y_label,rotation = 0, y = 1.1,horizontalalignment='right')
        plt.show()    
    print("Done")


#%%
'''
# Surgery
fname = "lys_all_light"
with open(DATA_FOLDER+fname +".pickle", "rb") as file:
    pulse_params, dmg_data,_,names,param_name_list = pickle.load(file) 
with open(DATA_FOLDER+"resolutions.pickle","rb") as file:
    _, _,resolutions,_,_ = pickle.load(file) 
save_data(fname +"Res", (pulse_params,dmg_data, resolutions, names, param_name_list))
'''

#%%
#------------Compare----------------------
contour_interval_delta_difference = 0.01 # 0.01
contour_interval_percentage_difference = 0.05
if __name__ == "__main__":
    batch_mode = True; mode = 1
    
    resolution = 2.45#2.2
    photon_min=1e10
    fname1 = "lys_full"; fname2 = "lys_all_light"
    plot_2D_constants = dict(energy_key = 15000, photon_key = 1e12,fwhm_key = 15)

    statistics_calculation_R_cutoff = 0.1
    #for percentage_difference in [False,True]:
    for percentage_difference in [True]:
        if percentage_difference:
            contour_interval = contour_interval_percentage_difference
            contour_colour = "viridis"
            vmax =  0.4
            mask_below_value = 0.1
        else:
            contour_interval = contour_interval_delta_difference
            contour_colour = "plasma"
            vmax = 0.1

        #name_of_set = "Difference" + "_"+fname1 + "_" + fname2
        name_of_set = "Sulfur_impact2"
        df1 =  load_df(fname1, resolution, check_batch_nums=batch_mode)
        df2 = load_df(fname2, resolution, check_batch_nums=batch_mode)
        neutze = False
        df_diff = copy.deepcopy(df1)
        for col in df1.columns:
            if col != "R" and col != "cc" and col!="name":
                assert np.array_equal(df_diff[col], df2[col]) and np.array_equal(df1[col], df_diff[col])
            elif col != "name":
                if percentage_difference:
                    df_diff[col] = df1[col]/df2[col] - 1
                    #mask = df1[col] < statistics_calculation_R_cutoff
                    #df_diff[col][mask] = None 
                else:
                    df_diff[col] -= df2[col]
                    if col == "R":
                        num_non_pos_count = np.sum(df_diff[col] <= 0)
                        print("NUMBER OF R_diff <= 0:", num_non_pos_count)  
                        print("AVERAGE FACTOR",np.sum(df1[col]/df2[col])/len(df1[col]))
                        print("MIN/MAXIMUM FACTOR",np.min(df1[col]/df2[col]),"|",np.max(df1[col]/df2[col]))
                        print("sample STD",np.std(df1[col]/df2[col],ddof=1))
                        print("==========")
                        mask = np.array(df1[col] > statistics_calculation_R_cutoff,dtype=bool)
                        print(np.sum(mask))
                        print("VALUES WHERE R FACTOR OF TARGET1 ABOVE", statistics_calculation_R_cutoff,":")
                        print("AVERAGE FACTOR",np.nansum(np.array(df1[col]/df2[col]/np.sum(mask))[mask==True]))
                        print("MIN/MAXIMUM FACTOR",np.nanmin((df1[col]/df2[col])[mask]),"|",np.max((df1[col]/df2[col])[mask]))
                        print("sample STD",np.std(np.array(df1[col]/df2[col])[mask],ddof=1))                        
                        
        if percentage_difference:
            dmg_measure = "R1/R2-1"
        else:
            dmg_measure = "R1- R2"
        df_diff = df_diff.rename(columns={"R":dmg_measure})
        df_diff.resolution = df1.resolution
        assert df_diff.resolution == df2.resolution
        name_of_set += "-res="+str("{:.1f}".format(df_diff.resolution)) # Currently there's an uncertainty in the 2nd decimal because I coded bad, so 2 Sig figs it shall be.
        if percentage_difference:
            name_of_set+="_PERC"
        else:
            name_of_set+="_DIFF"
        #print(df_diff)                                                 #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
        plot_that_funky_thing(df_diff,cmin=0,cmax=vmax,cmin_contour = 0, cmax_contour = vmax, contour_colour=contour_colour, contour_interval = contour_interval,name_of_set=name_of_set,**plot_2D_constants,template="plotly_dark",use_neutze_units = neutze,dmg_measure=dmg_measure,photon_min=photon_min) # 'electric' #"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
    print("Done")
# %%