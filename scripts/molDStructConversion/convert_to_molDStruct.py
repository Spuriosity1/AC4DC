#%%
import sys
import os
import pandas as pd
import os.path as path
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/scattering')
from scatter import XFEL,Crystal,stylin
from core_functions import get_sim_params,get_sim_elements,get_pdb_path,ATOMNO
import imaging_params as imaging_params
import numpy as np
from scipy import constants as C
import struct

# Converts AC4DC data to IONIZATION_DATA used for input to MolDStruct CR-MD.



crystal_params = dict(
    supercell_scale = 1,  ##3 # for SC: cell_scale^3 unit cells 
    num_supercells = 1,
    supercell_simulations = 1,        
    include_symmetries = None, ##True  # should unit cell contain symmetries?
    positional_stdv = 0, # Introduces disorder to positions. Note this is a deviation from the IDEAL structure, so is not a measure of similarity with undamaged and damaged structure but the ideal structure to recover and the dmaaged structure. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if quickly gauging serial crystallography R factor, as should somewhat average out.
    cell_packing = "SC",
    rocking_angle = 1.2,  # (approximating mosaicity, infinite crystal sim only)
)


def get_random_charge_states(element):
    SEEDED = False
    element.times_used = element.crystal.ff_calculator.get_times_used()
    if element.num_atoms != len(element.crystal.sym_rotations)*len(element.coords):
        raise Exception("num atoms was not same on set_stochastic_states call as when set by set_coord_deviation")
    if element.crystal.is_damaged:
        charges = np.empty(shape = (element.num_atoms,len(element.times_used)))  # (num atoms, times)   
        for idx in range(element.num_atoms):
            seed = None
            if SEEDED:
                seed = idx
            charges[idx] = element.crystal.ff_calculator.random_charge_snapshots(element.name,seed) 
            #assert np.all(charges.astype(np.ushort)[idx] <= ATOMNO[element.name])

    return charges

def binary(num):
    return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))

def create_charge_file(charges,element,save_dir,overwrite=False,csv=False):
    '''
    Generates a

    '''
    save_dir_csv = save_dir
    save_dir_bin = save_dir + "IONIZATION_DATA/"
    print(f"Creating {'full' if element is None else element} charge file for {save_dir.split('/')[-1]}")
    os.makedirs(save_dir_bin, exist_ok=True) 
    os.makedirs(save_dir_csv, exist_ok=True) 
    if element is None:
        save_path = f"{save_dir_bin}charges"
        save_path_csv = f"{save_dir_csv}charges"
    else:
        save_path = f"{save_dir_bin}{element}_charges"
        save_path_csv = f"{save_dir}{element}_charges"
    
    if os.path.isfile(save_dir): 
        if not overwrite:
            print("Cannot write, file already present at",save_dir)
            return
        os.remove(save_dir)

    # with open(save_path+".bin", "wb") as file:
    #     newFileByteArray = bytearray(charges)
    #     file.write(newFileByteArray)
    charges.astype(np.ushort).swapaxes(0,1).tofile(save_path+".bin") #uint16, shape = (num timesteps, num atoms)
    
    if csv:
        df = pd.DataFrame(charges)
        df.to_csv(save_path_csv+".csv", index=False)

    with open(save_path+".bin", "rb") as file:
        n=10
        print(f"First {n} binary data elements:")
        print(struct.unpack('H'*n, file.read(2*n)))


def create_data_file(data,tag,save_dir,overwrite=False,csv=False):
    '''
    Generates a
    '''
    save_dir_csv = save_dir
    save_dir_bin = save_dir + "IONIZATION_DATA/"
    
    data[0] = data[1]  # Patch.

    print(f"Creating {tag} file in {save_dir.split('/')[-1]}")
    os.makedirs(save_dir_bin, exist_ok=True) 
    os.makedirs(save_dir_csv, exist_ok=True) 
    save_path = save_dir+"IONIZATION_DATA/"+tag
    save_path_csv = save_dir+tag
    if os.path.isfile(save_dir): 
        if not overwrite:
            print("Cannot write, file already present at",save_dir)
            return
        os.remove(save_dir)

    # with open(save_path+".bin", "wb") as file:
    #     newFileByteArray = bytearray(data)
    #     file.write(newFileByteArray)
    data.astype('float32').tofile(save_path+".bin")
    
    if csv:
        df = pd.DataFrame(data)
        df.to_csv(save_path_csv+".csv", index=False)

    with open(save_path+".bin", "rb") as file:
        n=10
        print(f"First {n} binary data elements:")
        print(struct.unpack('f'*n, file.read(4*n)))

#PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"I3C") 

SCRIPTS_DIR = path.abspath(path.join(__file__ ,"../../")) + "/"

SCATTER_DIR = SCRIPTS_DIR +"scattering/"
#PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"tetra") 
OUTPUT_PATH =  path.abspath(path.join(__file__ ,"../"))+ "/output/"

TARGET_DIR = SCATTER_DIR+ "targets/"


MOLECULAR_PATH = path.abspath(path.join(SCRIPTS_DIR, "../output/__Molecular/")) + "/" # directory of damage sim output folders

SAVE_FOLDER = "test"


def charges(csv=False,individual_elements = False): 
    print("Beginning writing of charges...")

    out_folder = OUTPUT_PATH + SAVE_FOLDER + "/"
    
    num_steps = len(ff_calculator.get_times_used())
    species_charges = {}
    num_atoms = 0
    print(f"Processing charges at {num_steps} time steps for:")
    for element in crystal.species_dict.keys():
        print(f"{len(crystal.species_dict[element].serial_numbers)} {element} atoms")

    for element in crystal.species_dict.keys():
        print(f"{element}...")
        species_charges[element] = get_random_charge_states(crystal.species_dict[element])
        num_atoms += len(species_charges[element])
    
 
    # Order charges in order that matches the structure file. 
    combined_charges = np.empty(shape = (num_atoms,num_steps))  # (num atoms, times)   
    species_list = np.empty(shape = (num_atoms,),dtype=object)
    for element, charges in species_charges.items():
        for i, s_num in enumerate(crystal.species_dict[element].serial_numbers):
            combined_charges[s_num-1] = charges[i]
            species_list[s_num-1] = element

        if individual_elements:
            PDB_element = element.split("_")[0]
            create_charge_file(charges,PDB_element,out_folder,csv=csv)    # Each row is an atom. each column is a time step.
    create_charge_file(combined_charges,None,out_folder,csv=csv)    # Each row is an atom. each column is a time step.



def DebyeLength(csv=False):
    pl = ff_calculator
    

    tempList = []
    denseList = []
    for t in ff_calculator.get_times_used():
        tempList.append( pl.get_temp(t, 1000) ) # eV
        # denseList.append( pl.get_density(t) ) # per angstrom cube

    T = np.array(tempList)
    n = T [: ,1]
    T = T[:,0]
    print("n:",n)
    print("T:",T)

    # Can ignore if at step 0
    for i, elem in enumerate(T):
        if elem <= 0:
            print(f"Warning, T at step {i} = {elem}")
    for i, elem in enumerate(n):
        if elem <= 0:
            print(f"Warning, n at step {i} = {elem}")

    lambdaD=np.sqrt(C.epsilon_0 * C.nano * T *C.eV / n /C.e/C.e) # should have units nm
    #lambdaD=np.sqrt(C.epsilon_0 * C.angstrom * T *C.eV / n /C.e/C.e) # should have units Angstrom

    out_folder = OUTPUT_PATH + SAVE_FOLDER + "/"
    create_data_file(T*11606,"electron_temperature",out_folder,csv=csv) # K
    create_data_file(n/C.nano**3,"electron_density",out_folder,csv=csv) # nm^-3
    create_data_file(lambdaD,"debye_data",out_folder,csv=csv) # nm


sim_handle = "lys_nass_gauss_solvated_with_H_21"
num_steps = 3800
allowed_atoms = get_sim_elements(sim_handle)

target = "lys_conf.gro"




crystal = Crystal(TARGET_DIR + target,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=False,**crystal_params)



param_dict,_,_ = get_sim_params(sim_handle)
start_time = param_dict["start_t"]
end_time = param_dict["end_t"]
energy = param_dict["energy"]



# Assign plotter object to calculate charges (because code debt)
xfel = XFEL("dummy",energy,t_fineness=num_steps)
ff_calculator = xfel.get_ff_calculator(start_time,end_time,sim_handle,MOLECULAR_PATH)   
ff_calculator.allow_select_same_times = True  
crystal.set_ff_calculator(ff_calculator)    



DebyeLength(csv=True)
charges(csv=True)

# class Target(Enum):
#     UNIT = imaging_params.goldilocks_dict_unit
#     NINE = imaging_params.goldilocks_dict_3x3x3


#im_params = Target.NINE
#crystal = Crystal(PDB_STRUCTURE,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=True, **im_params["crystal"])

# %%






