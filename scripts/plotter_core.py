import stringprep
import matplotlib.rcsetup as rcsetup
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
# from scipy.interpolate import BSpline
from math import log
import os.path as path
import os
import matplotlib.colors as colors
import sys
# import glob
import csv
import subprocess
from matplotlib.ticker import LogFormatter 
import random
from scipy.optimize import curve_fit
from scipy.stats import linregress
from core_functions import get_mol_file, parse_elecs_from_latex, ATOMS, ATOMNO
from scipy.interpolate import splrep, splev
from scipy.signal import savgol_filter

#plt.rcParams.update(plt.rcParamsDefault)
#plt.style.use('seaborn-muted')
#plt.style.use('turrell-style')


def get_colors(num, seed):
    idx = list(np.linspace(0, 1, num))[1:]
    random.seed(seed)
    # random.shuffle(idx)
    idx.insert(0,0)
    C = plt.get_cmap('nipy_spectral')
    return C(idx)

class Plotter:
    # Example initialisation: Plotter(water,molecular/path)
    # --> Data is contained in molecular/path/water. 
    # Will use mol file within by default, or (with a warning) search input for matching name if none exists.  
    def __init__(self, data_folder_name, abs_molecular_path = None, num_subplots=None,use_electron_density = False,out_prefix_text = None,end_t=None,load_specific_atoms=None,sample_end_points_only=False):
        '''
        abs_molecular_path: The path to the folder containing the simulation output folder of interest.
        use_electron_density: If True, plot electron density rather than energy density
        '''
        self.end_t_plotting = end_t
        self.sample_end_points = sample_end_points_only
        self.use_electron_density = use_electron_density
        self.molecular_path = abs_molecular_path
        if self.molecular_path is None:
            self.molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
        AC4DC_dir = path.abspath(path.join(__file__ ,"../../"))  + "/"
        self.input_path = AC4DC_dir + 'input/'
        if out_prefix_text is None:
            out_prefix_text = "Initialising plotting with"
        molfile = get_mol_file(self.input_path,self.molecular_path,data_folder_name,"y",out_prefix_text = out_prefix_text) 

        self.mol = {'name': data_folder_name, 'infile': molfile, 'mtime': path.getmtime(molfile)}        

        # Stores the atomic input files read by ac4dc
        self.atomdict = {}
        self.statedict = {}

        # Outputs
        self.outDir = self.molecular_path + data_folder_name
        self.freeFile = self.outDir +"/freeDist.csv"
        self.freeFiles = []  # Distributions for each element's cascdes.
        self.intFile = self.outDir + "/intensity.csv"
        self.gridFile = self.outDir + "/knotHistory.csv"

        self.boundData={}
        self.photoData={}
        self.chargeData={}
        self.freeData=None
        self.intensityData=None
        self.energyKnot=None
        self.timeData=None

        self.get_atoms(load_specific_atoms)
        self.update_outputs()
        self.autorun=False
        if num_subplots is not None:
            self.setup_axes(num_subplots)
    
    def find_mol_file_from_directory(self, input_directory, mol):
        # Get molfile from all subdirectories in input folder.
        molfname_candidates = []
        for dirpath, dirnames, fnames in os.walk(input_directory):
            #if not "input/" in dirpath: continue
            for molfname in [f for f in fnames if f == mol+".mol"]:
                molfname_candidates.append(path.join(dirpath, molfname))
        # No file found
        if len(molfname_candidates) == 0:
            print('\033[91m[ Error: File not found ]\033[0m ' + "No '"+ mol + ".mol' input file found in input folders, trying default path anyway.")
            molfile = input_directory+mol+".mol"
        elif len(molfname_candidates) == 1:
            molfile = molfname_candidates[0]
        # File found in multiple directories
        else:
            print('\033[95m' + "Multiple mol files with given name detected. Please input number corresponding to desired directory." + '\033[0m' )
            for idx, val in enumerate(molfname_candidates):
                print(idx,val)
            selected_file = False
            # Loop until user confirms file
            while selected_file == False: 
                molfile = molfname_candidates[int(input("Input directory number: "))]
                print('\033[95m' + molfile + " selected." + '\033[0m')
                y_n_check = input("Input 'y'/'n' to continue/select a different file: ")
                if y_n_check.casefold() not in map(str.casefold,["y","yes"]):
                    if y_n_check.casefold() not in map(str.casefold,["n","no"]):
                        print("Unknown response, using 'n'.")
                    continue
                else:
                    print("Continuing...")
                    selected_file = True
        return molfile    
        

    # Reads the control file specified by self.mol['infile']
    # and populates the atomdict data structure accordingly
    def get_atoms(self,atoms_to_load = None):
        self.atomdict = {}
        with open(self.mol['infile'], 'r') as f:
            reading = False
            for line in f:
                if line.startswith("#ATOMS"):
                    reading=True
                    continue
                elif line.startswith("#") or line.startswith("//"):
                    reading=False
                    continue
                if line.startswith("####END####"):
                    break
                if reading:
                    a = line.split(' ')[0].strip()
                    if len(a) != 0 and (atoms_to_load is None or a in atoms_to_load):
                        file = self.input_path + 'atoms/' + a + '.inp'
                        self.atomdict[a]={
                            'infile': file,
                            'mtime': path.getmtime(file),
                            'outfile': self.outDir+"/dist_%s.csv"%a,
                            'photofile': self.outDir + "/photo_%s.csv"%a}
                        if atoms_to_load is not None:
                            self.freeFiles.append(self.outDir +"/freeDist_"+a+".csv")
        if atoms_to_load is None:
            self.freeFiles.append(self.freeFile)
    def update_inputs(self,load_specific_atoms=None):
        self.get_atoms(load_specific_atoms)
        self.mol['mtime'] = path.getmtime(self.mol['infile'])

    # def rerun_ac4dc(self):
    #     cmd = self.p+'/bin/ac4dc2 '+self.mol['infile']
    #     print("Running: ", cmd)
    #     subprocess.run(cmd, shell=True, check=True)

    def check_current(self):
        # Pull most recent atom mod time
        self.update_inputs()
        # Outputs are made at roughly the same time: Take the oldest.
        out_time = path.getmtime(self.freeFile)
        for atomic in self.atomdict.values():
            out_time = min(out_time, path.getmtime(atomic['outfile']))

        if self.mol['mtime'] > out_time:
            print("Master file %s is newer than most recent run" % self.mol['infile'])
            return False

        for atomic in self.atomdict.values():
            if atomic['mtime'] > out_time:
                print("Dependent input file %s is newer than most recent run" % atomic['infile'])
                return False
        return True

    def get_free_energy_spec(self):
        erow = []
        with open(self.freeFile) as f:
            r = csv.reader(f, delimiter=' ')
            for row in r:
                if row[0] != '#':
                    raise Exception('parser could not find energy grid specification, expected #  |')
                reading=False
                for entry in row:
                    # skip whitespace
                    if entry == '' or entry == '#':
                        continue
                    # it's gotta be a pipe or nothin'
                    if reading:
                        erow.append(entry)
                    elif entry == '|':
                        reading=True
                    else:
                        break
                if reading:
                    break
        return erow

    def get_bound_config_spec(self, a):
        # gets the configuration strings corresponding to atom a
        specs = []
        with open(self.atomdict[a]['outfile']) as f:
            r = csv.reader(f, delimiter=' ')
            for row in r:
                if row[0] != '#':
                    raise Exception('parser could not find atomic state specification, expected #  |')
                reading=False
                for entry in row:
                    # skip whitespace
                    if entry == '' or entry == '#':
                        continue
                    # it's gotta be a pipe or nothin'
                    if reading:
                        specs.append(entry)
                    elif entry == '|':
                        reading=True
                    else:
                        break
                if reading:
                    break
        return specs


    def get_atomic_numbers(self):
        # Get atomic number by considering states at time 0.
        return_dict = {}
        for a in self.atomdict:
            states = self.statedict[a]  
            atomic_number = 0
            for val in parse_elecs_from_latex(states[0]).values():
                atomic_number += val
            return_dict[a] = atomic_number
        return return_dict
    


    def initialise_form_factor_params(self, start_t,end_t, q_max, photon_energy, q_fineness=50,t_fineness=50,naive=False):
        args = locals()
        self.__dict__ .update(args)  # Is this a coding sin?  ...yes.
 
    
    # q in units of bohr^-1
    def theta_to_q(self,scattering_angle_deg,photon_energy=None):
        if photon_energy == None:
            photon_energy = self.photon_energy
        theta = scattering_angle_deg*np.pi/180
        bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h = 6.582119569e-16 *2 *np.pi
        wavelength =  h*c_nm/photon_energy
        wavelength *= bohr_per_nm
        return 4*np.pi*np.sin(theta)/wavelength
    # q in units of bohr^-1
    def q_to_theta(self,q,photon_energy=None):
        if photon_energy == None:
            photon_energy = self.photon_energy        
        bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h = 6.582119569e-16 *2 *np.pi
        wavelength =  h*c_nm/photon_energy   
        wavelength *= bohr_per_nm
        return 180/np.pi*np.arcsin(q*wavelength/(4*np.pi))     


    
    def plot_R_factor(self, atoms_subset=None):   #TODO Would be best to integrate with plot_ffactor_get_R_AC4DC since that uses the simulation's calculated struct. factors which are presumably more accurate.
        atoms = atoms_subset
        if atoms == None:
            atoms = self.atomdict

        # Plot R over different k ranges, up to k_max which corresponds to the best resolution but should have the lowest R-factor.

        q_points = np.linspace(0,self.q_max,self.q_fineness-1)  #TODO make better and don't try and do in head...     

        R_factor = []
        for q in q_points:
           self.get_average_form_factor() 

        F_i = []
        F_r = []
        for i,k in enumerate(k_points):
            F_i.append(np.sqrt(I_ideal[i][i]))
            F_r.append(np.sqrt(I_real[i][i]))
        print("real",F_r)
        print("ideal",F_i)

      
        
        #q_resolution = np.array([[3.2,1],[0.032,100]])  #TODO   
        q_theta = []
        for i, q in enumerate(q_points):
            if i%(int(self.q_fineness/5)) == 1:
                q_theta.append([q,self.q_to_theta(q)])
        q_theta.reverse()
        q_theta = np.array(q_theta)
        
        dual_axis = q_theta
        #top_label = r"Resolution (\AA)"
        top_label = r"Theta"         
        
        fig_size = 5
        self.fig = plt.figure(figsize=(fig_size,fig_size))
        ax = self.fig.add_subplot(111)
        ax.set_xlabel("q length",fontsize = 12)
        ax.set_xscale("log")
        ax.set_xlim(dual_axis[0][0],dual_axis[-1][0])
        ax.set_ylim(0,1) 
        #ax.set_xlim(0.8,0.2)
        #ax.set_ylim(0.5,0.7) 
        ax.set_ylabel("R")
        ax.plot(q_points,R_factor)

        # Resolution or angle axis 
        top_tick_locations = dual_axis[:,0]
        top_tick_labels = dual_axis[:,1]    
        top_tick_labels = ['%.2f' % l for l in top_tick_labels]
         
        ax2 = ax.twiny()              
        #ax2.set_xlim(q_resolution[1][1],q_resolution[0][1])
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticklabels(top_tick_labels)
        ax2.set_xticks(top_tick_locations)
        ax2.set_xlabel(top_label)
        

    def f_snapshots(self,q,atom,stochastic=False):
        '''
        Get t_fineness form factors, distributed as evenly between plotter_obj.start_t and plotter_obj.end_t as possible.
        t_fineness = number of snapshots
        f is modified from the typical definition is multiplied by the scalar sqrt(I(t)) to deal with gaussian case.
        
        Returns array with element f_sqrtI(t_i, q_j) given by f_sqrtI_{i,j}
        '''
        def snapshot(idx):
            # We pass in indices from 0 to fineness-1, transform to time:
            t = self.start_t + idx/self.t_fineness*(self.end_t-self.start_t)
            if stochastic:
                val, time_step = self.get_stochastic_form_factor(q,atom,t)
            else:
                val,time_step = self.get_average_form_factor(q,[atom],t)
            idx = np.searchsorted(self.timeData,t)  # SAME AS IN get_average_form_factor()
            val *= np.sqrt(self.intensityData[idx])
            return val, time_step
        form_factors_sqrt_I,time_steps = np.fromfunction(snapshot,(self.t_fineness,))   # (Need to double check working as expected - not using np.vectorise)
        if len(time_steps) != len(np.unique(time_steps)):
            print("times used:", time_steps)
            raise Exception("Error, used same times! Choose a different fineness or a larger range.")
        return form_factors_sqrt_I, time_steps
    

    #############Important bit#######################
    def random_state_snapshots(self,atom,seed=None):
        '''
        Get t_fineness form factors, distributed as evenly between plotter_obj.start_t and plotter_obj.end_t as possible.
        t_fineness = number of snapshots
        f is modified from the typical definition is multiplied by the scalar sqrt(I(t)) to deal with gaussian case.
        '''
        occ_dict = {} # dictionary that converts state index to subshell occupation list
        try:
            states = self.statedict[atom]
        except:
            raise Exception("species "+str(atom)+" data is is missing - data is only present for "+''.join([k + " " for k in self.statedict.keys()]))
        for i in range(len(states)):
            orboccs = parse_elecs_from_latex(states[i])
            occ_list = [-99]*10
            for orb, occ in orboccs.items():
                l = int(orb[0]) - 1
                if occ_list[l] == -99:
                    occ_list[l] = 0
                occ_list[l] += occ
            occ_dict[i] = occ_list[:len(occ_list)-occ_list.count(-99)] 

        time_steps = self.get_times_used()
        orb_occs,time_steps = self.get_random_states(time_steps,atom,seed)

        return orb_occs, time_steps, occ_dict
        
    def random_states_to_f_snapshots(self,time_steps,orb_occ_arr, q, atom,orb_occ_dict):
        idx = np.searchsorted(self.timeData,time_steps)  # SAME AS IN get_average_form_factor()
        if type(q) is np.ndarray:
            if len(q.shape) == 1:
                form_factors_sqrt_I = self.ff_from_state(time_steps,orb_occ_arr,q,atom,orb_occ_dict) * np.sqrt(self.intensityData[idx][...,None])   # (Need to double check working as expected - not using np.vectorise)
            elif len(q.shape) == 2:
                form_factors_sqrt_I = self.ff_from_state(time_steps,orb_occ_arr,q,atom,orb_occ_dict) * np.sqrt(self.intensityData[idx][...,None,None])   # (Need to double check working as expected - not using np.vectorise)
            else:
                raise Exception("unexpected q shape",q.shape)
        else:
            form_factors_sqrt_I = self.ff_from_state(time_steps,orb_occ_arr,q,atom,orb_occ_dict) * np.sqrt(self.intensityData[idx][...])   # (Need to double check working as expected - not using np.vectorise)
        return form_factors_sqrt_I, time_steps        
    
    def ff_from_state(self,time,orb_occ_arr,k,atom,orb_occ_dict):
        shielding = SlaterShielding(self.atomic_numbers[atom])    
        return shielding.get_ff(orb_occ_arr,k,orb_occ_dict)       

    def get_random_states(self,time_steps,atom,seed):  
        '''
        returns single list of random states corresponding to single atom throughout times.   
        '''
        rng = np.random.default_rng(seed)
        idx = np.searchsorted(self.timeData,time_steps)      
        if not np.array_equal(time_steps,self.timeData[idx]):
            raise Exception ("inconsistent_times")   
        num_times = len(time_steps)
        
        a = atom
        states = self.statedict[a]   
        atomic_density = np.sum(self.boundData[a][0, :])

        # Sample configurations from the discrete set
        orboccs = np.empty(shape=num_times)
        roll = rng.random(num_times)  
        for j in range(num_times):  # array compatibility
            cumulative_chance = 0 
            for i in range(len(states)):
                cumulative_chance += max(0,self.boundData[a][idx[j], i]/atomic_density) #TODO figure out why ground state goes slightly neg (~-1/300 of atomic density).
                if cumulative_chance >= roll[j]:
                    orboccs[j] = i 
                    break
                if i == len(states) - 1:
                    print("WARNING, cumulative chance topped below 1 (likely insignificant if only happened a few times). Atomic density:", atomic_density, "state densities:",self.boundData[a][idx[j], :])
                    orboccs[j] = i      
        return orboccs,time_steps

    

    def get_stochastic_form_factor(self,k,atom,time):
        '''
        Returns a stochastic form factor for each momentum transfer 'k' at 'time' [fs], and time used.
        Returns array with element f(t_i, q_j) given by f_{i,j}
        '''        
        orb_occs, time_steps = self.get_random_states(self.get_times_used(),atom)

        ff = self.ff_from_state(time,orb_occs,k,atom) 

        return ff,time_steps      

    def f_undamaged(self,q,atom,ground_state):
        '''        
        Compared to f_average, we now include the stochastic contribution by picking the atomic state from the distribution.
        
        Since we aren't using an average, we have to individually calculate the overall intensity (i.e. contributed to by all atoms) 
        at each point in time.
        I = int{F.F*}dt     (for const intensity) 
        where f(q)_iT(q)_i is the contribution to F by an atom, and f is the time-integrated component (atomic form factor), T(q) is the spatial component.
        
        F(t) = SUM(f(q,t)_i T(q)_i). 
        

        Returns row matrix, with element f(t_i,q_j) = f_{i,j}
        '''
        # return undamaged form factor, multiplied by sqrt of average pulse intensity
        I_avg, time_steps = self.I_avg()
        shielding = SlaterShielding(self.atomic_numbers[atom])             
        ff = shielding.get_ff("dummy",q,{"dummy":ground_state})
        f_sqrt_I = ff[...,None] * np.array(np.sqrt(I_avg)*np.ones(len(time_steps)))
        return np.moveaxis(f_sqrt_I,len(f_sqrt_I.shape)-1,0), time_steps                     
    # why did I do thissss
    def get_times_used(self):
        def snapshot(idx):
            # We pass in indices from 0 to fineness-1, transform to time:
            idx = np.searchsorted(self.timeData, self.start_t + idx/self.t_fineness*(self.end_t-self.start_t))            
            try:
                return self.timeData[idx]
            except:
                raise Exception("Start and end times provided seem to be outside the range of the output data.")
        time_steps = np.fromfunction(snapshot,(self.t_fineness,))
        if len(time_steps) != len(np.unique(time_steps)):
            print("times used:", time_steps)
            raise Exception("Error, used same times! Choose a different fineness or a larger range.")
        return time_steps      
    def I_avg(self): # average intensity for pulse 
        def snapshot(idx):
            # We pass in indices from 0 to fineness-1, transform to time:
            t = self.start_t + idx/self.t_fineness*(self.end_t-self.start_t)
            idx = np.searchsorted(self.timeData,t)  # SAME AS IN get_average_form_factor()
            time_step = self.timeData[idx]
            val = self.intensityData[idx]
            return val, time_step      
        I,time_steps = np.fromfunction(snapshot,(self.t_fineness,))
        I_avg = np.trapz(I,time_steps)/(time_steps[-1]-time_steps[0])
        return I_avg,time_steps

    ####################################
    def f_average(self,q,atom): 
        ''' 
        If we use an average form factor, then we are making the approximation that for a const intensity I ~= int{F.F*}dt =int{f^2}dt T.T*, 
        This function returns f, modified by the relative intensity so f-> sqrt(I_tot(t)) f
        where F(q)_i = f(q)T(q)_i is the contribution to F by an atom, and f is the time-integrated component, T(q) is the spatial component.
        def f_average(self,q,atom):
            if self.end_t == self.start_t:
                return self.get_average_form_factor(q,[atom],self.end_t)[0]

        
        '''
        if self.end_t == self.start_t:
            return self.get_average_form_factor(q,[atom],self.end_t)[0]

        form_factors_sqrt_I, time_steps = self.f_snapshots(q,atom)
        #Approximate integral with composite trapezoidal rule.
        f_avg = np.trapz(form_factors_sqrt_I,time_steps)/(time_steps[-1]-time_steps[0])
        return f_avg         

    def plot_form_factor(self,num_snapshots = 8, atoms=None,resolution=False,q_min=0,fig_width=5,fig_height=6,show_average=True,percentage_change=False):  # k = 4*np.pi = 2*q. c.f. Sanders plot.
        
        times = np.linspace(self.start_t,self.end_t,num_snapshots)
        if atoms == None:
            atoms = self.atomdict  # All atoms

        #lattice_const = 18.11**(1/3)


        #TODO figure out overlaying.
        q = np.linspace(q_min,self.q_max,self.q_fineness)
        if resolution:
            q = q[q!=0]
        self.fig = plt.figure(figsize=(3,2.5))
        ax = self.fig.add_subplot(111)

        cmap=plt.get_cmap('plasma')
        min_time = np.inf
        max_time = -np.inf 
        for time in times:
            min_time = min(min_time,time)
            max_time = max(max_time,time) 

        f_avg =[]
        print(times)
        ang_per_bohr = 0.529177
        angstrom_mom = [k/ang_per_bohr for k in q]
        X = np.array(angstrom_mom)
        if resolution:
            X = 2*np.pi/X       
        if percentage_change:
            f_init = np.array([self.get_average_form_factor(x,atoms,time=self.timeData[0])[0] for x in q])
        for time in times:
            f = np.array([self.get_average_form_factor(x,atoms,time=time)[0] for x in q])
            if percentage_change:
                f = -100*(1-f/f_init)
            f_avg.append(f)
            #ax.plot(X, f,label="t = " + str(time)+" fs",color=cmap((time-min_time)/(0.0001+max_time-min_time)))     
            ax.plot(X, f,label="%.1f"%time+" fs",color=cmap((time-min_time)/(0.0001+max_time-min_time)))     
            # Percentage difference from initital state
        if show_average:
            f_avg = np.average(f_avg,axis=0)
            ax.plot(X, f_avg,'k--',label="Average")  

        ax.set_title("")
        ax.set_xlabel("Resolution ($\\AA^{-1}$)")
        ax.set_ylabel("Form factor")
        if percentage_change:
            ax.set_ylabel("Change in $f(q)$ (\\%)")
        #ax.set_xlim(0,self.q_max/ang_per_bohr)   
        ax.set_xlim(np.min(X),np.max(X))
        if resolution:
            ax.invert_xaxis()
        loc = "upper right"
        if percentage_change:
            loc = "lower right"
        ax.legend(loc=loc)
        self.fig.set_size_inches(fig_width,fig_height)           

    
    # Note this combines form factors when supplied multiple species, which is not necessarily interesting.
    # TODO overly complicated..
    def get_average_form_factor(self,k,atoms,time=-7.5,n=None):
        '''
        Returns the form factor(s) at 'time' [fs], and time(s) used.
        use n for average form factor of a specific shell.
        '''        
        idx = np.searchsorted(self.timeData,time)      
        time_step = self.timeData[idx]  # The actual time used

        # Get relevant prevalence of each species for atoms provided.
        tot_density = 0
        for a in atoms:   
            tot_density += self.boundData[a][0, 0]
        # Iterate through each a passed. 
        ff = 0
        for a in atoms:
            states = self.statedict[a]   
            atomic_density = self.boundData[a][0, 0]
            for i in range(len(states)):
                if n != None:
                    if n != i+1: continue
                orboccs = parse_elecs_from_latex(states[i])             
                state_density = self.boundData[a][idx, i]       
                # Get the form factors for each subshell
                occ_list = [-99]*10
                for orb, occ in orboccs.items():
                    l = int(orb[0]) - 1
                    if occ_list[l] == -99:
                        occ_list[l] = 0
                    occ_list[l] += occ
                occ_list = occ_list[:len(occ_list)-occ_list.count(-99)]
                occ_dict = {1:occ_list}
                shielding = SlaterShielding(self.atomic_numbers[a])             
                ff += state_density/tot_density * shielding.get_ff(1,k,occ_dict)
                #print("FORM FACTOR IDEAL",ff)
        return ff,time_step    
    
    # Used by scatter code. Returns ground state's shell occupancies. (Assumes only s and p orbitals).
    def get_ground_state_shells(self,atom):
        orboccs = parse_elecs_from_latex(self.statedict[atom][0])
        occ_list = [-99]*10
        for orb, occ in orboccs.items():
            l = int(orb[0]) - 1
            if occ_list[l] == -99:
                occ_list[l] = 0
            occ_list[l] += occ     
        return occ_list[:len(occ_list)-occ_list.count(-99)]


    def print_bound_slice(self,time=-7.5):
        idx = np.searchsorted(self.timeData,time)
        for a in self.atomdict:
            states = self.statedict[a]

            print("Atom:",a)

            atomic_number = 0

            # Get atomic number and density by considering states at time 0.
            
            for val in parse_elecs_from_latex(states[0]).values():
                atomic_number += val
            atomic_density = self.boundData[a][0, 0]
            print("atomic number, density",atomic_number,atomic_density)
            
            average_occupancy = dict.fromkeys(parse_elecs_from_latex(states[0]), 0)
            for i in range(len(states)):
                orboccs = parse_elecs_from_latex(states[i])              
                state_density = self.boundData[a][idx, i]
                print("Occs: ", orboccs ,", time: ",self.timeData[idx], ", density: ",state_density)
                # Get the average number of orbitals in each state
                for orb, occ in orboccs.items(): 
                    average_occupancy[orb] += occ*(state_density/atomic_density)

            print("Average occupancy:",average_occupancy)
            print("Atomic density:", atomic_density)

    def aggregate_charges(self,charge_difference=False):
        # populates self.chargeData based on contents of self.boundData
        for a in self.atomdict:
            states = self.statedict[a]
            if len(states) != self.boundData[a].shape[1]:
                msg = 'states parsed from file header disagrees with width of data width'
                msg += ' (got %d, expected %d)' % (len(states), self.boundData[a].shape[1])
                raise RuntimeError(msg)

            initial_charge = 0
            if charge_difference:
                initial_charge = ATOMNO[a] - sum(parse_elecs_from_latex(states[0]).values()) 
            self.chargeData[a] = np.zeros((self.boundData[a].shape[0], ATOMNO[a]+1))
            for i in range(len(states)):
                orboccs = parse_elecs_from_latex(states[i])
                charge = ATOMNO[a] - sum(orboccs.values()) - initial_charge 
                self.chargeData[a][:, charge] += self.boundData[a][:, i]




    def update_outputs(self):
        
        num_sample_lines = 10
        if self.sample_end_points:
            with open(self.intFile,'rb') as f:
                lines = f.readlines()
                raw = np.genfromtxt(lines[-num_sample_lines:], comments='#', dtype=np.float64)
        else:
            raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        self.timeData = raw[:, 0]
        self.intensityData = raw[:,1]
        self.energyKnot = np.array(self.get_free_energy_spec(), dtype=np.float64)
        if self.sample_end_points:
            with open(self.freeFile,'rb') as f:
                lines = f.readlines()
                raw = np.genfromtxt(lines[-num_sample_lines:], comments='#', dtype=np.float64)
        else:
            tmp = []
            try:
                for elemContinuum in self.freeFiles:
                        tmp.append(np.genfromtxt(elemContinuum, comments='#', dtype=np.float64))
                raw = tmp[0]
            except:
                # In case of not tracking split continuums.
                tmp = [np.genfromtxt(self.freeFile, comments='#', dtype=np.float64)]
            if len(tmp) > 1:
                for r in tmp[1:]:
                    raw += r


        self.freeData = raw[:,1:]
        photo_data_present = False
        for a in self.atomdict:
            if self.sample_end_points:
                with open(self.atomdict[a]['outfile'],'rb') as f:
                    lines = f.readlines()
                    raw = np.genfromtxt(lines[-num_sample_lines:], comments='#', dtype=np.float64)
            else:
                raw = np.genfromtxt(self.atomdict[a]['outfile'], comments='#', dtype=np.float64)

            self.boundData[a] = raw[:, 1:]
            self.statedict[a] = self.get_bound_config_spec(a)
            # rates
            try:
                if self.sample_end_points:
                    with open(self.atomdict[a]['photofile'],'rb') as f:
                        lines = f.readlines()
                        raw = np.genfromtxt(lines[-num_sample_lines:], comments='#', dtype=np.float64)
                else:
                    raw = np.genfromtxt(self.atomdict[a]['photofile'], comments='#', dtype=np.float64)
                
                self.photoData[a] = raw[:, 1] 
                photo_data_present = True
            except:
                print("Warning: Missing '" + self.atomdict[a]['photofile'] + "'.")
        # Truncate data to time specified.
        if self.end_t_plotting is not None: 
            last_idx = np.searchsorted(self.timeData,self.end_t_plotting)
            self.timeData = self.timeData[0:last_idx]
            self.intensityData = self.intensityData[0:last_idx]
            self.freeData = self.freeData[0:last_idx]
            for a in self.atomdict:
                self.boundData[a] = self.boundData[a][0:last_idx] 
                if photo_data_present:
                    self.photoData[a] = self.photoData[a][0:last_idx]
                

        self.atomic_numbers = self.get_atomic_numbers()
        self.grid_update_time_Data = []
        self.grid_point_Data = []
        with open(self.gridFile) as file:
            for n, row in enumerate(csv.reader(file)):
                if n < 2:
                    continue
                self.grid_update_time_Data.append(float(row[0].split()[0])) 
                self.grid_point_Data.append([float(elem) for elem in row[0].split()[1:]])
        self.grid_update_time_Data = np.array(self.grid_update_time_Data,dtype=np.float64)
        self.grid_point_Data = np.array(self.grid_point_Data,dtype=object)

    def go(self):
        if not self.check_current():
            self.rerun_ac4dc()
            self.update_outputs()

    # makes a blank plot showing the intensity curve
    def setup_intensity_plot(self,ax,show_pulse_profile=True,col='black'):
        ax2 = ax.twinx()
        if show_pulse_profile:
            ax2.plot(self.timeData, self.intensityData, lw = 2, c = col, ls = ':', alpha = 0.7)
        # ax2.set_ylabel('Pulse Intensity (photons cm$^{-2}$ s$^{-1}$)')
        ax.set_xlabel("Time (fs)")
        ax.tick_params(direction='in')
        ax2.axes.get_yaxis().set_visible(False)
        # ax2.tick_params(None)
        ax.get_xaxis().get_major_formatter().labelOnlyBase = False
        ax2.set_ylim([0,ax2.get_ylim()[1]-ax2.get_ylim()[0]])
        #self.fig.subplots_adjust(left=0.11, right=0.81, top=0.93, bottom=0.1)   
        return (ax, ax2)
    
    def setup_axes(self,num_subplots):
        self.num_plotted = 0 # number of subplots plotted so far.
        width, height = 6, 3.3333  # 4.5,2.5 ~ abdallah

        if num_subplots >= 3:
            self.fig, self.axs = plt.subplots(int(0.999+(num_subplots**0.5)),int(0.999+(num_subplots**0.5)),figsize=(width*int((1+num_subplots)/2),height*int((1+num_subplots)/2)))
        else:
            self.fig, self.axs = plt.subplots(num_subplots,figsize=(width,height*num_subplots))
            if type(self.axs) is not np.ndarray:
                self.axs = np.array([self.axs])
        self.num_subplots = num_subplots

    def get_next_ax(self):
        ax = self.axs.flat[self.num_plotted]
        self.num_plotted+=1
        return ax
    def delete_remaining_axes(self):
        if self.num_subplots == 1:
            return
        num_plot_spaces = self.axs.shape[0]
        if self.num_subplots > 2:
            num_plot_spaces*=self.axs.shape[1]
        while self.num_plotted < num_plot_spaces:
            self.axs.flat[self.num_plotted].remove()
            self.num_plotted+=1

        
    def plot_atom_total(self, a):
        ax, ax2 = self.setup_intensity_plot(self.get_next_ax())
        tot = np.sum(self.boundData[a], axis=1)
        ax.plot(self.timeData, tot)
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        self.fig.figlegend(loc = (0.11, 0.43))
        #plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)
        
    def plot_atom_raw(self, a):
        ax, ax2 = self.setup_intensity_plot(self.get_next_ax())
        for i in range(self.boundData[a].shape[1]):
            ax.plot(self.timeData, self.boundData[a][:,i], label = self.statedict[a][i])
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        #self.fig.subplots_adjust(left=0.2, right=0.92, top=0.93, bottom=0.1)

    def plot_charges(self, a, ion_fract = True, rseed=404,plot_legend=True,show_pulse_profile=True,xlim=[None,None],ylim=[0,1],**kwargs):
        ax, ax2 = self.setup_intensity_plot(self.get_next_ax(),show_pulse_profile=show_pulse_profile)
        self.aggregate_charges()
        #print_idx = np.searchsorted(self.timeData,-7.5)
        ax.set_prop_cycle(rcsetup.cycler('color', get_colors(self.chargeData[a].shape[1],rseed)))
        max_at_zero = np.max(self.chargeData[a][0,:]) 
        norm = 1
        if ion_fract:
            norm = 1/max_at_zero            
        num_traces = 0 
        for i in range(self.chargeData[a].shape[1]):
            mask = self.chargeData[a][:,i] > max_at_zero*2
            mask |= self.chargeData[a][:,i] < -max_at_zero*2
            Y = np.ma.masked_where(mask, self.chargeData[a][:,i])     
            if np.sum(Y) == 0:
                continue                 
            ax.plot(self.timeData, Y*norm, label = "%d+" % i,**kwargs)
            num_traces+=1 
            #print("Charge: ", i ,", time: ",self.timeData[print_idx], ", density: ",self.chargeData[a][print_idx,i])
        # ax.set_title("Charge state dynamics")
        ax.set_ylabel(r"Density (\AA$^{-3}$)")
        if ion_fract:
            ax.set_ylabel(r"Ion Fractions")
        old_ytop = ax.get_ylim()[1]
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        ax2.set_ylim([0,ax2.get_ylim()[1]*ax.get_ylim()[1]/old_ytop])

        num_cols = 1+round(num_traces/30-num_traces%30/30)
        #self.fig.subplots_adjust(left=0.11, right=0.81, top=0.93, bottom=0.1)
        if plot_legend:
            ax.legend(loc='upper left',bbox_to_anchor=(1, 1),fontsize=4,ncol=num_cols)
    # TODO Remove this function (not the royle one though it's different enough)
    def plot_charges_leonov_style(self, a, ion_fract = True, rseed=404,plot_legend=True,show_pulse_profile=True,xlim=[None,None],ylim=[0,1],**kwargs):
        ax, ax2 = self.setup_intensity_plot(self.get_next_ax(),show_pulse_profile=show_pulse_profile)
        self.aggregate_charges()
        #print_idx = np.searchsorted(self.timeData,-7.5)
        max_at_zero = np.max(self.chargeData[a][0,:]) 
        norm = 1
        if ion_fract:
            norm = 1/max_at_zero            
        num_traces = 0 
        colours = ["#000000","#CF5346","#91E95C","#3E5CE7","#BDF9E1","#9D83CC","#EEE483","#858F2F","#3C44A1","#7541D9","#A77375","#59B535","#64A68E","#5465ED","#F0A078"]
        for i in range(self.chargeData[a].shape[1]):
            mask = self.chargeData[a][:,i] > max_at_zero*2
            mask |= self.chargeData[a][:,i] < -max_at_zero*2
            Y = np.ma.masked_where(mask, self.chargeData[a][:,i])     
            if np.sum(Y) == 0:
                continue                 
            ax.plot(self.timeData, Y*norm, label = a+r"$^{%d+}$" % i,color=colours[i],**kwargs)
            num_traces+=1 
            #print("Charge: ", i ,", time: ",self.timeData[print_idx], ", density: ",self.chargeData[a][print_idx,i])
        # ax.set_title("Charge state dynamics")
        ax.set_ylabel(r"Density (\AA$^{-3}$)")
        if ion_fract:
            ax.set_ylabel("Ion Fractions")
        old_ytop = ax.get_ylim()[1]
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        ax2.set_ylim([0,ax2.get_ylim()[1]*ax.get_ylim()[1]/old_ytop])
        num_cols = 1+round(num_traces/30-num_traces%30/30)
        #self.fig.subplots_adjust(left=0.11, right=0.81, top=0.93, bottom=0.1)
        if plot_legend:
            ax.legend(loc='upper left',bbox_to_anchor=(1, 1),fontsize=4,ncol=num_cols)  

    def plot_charges_royle_style(self, a, ion_fract = True,max_charge=11):
        ax, ax2 = self.setup_intensity_plot(self.get_next_ax(),show_pulse_profile=False)
        self.aggregate_charges()
        max_at_zero = np.max(self.chargeData[a][0,:]) 
        norm = 1
        if ion_fract:
            norm = 1/max_at_zero            
        num_traces = 0
        under_three_charge = None
        colours = ["#000000","#942192","#FF2600","#FF7C00","#F6D800","#008F00","#4180FF","#0433FF","#021CA1"]
        linestyle = ["dashed"] + ["solid"]*8
        for i in range(self.chargeData[a].shape[1]):
            if i >max_charge:
                break
            mask = self.chargeData[a][:,i] > max_at_zero*2
            mask |= self.chargeData[a][:,i] < -max_at_zero*2
            Y = np.ma.masked_where(mask, self.chargeData[a][:,i])     
            if np.sum(Y) == 0:
                continue                 
            if i == 0:
                under_three_charge = Y
                continue
            elif i <=3:
                under_three_charge+=Y
                if i != 3:
                    continue
                Y = under_three_charge

            ax.plot(self.timeData, Y*norm, label = "%d+" % i,color=colours[i-3],linestyle=linestyle[i-3])
            ax.set_ylim([0,0.8])
            ax.set_xlim([-100,100])
            num_traces+=1 
            #print("Charge: ", i ,", time: ",self.timeData[print_idx], ", density: ",self.chargeData[a][print_idx,i])
        # ax.set_title("Charge state dynamics")
        ax.set_ylabel(r"Density (\AA$^{-3}$)")
        if ion_fract:
            ax.set_ylabel(r"Ion Fraction")

        num_cols = 1+round(num_traces/30-num_traces%30/30)
        ax.legend(loc='upper left',bbox_to_anchor=(1, 1),fontsize=4,ncol=num_cols)

    def plot_charges_bar(self, a, ion_fract = True, rseed=404,plot_legend=True,show_pulse_profile=True,xlim=[None,None],ylim=[0,1],**kwargs):
        ax = self.get_next_ax()
        self.aggregate_charges()
        #print_idx = np.searchsorted(self.timeData,-7.5)
        ax.set_prop_cycle(rcsetup.cycler('color', get_colors(self.chargeData[a].shape[1],rseed)))
        max_at_zero = np.max(self.chargeData[a][0,:]) 
        norm = 1
        if ion_fract:
            norm = 1/max_at_zero            
        num_traces = 0             

        Y = []
        Z = []
        for i in range(self.chargeData[a].shape[1]):
            Y.append(i)   
            Z.append(norm*self.chargeData[a][:,i])  
            num_traces+=1 

        ax.set_facecolor('black')
        cm = ax.pcolormesh(self.timeData, Y, Z, shading='inferno',cmap="inferno",rasterized=True)
        cbar = self.fig.colorbar(cm,ax=ax,label="State density")

        ax.set_xlabel("Time (fs)")            

        ax.set_ylabel(r"Density (\AA$^{-3}$)")
        if ion_fract:
            ax.set_ylabel(a+r" charge")
        old_ytop = ax.get_ylim()[1]
        ax.set_ylim([-0.5,len(Y)-0.5])

    def plot_orbitals_bar(self, atoms = None, rseed=404,plot_legend=True,show_pulse_profile=True,xlim=[None,None],orbitals=None,normalise=False,atoms_excluded = None,label_all_yaxes=False,**kwargs):
        if atoms is None: 
            atoms = self.atomdict
            if atoms_excluded is not None:
                        for a in atoms_excluded:
                            atoms.pop(a)          
        else:
            assert atoms_excluded is None
        
        yaxis_plotted = False
        for a in atoms:
            # if a not in self.atomdict:
            #     continue
            if show_pulse_profile:  
                ax, ax2 = self.setup_intensity_plot(self.get_next_ax(),col="white")
            else:
                ax = self.get_next_ax()
            self.aggregate_charges()
            #print_idx = np.searchsorted(self.timeData,-7.5)
            ax.set_prop_cycle(rcsetup.cycler('color', get_colors(self.chargeData[a].shape[1],rseed)))                        

            Z,labels = self.get_orbital_data(a,orbitals)
            Y = range(len(Z))
            if normalise:
                Z = np.array(Z)[:]/np.array(Z)[:,0][:,None]
            ax.set_facecolor('black')
            vmax = None
            cm = ax.pcolormesh(self.timeData, Y, Z,cmap="Spectral",rasterized=True,vmin=0)
            z_label = "Avg. "+a.split("_")[0]+" orbital occupancy"
            if normalise:
                z_label = a.split("_")[0] + " orbital density"
            ax.set_yticks(ticks=Y,labels=labels)
            
            ax.set_xlabel("Time (fs)")            
            if label_all_yaxes or yaxis_plotted is False:
                ax.set_ylabel(r"Orbital")
                yaxis_plotted = True
            old_ytop = ax.get_ylim()[1]
            ax.set_ylim([-0.5,len(Y)-0.5])  
            #ax.set_xlim([None,0])
            #ax.set_xticks(np.arange(-15,0.1,5))
            if show_pulse_profile:   
                ax2.set_ylim([0,ax2.get_ylim()[1]*ax.get_ylim()[1]/old_ytop])
        if len(atoms) == 2:
            cbar_ax = self.fig.add_axes([0.9, 0.492, 0.02, 0.433])
            cbar = self.fig.colorbar(cm, cax=cbar_ax)
            plt.subplots_adjust(left=0.105, right=0.875, top=0.93, bottom=0)
        elif len(atoms) == 1:
            cbar_ax = self.fig.add_axes([0.9, 0.492, 0.02, 0.433])
            cbar = self.fig.colorbar(cm, cax=cbar_ax)
            plt.subplots_adjust(left=0.105, right=0.875, top=0.93, bottom=0)

        else:
            cbar = self.fig.colorbar(cm,ax=ax,label=z_label,format="%.1f",)        
        

    def get_orbital_data(self,a,orbitals):
        self.aggregate_charges(False)
        max_at_zero = np.max(self.chargeData[a][0,:]) 
        norm = 1/max_at_zero        

        colour = None
        states = self.statedict[a]
        orb_dict = parse_elecs_from_latex(states[0])
        if orbitals is not None:
            orb_dict = {k: v for k, v in orb_dict.items() if k in orbitals}

        orbital_idx = []
        orb_density = []
        labels = []
        for j, subshell in enumerate(orb_dict.keys()):
            orbital_idx.append(j)
            labels.append(subshell)
            tot = np.zeros_like(self.boundData[a][:,0])
            for i in range(len(states)):
                orboccs = parse_elecs_from_latex(states[i])
                tot += orboccs[subshell]*self.boundData[a][:,i]*norm
            orb_density.append(tot)      
        return orb_density, labels             
    
    def plot_photoionisation(self,atoms,show_pulse_profile=True):
        if atoms is None: 
            atoms = self.atomdict        
        for a in atoms:
            if show_pulse_profile:  
                ax, ax2 = self.setup_intensity_plot(self.get_next_ax())
            else:
                ax = self.get_next_ax()        
            ax.scatter(self.timeData,self.photoData[a])
            ax.ticklabel_format(style='sci',scilimits=(0,0))
            ax.set_ylabel(a+" species-wide cumulative photoionisations")
            

    def plot_subshell(self, a, subshell='1s',rseed=404):
        if not hasattr(self, 'ax_subshell'):
            self.ax_subshell, _ax2 = self.setup_axes()

        states = self.statedict[a]
        tot = np.zeros_like(self.boundData[a][:,0])
        for i in range(len(states)):
            orboccs = parse_elecs_from_latex(states[i])
            tot += orboccs[subshell]*self.boundData[a][:,i]

        self.ax_subshell.plot(self.timeData, tot, label=subshell)
        # ax.set_title("Charge state dynamics")
        self.ax_subshell.set_ylabel(r"Electron density (\AA$^{-3}$)")

        self.fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)


    def plot_tot_charge(self, every=1,densities = False,colours=None,atoms=None,plot_legend=True,charge_difference=True,xlim=[None,None],ylim=[None,None],plot_derivative=False,legend_loc='upper left',**kwargs):
        '''
        plot_derivative (bool), if True, plots average ionisation rate instead of average charge. 
        '''
        ax, ax2 = self.setup_intensity_plot(self.get_next_ax())
        #self.fig.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.17)
        T = self.timeData[::every]
        T_start = 0
        if ylim[0] is None:
            ylim[0] = 0
        if xlim[0] is not None:
            T_start = np.searchsorted(T,xlim[0])
        T_end = len(T)
        if xlim[1] is not None:
            T_end = np.searchsorted(T,xlim[1])
        T = T[T_start:T_end]
        
        self.aggregate_charges(charge_difference)
        self.Q = np.zeros(T.shape[0]) # total charge
        colour = None
        if atoms is None:
            atoms = self.atomdict
        for j,a in enumerate(atoms):
            if atoms is not None and a not in atoms:
                continue
            if colours != None:
                colour = colours[j]
            kwargs["label"] = a
            kwargs["color"] = colour
            # Plot trace for the atom
            atomic_charge = np.zeros(T.shape[0])
            for i in range(self.chargeData[a].shape[1]):
                atomic_charge += self.chargeData[a][::every,i][T_start:T_end]*i
            if not densities:
                atomic_charge /= np.sum(self.chargeData[a][0])
            if not plot_derivative:
                ax.plot(T,atomic_charge,**kwargs)
            else:
                smooth = False
                dydx = np.gradient(atomic_charge,T)
                if not smooth:
                    ax.plot(T, dydx,**kwargs)
                else:
                    w = np.zeros(len(T))
                    w+=1
                    from scipy.ndimage import uniform_filter1d
                    spikiness = np.abs(uniform_filter1d(splev(T,splrep(T,atomic_charge,k=3,s=0),der=3),size=10))
                    w/=np.sqrt(spikiness)
                    w /= np.median(w)
                    w[0] = 100
                    smoothed = splev(T,splrep(T,atomic_charge,w=w,k=3,s=0.015),der=1)
                    smoothed *= np.max(dydx)/np.max(smoothed)
                    
                    #smoothed = uniform_filter1d(atomic_charge,size=200,mode="reflect")
                    #smoothed = splev(T,splrep(T,smoothed,k=3,s=0),der=1)
                    #smoothed *= np.max(dydx)/np.max(smoothed)

                    YY = smoothed

                    ax.plot(T, YY,**kwargs)
                
                
                #dt = np.append(T, T[-1]*2 - T[-2])   
                #dt = dt[1:] - dt[:-1]
                #dydx = savgol_filter(atomic_charge, window_length=11, polyorder=4, deriv=1)*10
                #ax.plot(T,dydx,**kwargs)
                

            self.Q += atomic_charge

        # Free data
        if densities:
            de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])   
            de = de [1:] - de[:-1]
            tot_free_Q =-1*np.dot(self.freeData, de)
            ax.plot(T, tot_free_Q[::every], label = 'Free')
            self.Q += tot_free_Q[::every]
            # ax.set_title("Charge Conservation")
            if not plot_derivative:
                ax.plot(T, self.Q, label='total')
                ax.set_ylabel("Charge density ($e$ \AA$^{-3}$)")
            else:
                ax.plot(T, np.gradient(self.Q,T), label='total')
                ax.set_ylabel("dQ/dt ($e$ \AA$^{-3} \cdot fs^{-1}$)")
        else:
            ax.set_ylabel("Average charge")
            if plot_derivative:
                ax.set_ylabel("Average ionisation rate ($e \cdot fs^{-1}$)")
            elif charge_difference:
                #ax.set_ylabel("Avg. charge difference")  # Sometimes we start with charged states.
                ax.set_ylabel("Charge gain")  # Sometimes we start with charged states.
        ax.set_xlim([T[0],T[-1]])
        ax.set_ylim(ylim)
        #ax.yaxis.set_ticks([0,0.2,0.4,0.6,0.8,1])
        #ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        # ax.set_xlim(None,-18.6)
        # ax.set_ylim(0,5)
        if plot_legend:
            ax.legend(loc = legend_loc)
        return ax
    
    def plot_orbitals_charge(self, every=1,densities = False,cmap=None,atom=None,orbitals = None,plot_legend=True,xlim=[None,None],ylim=[None,None],plot_derivative=False,legend_loc='lower left',custom_legend=None,**kwargs):
        '''
        plot_derivative (bool), if True, plots average ionisation rate instead of average charge. 
        '''
        ax, ax2 = self.setup_intensity_plot(self.get_next_ax())
        #self.fig.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.17)
        T = self.timeData[::every]
        T_start = 0
        if xlim[0] is not None:
            T_start = np.searchsorted(T,xlim[0])
        T_end = len(T)
        if xlim[1] is not None:
            T_end = np.searchsorted(T,xlim[1])
        T = self.timeData[T_start:T_end]        
        
        avg_occupancies,labels = self.get_orbital_data(atom,orbitals)

        if cmap is None:
            cmap = "tab20"

        colours = [plt.get_cmap(cmap)((i)/len(avg_occupancies)) for i in range(len(avg_occupancies)) ]

        for Y, label, col in zip(avg_occupancies,labels,colours):
            Y = Y[::every][T_start:T_end]
            if not plot_derivative:
                ax.plot(T,Y,label=label,color=col,**kwargs)
            else:
                smooth = False
                dydx = np.gradient(Y,T)
                if not smooth:
                    ax.plot(T, dydx,label=label,**kwargs)
                else:
                    w = np.zeros(len(T))
                    w+=1
                    from scipy.ndimage import uniform_filter1d
                    spikiness = np.abs(uniform_filter1d(splev(T,splrep(T,Y,k=3,s=0),der=3),size=10))
                    w/=np.sqrt(spikiness)
                    w /= np.median(w)
                    w[0] = 100
                    smoothed = splev(T,splrep(T,Y,w=w,k=3,s=0.015),der=1)
                    smoothed *= np.max(dydx)/np.max(smoothed)
                    
                    #smoothed = uniform_filter1d(atomic_charge,size=200,mode="reflect")
                    #smoothed = splev(T,splrep(T,smoothed,k=3,s=0),der=1)
                    #smoothed *= np.max(dydx)/np.max(smoothed)

                    YY = smoothed

                    ax.plot(T, YY,label=label**kwargs)
                
                
                #dt = np.append(T, T[-1]*2 - T[-2])   
                #dt = dt[1:] - dt[:-1]
                #dydx = savgol_filter(atomic_charge, window_length=11, polyorder=4, deriv=1)*10
                #ax.plot(T,dydx,**kwargs)
                


        # Free data
        if densities:
            assert False,"Densities not implemented with orbitals charge plot"
        #     de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])   
        #     de = de [1:] - de[:-1]
        #     tot_free_Q =-1*np.dot(self.freeData, de)
        #     ax.plot(T, tot_free_Q[::every], label = 'Free')
        #     self.Q += tot_free_Q[::every]
        #     # ax.set_title("Charge Conservation")
        #     if not plot_derivative:
        #         ax.plot(T, self.Q, label='total')
        #         ax.set_ylabel("Charge density ($e$ \AA$^{-3}$)")
        #     else:
        #         ax.plot(T, np.gradient(self.Q,T), label='total')
        #         ax.set_ylabel("dQ/dt ($e$ \AA$^{-3} \cdot fs^{-1}$)")
        # else:
        #     ax.set_ylabel("Average charge")
        #     if plot_derivative:
        #         ax.set_ylabel("Average ionisation rate ($e \cdot fs^{-1}$)")
        #     elif charge_difference:
        #         ax.set_ylabel("Avg. charge difference")  # Sometimes we start with charged states.
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel("Average occupancy")
        # ax.set_xlim(None,-18.6)
        # ax.set_ylim(0,5)
        if plot_legend:
            #ax.legend(loc = legend_loc)
            if custom_legend is None: 
                ax.legend(bbox_to_anchor=(1.02, 1),loc='upper left', ncol=1,handlelength=1)  # Top right legend.
            else:
                handles,_ = ax.get_legend_handles_labels()
                handles = list(handles)
                assert len(custom_legend)==len(handles)
                ax.legend(handles,custom_legend,bbox_to_anchor=(1.02, 1),loc='upper left', ncol=1,handlelength=1)
        return ax    
    
 
    
    def get_total_charge(self,every=1,densities=False,atoms=None,charge_difference=False):
        T = self.timeData[::every]
        self.aggregate_charges(charge_difference)
        self.Q = np.zeros(T.shape[0]) # total charge
        for j,a in enumerate(self.atomdict):
            if atoms is not None and a not in atoms:
                continue
            atomic_charge = np.zeros(T.shape[0])
            for i in range(self.chargeData[a].shape[1]):
                atomic_charge += self.chargeData[a][::every,i]*i
            if not densities:
                atomic_charge /= np.sum(self.chargeData[a][0])
            self.Q += atomic_charge
        # Free data
        if densities:
            de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])   
            de = de [1:] - de[:-1]
            tot_free_Q =-1*np.dot(self.freeData, de)
            self.Q += tot_free_Q[::every]
        return self.Q
            
        

    def plot_all_charges(self, ion_fract = True, rseed=404,plot_legend=True,show_pulse_profile=True,xlim=[None,None],ylim=[None,None],**kwargs):
        for a in self.atomdict:
            self.plot_charges(a, ion_fract, rseed,plot_legend,show_pulse_profile=show_pulse_profile,xlim=xlim,ylim=ylim,**kwargs)
        
    def plot_free(self, N=100, log=True, cmin = 1e-9, cmax=None, every = None,mask_below_min=True,cmap='magma',ylim=[None,None],ymax=np.Infinity,leonov_style = False,keV=False,ylog=False):
        ax = self.get_next_ax()
        #self.fig_free.subplots_adjust(left=0.12, top=0.96, bottom=0.16,right=0.95)

        if every is None:
            Z = self.freeData.T
            T = self.timeData
        else:
            Z = self.freeData.T [:, ::every]
            T = self.timeData [::every]
        
        if mask_below_min:
            min_col = 0 
            if leonov_style:  
                min_col = 20
                Z = np.ma.masked_where(Z <= 0,Z)
                Z = np.ma.masked_where(Z*4200*np.linalg.norm(np.log10(Z[Z.mask==False]))<min_col,Z) # 2550  4200 no time for math time for eye. Desperate times...
            else:
                Z = np.ma.masked_where(Z <= cmin, Z)
            cmap = plt.get_cmap(cmap)
            cmap.set_over(cmap(np.inf))
            cmap.set_under(cmap(min_col))    
            cmap.set_extremes(bad=cmap(min_col))        
            if cmax is not None:
                pass
            #     Z = np.ma.masked_where(Z > max, Z)
            else:
                cmax = Z.max()
        if cmin == 0 and log:
            cmin = Z.min()

        # if log:
        #     Z = np.log10(Z)

        ax.set_facecolor('black')
        
        norm = colors.LogNorm(vmin=cmin, vmax=cmax)
        if log:
            scale = 1
            if keV:
                scale = 1e-3
            # Removing erroneous data points hack
            Z = Z[:,(T!=-2.778)&(T!=4.422)]
            self.intensityData = self.intensityData[(T!=-2.778)&(T!=4.422)]
            T = T[(T!=-2.778)&(T!=4.422)]
            cm = ax.pcolormesh(T, self.energyKnot[self.energyKnot<ymax]*scale, Z[self.energyKnot<ymax], shading='gouraud',norm=norm,cmap=cmap,rasterized=True)
            cbar = self.fig.colorbar(cm,ax=ax)
        else:
            if cmap != None:
                shading='gouraud'            
            cm = ax.contourf(T, self.energyKnot, Z, N, cmap=cmap,rasterized=True)
            cbar = self.fig.colorbar(cm,ax=ax)

        if keV:
            ax.set_ylabel("Energy (keV)")
        else:
            ax.set_ylabel("Energy (eV)")
            
        ax.set_xlabel("Time (fs)")
       
        if ylog:
            ax.set_yscale('log')
        if keV:
            ylim = [y*1e-3 for y in ylim]
        ax.set_ylim(ylim)
        if log:
            # minval = np.floor(np.min(Z, axis=(0,1)))
            # maxval = np.ceil(np.max(Z,axis=(0,1)))
            # formatter = LogFormatter(10, labelOnlyBase=True) 
            
            # cbar.ax.yaxis.set_ticks(vals)
            # cbar.ax.yaxis.set_ticklabels(10**vals)
            cbar.ax.yaxis.set_ticks((1e-8,1e-7,1e-6,1e-5,1e-4))
        else:
            cbar = self.fig_free.colorbar(cm)
        #cbar.ax.set_ylabel('Free Electron Density, Å$^{-3}$', rotation=270,labelpad=20)
        #cbar.ax.set_ylabel('Energy density (eV/Å$^{3}$)', rotation=270,labelpad=20)
        cbar.ax.set_ylabel('Energy density (arb. u.)', rotation=270,labelpad=20)

        # plot the intensity
        ax2 = ax.twinx()
        ax2.plot(T, self.intensityData[::every],'w:',linewidth=1)

        ax2.get_yaxis().set_visible(False)
        ax2.set_ylim([0,ax2.get_ylim()[1]-ax2.get_ylim()[0]])

    # Plots a single point in time.
    def initialise_step_slices_ax(self):
        self.ax_steps = self.get_next_ax()
        self.fig_steps = self.fig
    
    #TODO get rid of hacky multiply_by_netural_state.
    def plot_step(self, t, prefactor_function = None, prefactor_power = 1, prefactor_args = {}, normed=True, fitE=None, multiply_by_neutral_state=None, **kwargs):        
        self.ax_steps.set_xlabel('Energy (eV)')
        #self.ax_steps.set_ylabel('$f(\\epsilon) \\epsilon $') # \\Delta \\epsilon is implied now. Want to distinguish from Hau-Riege whose f(e) is our f(e)e
        self.ax_steps.set_ylabel('Energy density (eV/$\\AA^{3}$)') # \\Delta \\epsilon is implied now. Want to distinguish from Hau-Riege whose f(e) is our f(e)e
        if self.use_electron_density:
            self.ax_steps.set_ylabel('Electron density ($\\AA^{-3}$)')#self.ax_steps.set_ylabel('$f(\\epsilon)')
        self.ax_steps.loglog()
        # norm = np.sum(self.freeData[n,:])
        n = self.timeData.searchsorted(t)
        data = self.freeData[n,:]
        X = self.energyKnot

        if normed:
            tot = self.get_density(t)
            data /= tot
            data/=4*3.14
        
        prefactor = 1
        if prefactor_function != None or prefactor_power not in [1,None]:
            assert prefactor_function != None and prefactor_power != None
            prefactor_function = np.vectorize(prefactor_function,otypes=[np.double])
            prefactor = prefactor_function(X,**prefactor_args)
            data = data[prefactor!=0]
            X = X[prefactor!=0]
            prefactor = prefactor[prefactor!=0]
            prefactor **= prefactor_power
        if multiply_by_neutral_state:
            #Multiply by proportion of states that are neutral
            self.aggregate_charges(False)
            t_idx = self.timeData.searchsorted(t)
            prefactor*=self.chargeData[multiply_by_neutral_state][t_idx,0]/np.sum(self.chargeData[multiply_by_neutral_state][0])
        density_factor = X # energy density
        if self.use_electron_density:
            density_factor = 1 

        return self.ax_steps.plot(X, prefactor*data*density_factor, label='%1.1f fs' % t, **kwargs)

    def plot_the_knots(self,times,vert_anchors,colours,padding=0.01):
        assert len(self.grid_update_time_Data == len(self.grid_point_Data))
        ylims = self.ax_steps.get_ylim()
        xlims = self.ax_steps.get_xlim()
        # Add a background for the knots
        tangle_top = np.e**(np.log(ylims[0]) + ( np.max(vert_anchors) )*(np.log(ylims[1]) - (1+padding)*np.log(ylims[0])))
        tangle_bot = np.e**(np.log(ylims[0]) - ( np.min(vert_anchors) )*(np.log(ylims[1]) - (1+padding)*np.log(ylims[0])))
        under_tangle =  patches.Rectangle((xlims[0]*0.1,tangle_bot), xlims[1], tangle_top, lw=0, edgecolor='none', facecolor='white',alpha=0.8,zorder=98)
        self.ax_steps.add_patch(under_tangle)
        # Plot the knots
        for t,y_anchor,col in zip(times,vert_anchors,colours):
            if self.grid_update_time_Data[-1] <= t:
                update_idx = len(self.grid_update_time_Data) - 1
            else:
                update_idx = self.grid_update_time_Data.searchsorted(t+1e-12) -1 # we increase the time slightly, as the indices of the knot updates correspond to the first step in the new basis.
            knots_to_plot = self.grid_point_Data[update_idx]
            y = [np.e**(np.log(ylims[0]) + y_anchor*(np.log(ylims[1]) - np.log(ylims[0])))]* len(knots_to_plot)
            self.ax_steps.scatter(knots_to_plot,y,color = col,s=80,zorder=99,lw=0.8, marker="|") #marker="|" or 10 work well.       
    def plot_fit(self, t, fitE, normed=True, **kwargs):
        t_idx = self.timeData.searchsorted(t)
        fit = self.energyKnot.searchsorted(fitE)
        data = self.freeData[t_idx,:]
        if normed:
            tot = self.get_density(t)
            data /= tot
            data/=4*3.14

        density_factor = self.energyKnot # energy density
        if self.use_electron_density:
            density_factor = 1 

        Xdata = self.energyKnot[:fit]
        Ydata = data[:fit]
        mask = np.where(Ydata > 0)
        T, n = fit_maxwell(Xdata, Ydata)
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, T, n)*density_factor,
            '--',label='%3.1f eV' % T, **kwargs)

    def plot_maxwell(self, kT, n, **kwargs):
        density_factor = self.energyKnot # energy density
        if self.use_electron_density:
            density_factor = 1         
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, kT, n)*density_factor,
            '--',label='%3.1f eV' % kT, **kwargs)


    def get_temp(self, t, fitE):
        t_idx = self.timeData.searchsorted(t)
        fit = self.energyKnot.searchsorted(fitE)
        Xdata = self.energyKnot[:fit]
        Ydata = self.freeData[t_idx,:fit]
        T, n = fit_maxwell(Xdata, Ydata)
        return (T, n)

    def get_density(self, t):
        t_idx = self.timeData.searchsorted(t)
        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2]) 
        de = de [1:] - de[:-1]
        return np.dot(self.freeData[t_idx, :], de)
    
    # Seems to plot ffactors through time, "timedata" is where the densities come from - the form factor text file corresponds to specific configurations (vertical axis) at different k (horizontal axis)
    # In any case, my ffactor function is probably bugged since it doesn't match this. - S.P.
    def plot_ffactor_get_R_AC4DC(self, a, num_tsteps = 10, timespan = None, show_avg = True, plot=True,**kwargs):

        if timespan is None:
            timespan = (self.start_t,self.end_t)#(self.timeData[0], self.timeData[-1])

        start_idx = self.timeData.searchsorted(timespan[0])
        stop_idx = self.timeData.searchsorted(timespan[1])

        ff_path = path.abspath(path.join(__file__ ,"../../output/"+a+"/Xsections/Form_Factor.txt"))
        fdists = np.genfromtxt(ff_path)
        # These correspond to the meaning of the FormFactor.txt entries themselves
        KMIN = 0
        KMAX = 2
        dim = len(fdists.shape)
        kgrid = np.linspace(KMIN,KMAX,fdists.shape[0 if dim == 1 else 1])
        if plot:
            fig2 = plt.figure()
            ax = fig2.add_subplot(111)
            ax.set_xlabel('$u$ (spatial frequency, atomic units)')
            ax.set_ylabel('Form factor (arb. units)')
        
        timedata = self.boundData[a][:,:-1] # -1 excludes the bare nucleus
        dynamic_k = np.tensordot(fdists.T, timedata.T,axes=1)   # Getting all k points? This has equal spacing -S.P. 
        step = (stop_idx - start_idx) // num_tsteps
        if plot:
            cmap=plt.get_cmap('plasma')
            fbar = np.zeros_like(dynamic_k[:,0])

        n=0
        times_used = []
        for i in range(start_idx, stop_idx, step):
            times_used.append(self.timeData[i])
            if plot:
                ax.plot(kgrid, dynamic_k[:,i], label='%1.1f fs' % self.timeData[i], color=cmap((i-start_idx)/(stop_idx - start_idx)))
            fbar += dynamic_k[:,i]
            n += 1

        print("Times used:",times_used)
        fbar /= n
        if show_avg and plot:
            ax.plot(kgrid, fbar, 'k--', label=r'Effective Form Factor')
        freal = dynamic_k[:,0]  # Real as in the real structure? - S.P.
        print("Warning, does not account for altering time step sizes.")
        print("R = ", np.sum(np.abs(fbar - freal))/np.sum(freal))    # Struct fact defn. of R (essentially equivalent to intensity definition)
        print("freal",freal)
        print("fbar",fbar)
        freal /= np.sum(freal)
        fbar /= np.sum(fbar)
        print("Normed R = ", np.sum(np.abs(fbar - freal))/np.sum(freal))
        # print("freal",freal)
        # print("fbar",fbar)
        return (fig2, ax)
        

def fit_maxwell(X, Y):
    guess = [200, 12]
    # popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess, sigma=1/(X+10))
    popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess)
    return popt

def maxwell(e, kT, n):
    if kT < 0:
        return 0 # Dirty silencing of fitting error - note we get negative values from unphysical oscillations, so this increases the average value around this point. -S.P.
    return n * np.sqrt(e/(np.pi*kT**3)) * np.exp(-e/kT)

def lnmaxwell(e, kT, n):
    return np.log(n) + 0.5*np.log(e/np.pi*kT**3) - e /kT

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Slater form factor
# TODO array of form factors for each state for each atom in species and each mom. transfer q


# Warning: Only for states with exclusively s, p orbitals
class SlaterShielding:
    # Z:  nuclear charge
    # shell_config: The number of electrons in each shell, passed as a list e.g. [1,8,8]
    def __init__(self,Z):
        self.Z = Z 
    def get_ff(self,occ_indices,k, occ_dict):
        
        all_possible_occs = occ_dict.values()
        ff_shape = () 
        if type(occ_indices) == np.ndarray:
            num_atoms = occ_indices.shape[0]
            num_times = occ_indices.shape[1]
            ff_shape += (num_atoms,num_times,)  # [atoms,times]       
        if type(k) == np.ndarray:
            ff_shape += k.shape  # [atoms,times, qx, qy]
        ff = 0
        if ff_shape != ():
            ff = np.zeros(ff_shape)
        if type(ff) == np.ndarray:
            for occ in occ_dict.keys():
                s = self.s_p_slater_shielding(occ,occ_dict)
                ff[occ_indices == occ] = self.calculate_config_ff(occ,k,s,occ_dict)[None,None]     # shell_occs -> [atoms,times,qx,qy]
        else:
            s = self.s_p_slater_shielding(occ_indices,occ_dict)
            ff = self.calculate_config_ff(occ_indices,k,s,occ_dict)
        return ff

    def s_p_slater_shielding(self,occ_index,occ_dict):
        shells_occ = occ_dict[occ_index]
        s_list = []
        for n in range(len(shells_occ)): 
            s = 0
            for m, N_e in enumerate(shells_occ): # N_e = num (screening) electrons
                if m > n:
                    continue
                if m == n:
                    s_factor = 0.35
                    if n == 0:
                        s_factor = 0.30
                    s += s_factor * max(0, N_e - 1)
                elif m == n - 1:
                    s += 0.85*N_e
                else:
                    s += N_e           
            s_list.append(s) 
        return s_list # s for each subshell
                    
    def get_shell_ff(self, k, shell_num,s):    #dim = scalar or [num_atoms]
        lamb = self.Z - s[shell_num-1]
        lamb/=shell_num
        D = 4*lamb**2+np.power(k,2)
        if shell_num == 1:  
            return 16*lamb**4/np.power(D,2)
        if shell_num == 2:
            return 64*lamb**6*(4*lamb**2-np.power(k,2))/np.power(D,4)
        if shell_num == 3:
            return 128/6*lamb**7*(192*lamb**5-160*lamb**3*np.power(k,2)+12*lamb*np.power(k,4))/np.power(D,6)
        else:
            raise Exception("ERROR! shell_num = "+str(shell_num))
        
    # Gets a configuration's form factor for array of k. 
    def calculate_config_ff(self,occ_index, k,  s,occ_dict):
        # We calculate the form factor for each subshell, then we multiply by the corresponding subshell occupancies to get the form factor.
        debug_old_ff = "0"
        occ = occ_dict[occ_index]
        num_subshells = len(occ)
        ff = 0
        for i in range(num_subshells):  
            ff += self.get_shell_ff(k,i+1,s)*occ[i]
            # check atomic ff is below Z.
            atomic_ff = "{:e}".format(np.array(ff).flatten()[-1])
            # Sanity check: each electron contributes at most 1 to the form factor (corresponding to free electron scattering)
            if np.array(ff).flatten()[-1] > self.Z:
                shell_ff = "{:e}".format(np.array(self.get_shell_ff(k,i+1,s)).flatten()[-1])
                if type(occ[i]) == int:
                    shell_occ = "{:e}".format(occ[i])
                else:
                    shell_occ = str(np.array(occ.astype(float)).reshape(-1,np.array(occ.astype(float)).shape[-1])[-1][i])
                raise Exception("Form factor (sample) of " + debug_old_ff + " + " + shell_ff +  " x " + shell_occ + " = " + atomic_ff + " above atomic number " + str(self.Z))
            debug_old_ff = atomic_ff
        return ff

if __name__ == "__main__":
    pl = Plotter(sys.argv[1])
    # pl.plot_free(log=True,min=1e-7)
    # pl.plot_all_charges()
    plt.show()
