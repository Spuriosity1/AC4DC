#%%
#TODO
## Important
# - make it so reflections don't overwrite same orientation, as stochastic now.
# - should have option to average out same miller indices be averaged out.
# - normalise damaged I rather than using neutze-style k factor in R factor.
## Not so important 
# - implement rhombic miller indices as the angle is actually 120 degrees on one unit cell lattice vector (or just do SPI)
# - Select times and integrate snapshots of intensities with gaussian quadrature or some better method than trapezoid method. (I started this with scatter_quad.py)

'''
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================
(C) Spencer Passmore 2023
'''

import os
os.getcwd()
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
######

import os.path as path
from Bio.PDB.vectors import Vector as Bio_Vect
from Bio.PDB.vectors import homog_trans_mtx, set_homog_trans_mtx
from Bio.PDB.vectors import rotaxis2m
#from Bio.PDB.PDBParser import PDBParser
#from Bio.PDB.PDBIO import PDBIO
#from Bio.PDB.StructureBuilder import StructureBuilder
from xpdb import sloppyparser as xPDBParser
from xpdb import SloppyPDBIO as xPDBIO
from xpdb import SloppyStructureBuilder as xStructureBuilder  # Hack, enables atom counts over 10,000
from Bio.PDB.Atom import Atom as PDB_Atom
#from sympy.utilities.iterables import multiset_permutations
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import cos
from numpy import sin
from plotter_core import Plotter
from scipy.spatial.transform import Rotation as Rotation
from matplotlib.colors import to_rgb
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
import copy
import pickle
import colorcet as cc; import cmasher as cmr
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.offline as pltly_offline
from IPython.display import display, HTML
from IPython import get_ipython
interactive = True
if interactive and __name__ == "__main__":
    get_ipython().run_line_magic('colors', 'nocolor')
    get_ipython().run_line_magic('matplotlib', 'widget')
    # pltly_offline.init_notebook_mode()
    # display(HTML(
    #     '<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_SVG"></script>'
    # ))    

plt.ioff()  # stops weird vscode stuff

DEBUG = False 
SEEDED = False# TODO fully implement for all random stuff
c_au = 137.036; eV_per_Ha = 27.211385; ang_per_bohr = 1/1.88973  # > 1 ang = 1.88973 bohr

RESULTS_LOCAL_PATH = "results/"

class Results():
    def __init__(self,num_points,image_index):
        self.phi = np.zeros(num_points)
        self.phi_aligned = np.zeros(num_points)
        self.I = np.zeros(num_points)
        self.q = np.zeros(num_points)
        self.X = np.zeros(num_points)
        self.image_index = image_index 

    def package_up(self,miller_indices):
        _, self.phi_mesh = np.meshgrid(self.q,self.phi)  
        self.X, self.phi_mesh = np.meshgrid(self.X, self.phi)
        self.q, self.phi_aligned_mesh = np.meshgrid(self.q,self.phi_aligned) 
        #self.q_scr, self.phi_aligned_mesh = np.meshgrid(self.q_scr,self.phi_aligned) 
        self.miller_indices = miller_indices 
    def diff(self,other):
        '''
        Get single-point R factors between two images
        '''
        self.R = np.zeros(self.I.shape)
        for i in range(len(self.q)):
            subtracted = False
            for j in range(len(other.q)):
                if self.phi[i] == other.phi[j] and self.q[0][i] == other.q[0][j]:   # Using phi is equivalent to using phi_aligned, as this function is only used for same-orientation comparisons.
                    #self.I[i] = np.abs(self.I[i]- other.I[j])/self.I[i] # Intensity
                    self.R[i] = np.abs(np.sqrt(self.I[i])- np.sqrt(other.I[j]))/np.sqrt(max(self.I[i],other.I[j]))# form factor  # TODO not sure why later times is sometimes larger but probs interference thing. Hopefully will disappear when we use large time gaps... or at least when we average over to get R factor.
                    subtracted = True
            if subtracted == False:
                self.R[i] = -1 
class Results_SPI():
    pass
class Results_Grid():
    pass    

class Crystal():
    def __init__(self, pdb_fpath, allowed_atoms, positional_stdv = 0, is_damaged=True, include_symmetries = True, rocking_angle = 0.3, cell_packing = "SC", CNO_to_N = False, supercell_scale = 1,num_supercells=1, supercell_simulations = 1, S_to_N=False,convert_excluded_elements_to_N=False):
        '''
        rocking_angle [degrees]
        cell_packing ("SC","BCC","FCC","FCC-D")
        '''
        self.cell_packing = cell_packing
        self.rocking_angle = rocking_angle * np.pi/180            
        self.is_damaged = is_damaged
        self.supercell_scale = supercell_scale  # TODO allow for non-cubic crystals and non SC cell packing.
        self.num_supercells = num_supercells
        self.supercell_simulations = supercell_simulations
        if positional_stdv == 0 and is_damaged == False:
            self.supercell_simulations = 1
        self.pdb_fpath = pdb_fpath
        self.positional_stdv = positional_stdv/ang_per_bohr # RMS error in coord positions, designed for SPI sim only but I guess it wouldn't be detrimental for crystal sim.
        
        assert self.supercell_simulations <= num_supercells
        ## Parameters to be parsed by custom function because I cannot understand Bio.PDB.PDBParser's documentation.
        self.sym_rotations = []; self.sym_translations = [];   # Symmetry for each asymmetric unit simulated. If non-SPI, this is defining the supercell.
        self.cell_dim = None   # unit cell length parameters. Angles not implemented yet.        
        self.parse_data_from_pdb() # All asymmetric units in unit cell

        if not include_symmetries:
            # Ignore parsed symmetries and use a single asymmetric unit per cell.
            self.sym_rotations = []; self.sym_translations = []; 
            self.add_symmetry_to_cells(np.identity(3),np.zeros((3)),"(X,Y,Z)")

        self.supercell_dim = self.cell_dim*supercell_scale 

        ## Dictionary for going from pdb to ac4dc names.
        # different names
        PDB_to_AC4DC_dict = dict(
            NA = "Sodion", CL = "Chloride",
        )
        # Same names
        for elem in ["H","He","C","N","O","P","S","Gd"]:  # pdb names
            PDB_to_AC4DC_dict[elem] = elem   # ac4dc name
        # Modify the values in the dictionary according to arguments.
        for k,v in PDB_to_AC4DC_dict.items():
            # Light atom approximation. 
            if CNO_to_N:
                if v in ["C","O"]: 
                    v = "N"        
            if S_to_N:
                if v == "S":
                    v = "N"  
            if convert_excluded_elements_to_N and v not in allowed_atoms: 
                    v = "N"    
            PDB_to_AC4DC_dict[k] = v       
        # Get structure using Bio.PDB's parser
        parser=copy.deepcopy(xPDBParser)
        structure_id = os.path.basename(self.pdb_fpath)
        structure = parser.get_structure(structure_id, self.pdb_fpath)     
        # Get those cheeky charge clusters
        species_dict = {}
        pdb_atoms = []
        pdb_atoms_ignored = ""
        #TODO change to just reading data folder and only excluding atoms that are specified, rather than requiring user to pass in all allowed atoms.       
        for ac4dc_atom in allowed_atoms:
            if ac4dc_atom + "_fast" in allowed_atoms:
                #TODO read data folder to resolve ambiguity.
                raise Exception("Ambiguity, both "+ac4dc_atom+" and "+ac4dc_atom+"_fast were given in allowed_atoms.") 
            if ac4dc_atom + "_faster" in allowed_atoms:
                #TODO read data folder to resolve ambiguity.
                raise Exception("Ambiguity, both "+ac4dc_atom+" and "+ac4dc_atom+"_faster were given in allowed_atoms.")            
        for atom in structure.get_atoms():
            # Get ze data
            R = atom.get_vector()
            name = PDB_to_AC4DC_dict.get(atom.element)
            # Pop into our desired format
            if name == None:
                pdb_atoms_ignored += atom.element + " "
                continue
            # TODO make this not suck.
            if name not in allowed_atoms: 
                name+= "_fast"
                if name not in allowed_atoms: 
                    name +="er"
                    if name not in allowed_atoms:
                        continue
            if name not in species_dict.keys():
                species_dict[name] = Atomic_Species(name,self) 
                pdb_atoms.append(atom.element)
            species_dict[name].add_atom(R)
        ac4dc_atoms_ignored = ""
        for string in allowed_atoms:
            if string not in species_dict.keys():
                ac4dc_atoms_ignored += string + " "
                continue
        print("The following atoms will be considered (AC4DC names ; pdb names):")
        for p,a in zip(pdb_atoms,[PDB_to_AC4DC_dict[x] for x in pdb_atoms]): 
            print("%-10s %10s" % (p, a))
        if pdb_atoms_ignored != "":
            print("The following pdb atoms were found but ignored:",pdb_atoms_ignored)
        if ac4dc_atoms_ignored != "":
            print("The following atoms were allowed but not found:",ac4dc_atoms_ignored)

        for species in species_dict.values():
            species.set_coord_deviation()
        self.species_dict = species_dict 


    def set_ff_calculator(self,ff_calculator):
        self.ff_calculator = ff_calculator                  
    
    def add_symmetry_to_cells(self,symmetry_factor,symmetry_translation,symmetry_label="Symmetry Operation:"):
        '''
        Takes in the symmetry factor and translation of the unit cell, and translates them across each unit cell
        '''
        symmetry_translation/= ang_per_bohr
        # Simple cubic packing. (actually it's all rectangular prisms, TODO)
        if self.cell_packing == "SC":  
            self.num_cells = self.supercell_scale**3
            # Generate the coordinates of the cube (performance: defining scaling matrix/cube coords outside of here would be more efficient). But only runs once  \_( '_')_/ ¯\_(ツ)_/¯.
            x, y, z= np.meshgrid(np.arange(0, self.supercell_scale), np.arange(0, self.supercell_scale), np.arange(0, self.supercell_scale))
            cube_coords = np.stack([ x.flatten(), y.flatten(), z.flatten()], axis = -1)
            # Construct the cube by adding the "unit cell translation" to the symmetry translation. 
            # (e.g. if crystal is centred on origin, in fractional crystallographic coordinates the index of the unit cell is equal to the unit cell's translation from the origin.) 
            print_thing = np.array(["][x]   [","][y] + [","][z]   ["],dtype=object)
            print(symmetry_label,"\n",np.c_[symmetry_factor, print_thing ,symmetry_translation],sep="")
            for coord in cube_coords:
                #print(coord)
                translation = coord*self.cell_dim + symmetry_translation
                self.sym_translations.append(translation)
                self.sym_rotations.append(symmetry_factor)
                 
        else:
            #Can just duplicate symmetries if not SPI since equivalent for crystal.
            raise Exception("Lacking implementation")
            
    def parse_data_from_pdb(target):
        '''
        Parses the symmetry transformations from the pdb file into ndarrays and adds each to the crystal.  
        No unit operations are performed.
        Also returns the parsed matrices. 
        '''
        at_SYMOP = False
        at_symmetry_xformations = False
        sym_factor = np.zeros((3,3))
        sym_trans = np.zeros((3))
        sym_labels = []
        sym_mtces_parsed = []        
        with open(target.pdb_fpath) as pdb_file:
            for line in pdb_file:
                line = line.strip()
                # End data - Check if left section of data; flagged by the line containing solely "REMARK ###". 
                if len(line) <= 10: 
                    if at_SYMOP:
                        at_SYMOP = False
                    if at_symmetry_xformations: 
                        at_symmetry_xformations = False
                if line[0:4] == "ATOM":
                    break
                # Symmetry operators (purely for terminal output - doesn't affect simulation)
                if at_SYMOP:
                    sym_labels.append(line[17:21]+": ("+line[21:].strip()+")")
                # Add each symmetry to target
                if at_symmetry_xformations:
                    entries = line.split()[2:]
                    x = int(entries[0][-1]) - 1  
                    entries = [float(f) for f in entries[2:]]
                    # x/y index rows/cols
                    for y, element in enumerate(entries[:-1]):
                        sym_factor[x,y] = element
                    sym_trans[x] = entries[-1]
                    if x == 2:
                        sym_mtces_parsed.append([sym_factor.copy(),sym_trans.copy()])
                        

                ## inline data
                if line[0:6] == "CRYST1":
                    if target.cell_dim != None:
                        raise Exception("Cannot handle multiple crystal entries")
                    entries = line.split()[1:]
                    target.cell_dim = [float(a) for a in entries[0:3]]
                    target.cell_dim = np.array(target.cell_dim)/ang_per_bohr 

                ## Start data - Check if this line marks the next as the beginning of a desired data section.                    
                # Symmetry operators
                if "NNNMMM   OPERATOR" in line:
                    at_SYMOP = True
                    print("Parsing symmetry operations...")
                # symmetry matrices
                if line == "REMARK 290 RELATED MOLECULES.":
                    at_symmetry_xformations = True
        for sym_factor, sym_trans in sym_mtces_parsed:
            target.add_symmetry_to_cells(sym_factor.copy(),sym_trans.copy(),sym_labels.pop(0))
        return sym_mtces_parsed
    
    def get_sym_xfmed_point(self,R,symmetry_index):
        '''
        R = 1D or 2D array of shape (3,N)
        '''
        i = symmetry_index
        # Transpose because we've stored the vectors in the 0th axis. 
        return (self.sym_rotations[i] @ R.T).T+ self.sym_translations[i]   # dim = [xyz,xyz] x [xyz,N] or dim = [xyz,xyz] x [xyz,1]
            
    def save_structure(self,dir="targets"):
        '''
        Saves the full structure in a pdb file format for use with Solvate1.0   
        Stdv not included.
        I have not tested if using it for over 1e5 atoms (when xpdb alters things so it doesn't break) works when plugged into SOLVATE
        TODO need to add boilerplate (replace HETATM with 'ATOM  ', add symmetries and cell dimensions to start of file.)
        '''
        #Load it
        # parser=PDBParser(PERMISSIVE=1)
        # structure_id = os.path.basename(self.pdb_fpath)
        # structure = parser.get_structure(structure_id, self.pdb_fpath)       

        #Initialsie structure
        structure = xStructureBuilder()
        structure.init_structure("full_struct")
        structure.init_model("M")
        structure.init_seg("")

        # Add atoms
        serial_number = 1
        structure.init_chain("C")
        for i in range(len(self.sym_rotations)):
            residue_name = "L"+str(i); r_args = (" ",i+1,"r")
            print(structure.chain.child_dict)
            #residue_id = (r_args[0]+"_"+residue_name,r_args[1],r_args[2])  # use if set r_args[0] to "H"
            residue_id = r_args
            structure.init_residue(residue_name,*r_args)
            print(structure.chain.child_dict)
            residue = structure.chain[residue_id] 
            for species in self.species_dict.values():
                points = np.array(species.coords)                    
                coord_list = self.get_sym_xfmed_point(points,i).tolist()
                for j, coord in enumerate(coord_list):
                    coord=tuple([c*ang_per_bohr for c in coord])
                    name = ' '+species.name+str(j)+' '
                    residue.add(PDB_Atom(name= name, coord=coord, bfactor=0., occupancy=1., altloc=' ', fullname=name, serial_number=serial_number,element=species.name))                
                    serial_number+=1
        # TODO boiler plate
        # structure.set_symmetry("P 1", "1")
        
        # Save it
        io=xPDBIO()
        io.set_structure(structure.get_structure())  # StructureBuilder object is not Structure object
        fname = path.basename(self.pdb_fpath)[:-4]+"_full_struct.pdb"
        io.save(dir+"/"+fname)     

    def plot_me(self,max_points = 100000,water_index = None,**layout_kwargs):
        if water_index != None:
            assert self.supercell_scale == 1 and len(self.sym_rotations) == 1, "Plotting water with an added intra-cell or crystal symmetry is not supported."
        '''
        plots the symmetry points. Origin is the centre of the base asymmetric unit (not the centre of a unit cell).
        RMS error in positions ('positional_stdv') is not accounted for
        ''' 
        plt.close()
        view_width = 1000
        view_height = 800

        num_atoms_avail = 0
        for species in self.species_dict.values():
             num_atoms_avail += len(species.coords)
        assert water_index is None or water_index < num_atoms_avail
        num_test_points = max(1,int(max_points/(len(self.sym_translations))))
        num_test_points = min(num_atoms_avail,num_test_points)         
        test_points = np.empty((num_test_points,3))
        qualifier = "sample of"
        if num_atoms_avail == len(test_points):
            qualifier = "all" 
        print("Number of atoms in supercell:",num_atoms_avail*len(self.sym_translations))
        print("Generating plot with",qualifier,len(test_points),"atoms per asymmetric unit (plotting "+str(len(test_points)*len(self.sym_translations))+" in total)")
        i = 0
        for species in self.species_dict.values():
            if i >= num_test_points:
                break                
            for R in species.coords:
                test_points[i] = R
                i+=1
                if i >= num_test_points:
                    break      
        plot_coords = []
        for i in range(len(self.sym_rotations)):
            coord_list = self.get_sym_xfmed_point(test_points,i).tolist()
            plot_coords.extend(coord_list)  

        #Colors - which to highlight (root atom is first atom of cell)
        c_first_root_atom = False # red
        c_first_unit = True # pink
        c_first_cell = True # aquamarine
        c_all_root_atoms = False # yellow

        top_atom_index = 0 # index of atom that has highest z
        max_height = -np.inf
        for i,coord in enumerate(test_points):
            if coord[2] > max_height:
                max_height = coord[2]
                top_atom_index = i 
        if water_index == None:
            alpha = None
            color = np.empty(len(plot_coords),dtype=object) 
            color.fill('green')
            if c_first_cell:
                for i in range(int(len(self.sym_rotations)/self.num_cells)):           
                    color[i*self.num_cells*len(test_points):i*self.num_cells*len(test_points)+len(test_points)] = 'aquamarine' #'c'  # atoms in one unit cell  
            if c_first_unit:
                color[:len(test_points)] = 'pink'  # atoms in one asym unit       
            #Cursed but it works. 
           
            for i in range(self.num_cells):
                if c_all_root_atoms:
                    for j in range(int(len(self.sym_rotations)/self.num_cells)):
                        color[i*int(len(test_points)*len(self.sym_rotations)/self.num_cells) + j*len(test_points) + top_atom_index] = 'yellow'      # same atom in every asym unit.  
                if c_first_root_atom:
                    color[i*len(test_points)+top_atom_index] = 'red'      # same atom in every same unit cell         
        else:
            kernel_color = np.empty(water_index+1,dtype=object);# kernel_alpha =  np.empty(water_index+1,dtype=np.double)
            bg_color = np.empty(len(plot_coords)-(water_index+1),dtype=object); #bg_alpha = np.empty(len(plot_coords)-(water_index+1),dtype=np.double)
            kernel_color.fill('rgba(92,169,4,1)'); #kernel_alpha.fill(1) 
            bg_color.fill('rgba(0,30,255,0.2)');  # water atom colours
            color = np.concatenate((kernel_color,bg_color))
            #alpha = np.concatenate((kernel_alpha,bg_alpha))
        
        plot_coords = np.array(plot_coords)*ang_per_bohr # convert to angstrom
        raise_non_unique_exception = False
        if np.unique(plot_coords,axis=0).shape != plot_coords.shape:
            a, unique_indices = np.unique(plot_coords,axis=0,return_index = True)
            print("WARNING", len(np.delete(plot_coords, unique_indices)), "non-unique coords found:")
            raise_non_unique_exception = True
            print(np.delete(plot_coords, unique_indices))  # delete anyway.            
                            
            #color[color!='y'] = '#0f0f0f00'   # example to see just one atom type
        x_range = np.array([np.min(plot_coords[:,0]),np.max(plot_coords[:,0])])
        y_range = np.array([np.min(plot_coords[:,1]),np.max(plot_coords[:,1])])
        z_range = np.array([np.min(plot_coords[:,2]),np.max(plot_coords[:,2])])
        print("---")
        print("x_min/max:",x_range[0],x_range[1])
        print("y_min/max:",y_range[0],y_range[1])
        print("z_min/max:",z_range[0],z_range[1])
        print("---")
        max_len = 0
        min_len = np.inf
        for elem in x_range,y_range,z_range:
            max_len = max(max_len,abs(elem[1]-elem[0]))      
            min_len = min(min_len,abs(elem[1]-elem[0]))      

        titles = [{"text": ax+'$ [\AA]$'} for ax in ['x','y','z']]
        xaxis,yaxis,zaxis = [{'title': title} for title in titles]
        ##
        custom_axes = False #(i.e. we won't use the axis variables above)
        ranges = np.array([np.max(plot_coords[:,i]) - np.min(plot_coords[:,i]) for i in range(3)])
        if custom_axes:
            # data for custom axis centred in the space
            avg_range = np.sum(ranges)/3
            step_size = max(1, 2.5*10**(np.log10(avg_range/10) - ((np.log10(avg_range/10)-1)%1)))
            # Translate to be centred on axis ticks
            plot_coords -= np.array([np.mean(plot_coords[:,i]) for i in range(3)])        
            centre_x = centre_y = centre_z = 0
            #centre_x = np.mean(plot_coords[:,0]); centre_y = np.mean(plot_coords[:,1]); centre_z = np.mean(plot_coords[:,2])

            #centre_x -= centre_x%step_size; centre_y -= centre_y%step_size; centre_z -= centre_z%step_size 
        origin_on_corner = True  
        # asym unit
        # size = 5.5*max(1,500/(min_len+max_len))   
        # angular_aperture = np.pi*0.6 
        # dot_lw = 1
        # unit cell
        # size = 7.5*max(1,500/(min_len+max_len))   
        # angular_aperture = np.pi*0.6
        # dot_lw = 0.1 # unit
        #2x2x2 Crystal
        # size = 10*max(1,500/(min_len+max_len))   #solvated crystal
        # angular_aperture = np.pi*0.7 # Solvated crystal        
        # dot_lw = 0 # cryst    

        # tetrapeptide:
        # size = 3*max(1,500/(min_len+max_len))
        # angular_aperture = np.pi*0.6
        # dot_lw = 0.1 

        # non-camera perspective (dots all same):
        size = 1.5*max(1,500/(min_len+max_len))
        dot_lw = 0
        if origin_on_corner:
            max_num_ticks = 7
            minima = np.array([np.min(plot_coords[:,i]) for i in range(3)])
            
            #Translate origin to corner.
            plot_coords -= minima  
            minima = np.array([np.min(plot_coords[:,i]) for i in range(3)])
            maxima = np.array([np.max(plot_coords[:,i]) for i in range(3)])
            tick_spacing = 10
            addendums = [dict(tickfont=dict(size=15), showbackground=False,mirror="all", nticks=1+max(1,round(m*1.4/tick_spacing)), range=[0,m*1.4]) for m in maxima] #nticks=max(2,round(m/max(maxima)*max_num_ticks))
            addendums[1]["range"] = [0,maxima[1]] #pull back wall in
            addendums[1]["nticks"] = 1+max(1,round(maxima[1]/tick_spacing))
            for i, axis in enumerate([xaxis,yaxis,zaxis]):
                axis |= addendums[i]
            
            # INCREDIBLY LAZY implementation to provide some sense of distance in png snapshots.
            hack_snapshots = False
            if hack_snapshots:
                # Calculate size relative to view from camera that is capturing full width of structure.
                # Hacky as I am assuming camera is at a certain position.
                r = size  # some arbitrary 'radius' of visible atom that is equivalent to the largest possible pixel width. 
                #width = view_width/r
                # a guess (increase to make fall off more slowly)

                target_length = np.sqrt(np.sum([range**2 for range in ranges]))
                distance_to_full_capture = np.tan(angular_aperture/2)*(target_length/2)
                #psi = -np.pi/4; thet=np.pi/4 #Camera position relative to target's centre
                psi = -np.pi/2; thet=np.pi/6 # cryst alternate angle
                camera_pos = [distance_to_full_capture*np.sin(thet)*np.cos(psi),distance_to_full_capture*np.sin(thet)*np.sin(psi),distance_to_full_capture*np.cos(thet)]
                max_atom_angle = np.arctan(r/(target_length/2))
                atom_distances = [np.sqrt(np.sum([(val-camera_pos[i])**2 for i, val in enumerate(coord)])) for coord in plot_coords]
                #r*angle/max_angle
                size = [r*np.arctan(r/(distance_to_full_capture + d))/max_atom_angle for d in atom_distances]
        #size = max(1,500/(min_len+max_len))  


        scatter_points = go.Scatter3d(
            x=plot_coords[:,0], 
            y=plot_coords[:,1], 
            z=plot_coords[:,2], 
            marker=go.scatter3d.Marker(color=color,size=size,line=dict(width=dot_lw,color='black'),opacity=1), 
            mode='markers',
        )
        fig=go.Figure(data=scatter_points)
        # def sphere(x, y, z, radius, resolution=4):
        #     """Return the coordinates for plotting a sphere centered at (x,y,z)"""
        #     u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
        #     X = radius * np.cos(u)*np.sin(v) + x
        #     Y = radius * np.sin(u)*np.sin(v) + y
        #     Z = radius * np.cos(v) + z
        #     return (X, Y, Z)
        # scatter_spheres = []
        # for point in plot_coords:
        #     (x_pns_surface, y_pns_surface, z_pns_surface) = sphere(*point,max(1,500/(min_len+max_len)))
        #     scatter_spheres.append(go.Surface(x=x_pns_surface, y=y_pns_surface, z=z_pns_surface, opacity=0.5))
        #fig=go.Figure(data=scatter_spheres)
        # Setup 3D scene stuff
        aspect = np.empty(3)
        for i, elem in enumerate((x_range*1.4,y_range,z_range*1.4)):
            aspect[i] = abs(elem[1]-elem[0])/max_len    
        zoom = 1.8 # inital zoom      
        fig.update_layout(
            margin=dict(l=20,r=20,t=20,b=20),
            width=view_width, height=view_height,
            scene = dict(
                xaxis = xaxis,
                yaxis = yaxis,
                zaxis = zaxis,
                aspectratio = dict(x=aspect[0]*zoom,y=aspect[1]*zoom,z=aspect[2]*zoom),
                #camera=dict(eye=dict(x=1,y=0,z=0.6))
            ),
        ) 
        if custom_axes:
            # Remove default axis
            axis_args = dict(showgrid = False, zeroline = False, showticklabels = False, title = dict(text = ""))
            fig.update_layout(scene = dict(zaxis=axis_args,yaxis=axis_args,xaxis=axis_args))      
            # Add custom axis centred in the space
            x_tickvals = np.append(
                np.flip(np.arange(centre_x,np.min(plot_coords[:,0])+abs(np.min(plot_coords[:,0]))%step_size-step_size*3/4,-step_size)),
                np.arange(centre_x,np.max(plot_coords[:,0])-abs(np.max(plot_coords[:,0]))%step_size+step_size*3/4,step_size)
            )
            y_tickvals = np.append(
                np.flip(np.arange(centre_y,np.min(plot_coords[:,1])+abs(np.min(plot_coords[:,1]))%step_size-step_size*3/4,-step_size)),
                np.arange(centre_y,np.max(plot_coords[:,1])-abs(np.max(plot_coords[:,1]))%step_size+step_size*3/4,step_size)
            )
            z_tickvals = np.append(
                np.flip(np.arange(centre_z,np.min(plot_coords[:,2])+abs(np.min(plot_coords[:,2]))%step_size-step_size*3/4,-step_size)),
                np.arange(centre_z,np.max(plot_coords[:,2])-abs(np.max(plot_coords[:,2]))%step_size+step_size*3/4,step_size)
            )
            x_tickvals = np.round(x_tickvals,0)
            y_tickvals = np.round(y_tickvals,0)
            z_tickvals = np.round(z_tickvals,0)
            print(np.flip(np.arange(centre_y,np.min(plot_coords[:,1])+abs(np.min(plot_coords[:,1]))%step_size,-step_size)))
            line_width = 10
            marker_size = 3
            fontsize = 20
            xaxis_line =go.Scatter3d(
                            x = x_tickvals,
                            y = (centre_y,)*len(x_tickvals),
                            z = (centre_z,)*len(x_tickvals),
                            mode = "lines+markers+text",
                            marker = dict(size=marker_size),
                            line = dict(width = line_width),
                            text=x_tickvals,
                            textfont=dict(size=fontsize),
                            )
            yaxis_line =go.Scatter3d(
                            y = y_tickvals,
                            z = (centre_z,)*len(y_tickvals),
                            x = (centre_x,)*len(y_tickvals),
                            mode = "lines+markers+text",
                            marker = dict(size=marker_size),
                            line = dict(width = line_width),
                            text=y_tickvals,
                            textfont=dict(size=fontsize),
                            )
            zaxis_line =go.Scatter3d(
                            z = z_tickvals,
                            x = (centre_x,)*len(z_tickvals),
                            y = (centre_y,)*len(z_tickvals),
                            mode = "lines+markers+text",
                            marker = dict(size=marker_size),
                            line = dict(width = line_width),
                            text=z_tickvals,
                            textfont=dict(size=fontsize),
                            )                
            centre_ball = go.Scatter3d(x = (centre_x,centre_x), 
                                    y = (centre_y,centre_y), 
                                    z = (centre_z,centre_z), 
                                    mode = "markers", 
                                    hoverinfo = "skip", 
                                    marker = dict(size = 8))
            fig.add_traces([xaxis_line,yaxis_line,zaxis_line,centre_ball])

        # Update the kwargs separately so that we can call arguments used above without conflict.
        fig.update_layout(**layout_kwargs)
        fig.show()       
        #fig.write_html("Structure.html")

        if raise_non_unique_exception:
            raise Exception(str(len(np.delete(plot_coords, unique_indices))) + " non-unique coords found:")


class Atomic_Species():
    def __init__(self,name,crystal):
        self.name = name 
        self.crystal = crystal 
        self.ff = 0 # form_factor
        self.coords = []  # coord of each atom in species in asymmetric unit

    def add_atom(self,vector):
        '''
        This function adds an atom to the asymmetric unit of the crystal. 
        We do not store additional coordinates, instead storing the symmetries, and an array of atomic states corresponding to each atom, for each symmetry. (so num symmetries * num atoms added)
        '''           

        self.coords.append(vector.get_array()/ang_per_bohr)
        
    def set_stochastic_states(self):
        '''
        We set a state for each atom, including symmetries, so that different q applied to the same atom at the same time corresponds to the same state. 
        When we are finding the atomic form factors, we call get_stochastic_f() on these states.
        The dimensionless nature of the model forces us to make the dubious approximation that an atom's state is independent of its prior states.
        With a hybrid molecular dynamics model informed by AC4DC, the nuclei's states could be tracked properly throughout time, and this function would be replaced
        by a call to the data of the atomic nuclei's states.
        '''
        print("Creating time-varying states for atom "+self.name+" from plasma simulation's data")
        self.times_used = self.crystal.ff_calculator.get_times_used()
        if self.num_atoms != len(self.crystal.sym_rotations)*len(self.coords):
            raise Exception("num atoms was not same on set_stochastic_states call as when set by set_coord_deviation")
        if self.crystal.is_damaged:
            # Initialise an array that tracks the individual atoms' form factors.
            orb_occs_shape = (self.num_atoms,len(self.times_used))  # [num atoms,times]
            # TODO instead of storing lists, replace with indices and a list with corresponding states. Also use index to get ff rather than orbocc list
            self.orb_occs= np.empty(orb_occs_shape,dtype=list)    # self.orb_occs[i] is an array of states corresponding to each time. We make the necessary approximation that an atom's state is independent of its prior states. This approximation is dubious at low unit cell numbers, but at higher numbers, because the contribution from an atom at the same relative cell coordinate and state as another atom will be equivalent, we get the same outcome so long as the probability distribution of states is representative of the actual distribution of states. i.e. tracking state history at the same global position is redundant at high unit cell counts where we can expect the distribution of states at a given coordinate to have a low deviation between species.  
            for idx in range(self.num_atoms):
                seed = None
                if SEEDED:
                    seed = idx
                self.orb_occs[idx],_dummy,self.orb_occ_dict = self.crystal.ff_calculator.random_state_snapshots(self.name,seed) 
        else:
           self.ground_state = self.crystal.ff_calculator.get_ground_state_shells(self.name)           
    def set_coord_deviation(self):
        self.num_atoms = len(self.crystal.sym_rotations)*len(self.coords)
        self.error = np.empty((self.num_atoms,3))
        for idx in range(self.num_atoms):
            # get random error in spherical coords based on RMS error in position, convert to cartesian.
            # there's probably a better way to do it
            err_phi,err_thet = np.random.random()*2*np.pi, np.random.random()*np.pi
            err_r = np.random.normal(0,self.crystal.positional_stdv, size = (3))
            self.error[idx] = err_r*(np.sin(err_phi)*np.cos(err_thet),np.sin(err_phi)*np.cos(err_thet),np.cos(err_thet))


    def get_stochastic_f(atom_idx,q_arr):
        '''
        (Is defined by set_scalar_form_factor).
        Returns the form factor multiplied by sqrt(I). f.shape = ( len(times) , ) + momenta.shape  
        
        atom_idx, int or int array
        q_arr = mom. transfer [1/a0], scalar or array
        '''
        raise Exception("Did not set_scalar_form_factor before calling stochastic f")
        
    def set_scalar_form_factor(self,stochastic=True):
        '''
        
        '''
        # F_i(q) = f(q)*T_i(q), where f = self.ff is the time-integrated average 
        if not stochastic:
            pass
            #self.ff = self.crystal.ff_calculator.f_average(q_arr,self.name)      # note that in taking the integral to get this ff, we included the relative intensity.
        else:
            # Undamaged case, no stochastic dynamics.
            if not self.crystal.is_damaged: 
                def tmp_func(atom_idx,q_arr): 
                    return self.crystal.ff_calculator.f_undamaged(q_arr,self.name,self.ground_state)[0]
            # Damaged, we 
            else:
                def tmp_func(atom_idx,q_arr): 
                    return self.crystal.ff_calculator.random_states_to_f_snapshots(self.times_used,self.orb_occs[atom_idx],q_arr,self.name,self.orb_occ_dict)[0]  # f has form  [times,momenta]
            self.get_stochastic_f = tmp_func
class XFEL():
    def __init__(self, experiment_name, photon_energy, detector_distance_mm=100, q_minimum = None, q_cutoff = None, max_miller_idx = None, screen_type = "hemisphere", num_orients_crys=1, orientation_axis_crys = None, x_orientations = 1, y_orientations = 1, pixels_per_ring = 400, num_rings = 50,t_fineness=100,SPI_y_rotation = 0,SPI_x_rotation = 0,SPI_z_rotation = 0,all_miller_indices=False, custom_cell_dims_for_miller_indices=None,override_max_q = False):
        """ #### Initialise the imaging experiment's controlled parameters
        experiment_name:
            String that the output folder will be named.        
        photon_energy [eV]:
            Should be the same as that given in the original input file!
        detector_distance_mm [mm];
            The distance in mm between the target and the centre of the detector 
        screen_type:
            "circle", "hemisphere", "sphere". Circle corresponds to a flat screen - i.e. it 'squashes ya dots'. Sphere makes no difference if a max q is specified below the hemisphere q range.
        pixels_per_ring: 
            The number of different values for phi to plot. phi is the angle that the incident light makes with 
            the y axis, (z-axis is firing line. y-axis is any arbitrary axis perpendicular to z-axis.)
        num_rings: 
            Determines the number of q to calculate. Note pixels are currently just points
        t_fineness:
            Number of time steps to calculate.
        alpha,beta,gamma:
            Angle of rotation about the x, y, and z axis.
        orientation_set:
            contains each set of cardan angles for each orientation to images. Crystal only. 
            Overriden for imagings with random orientations. TODO replace x_orientations and y_orientations with this.
        max_miller_idx:
            If specified, the maximum momentum transfer is set to that corresponding to (m,m,m), where m = max_miller_idx
        q_cutoff [1/angstrom]:
            If specified, only bragg points corresponding to momentum transfers at or below this value will be simulated.    
        """
        self.experiment_name = experiment_name
        self.detector_distance = detector_distance_mm*1e7/ang_per_bohr  # [converts to a0 (bohr)]
        self.photon_momentum = 2*np.pi/E_to_lamb(photon_energy)  # atomic units, [a0^-1]. (2pi=h)
        self.photon_energy = photon_energy  #eV Attention: not in atomic units
        self.pixels_per_ring = pixels_per_ring
        self.num_rings = num_rings
        self.x_orientations = x_orientations
        self.y_orientations = y_orientations
        self.max_miller_idx = max_miller_idx
        self.all_miller_indices = all_miller_indices
        self.custom_cell_dims_for_miller_indices = custom_cell_dims_for_miller_indices
        self.override_max_q = override_max_q

        if self.override_max_q:
            assert self.all_miller_indices, "Currently overriding maximum q is only supported if searching all Miller indices"
            assert self.max_miller_idx != None, "Require a maximum Miller index as maximum q is overridden." 
        if self.max_miller_idx is None:
            assert self.all_miller_indices == False
        if self.custom_cell_dims_for_miller_indices is not None:
            self.custom_cell_dims_for_miller_indices = np.array(custom_cell_dims_for_miller_indices)/ang_per_bohr
        self.min_q = 0
        if q_minimum!=None:
            self.min_q = q_minimum*ang_per_bohr

        self.hemisphere_screen = True
        eps = 0.00000000000001
        if screen_type == "flat": #well a circle
            self.hemisphere_screen = False
            self.max_q = (1-eps)*self.photon_momentum 
        elif screen_type == "hemisphere":
            # as q = 2ksin(theta), non-inclusive upper bound is 45 degrees(theoretically).  
            self.max_q = (2-eps)*self.photon_momentum/np.sqrt(2)  

        elif screen_type == "sphere":
            self.max_q = (2-eps)*self.photon_momentum
        else:
            raise Exception("Available screen types are 'flat, 'hemisphere', and 'sphere'")
        if q_cutoff != None:
            self.max_q = min(self.max_q,q_cutoff*ang_per_bohr) # convert to atomic units 

        self.t_fineness = t_fineness


        # SPI only... need to refactor
        self.phi_array = np.linspace(0,2*np.pi,self.pixels_per_ring,endpoint=False)  
        self.z_rotation = SPI_z_rotation *np.pi/180
        self.y_rotation = SPI_y_rotation *np.pi/180 # Current rotation of crystal (y axis currently) (did I mean z axis?)
        self.x_rotation = SPI_x_rotation *np.pi/180# Current rotation of crystal (y axis currently)
        self.z_rot_matrix = rotaxis2m(self.z_rotation,Bio_Vect(0, 0, 1))
        self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))     
        self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))   

        # crystal only...
        self.num_orientations =  num_orients_crys
        if orientation_axis_crys != None:
            if np.array(orientation_axis_crys).shape != (3,):
                raise Exception("invalid axis", orientation_axis_crys)            
            if orientation_axis_crys[0] == orientation_axis_crys[1] == orientation_axis_crys[2] == 0:
                raise Exception("axis vector has 0 length")      


        self.orientation_set = None
        if orientation_axis_crys != None:
            axis = Bio_Vect(*orientation_axis_crys)
            ori_set = [rotaxis2m(angle, axis) for angle in np.linspace(0,2*np.pi,self.num_orientations,endpoint=False)]
            self.set_orientation_set([Rotation.from_matrix(m).as_euler("xyz") for m in ori_set])
            print("orientation set set:", self.orientation_set)
    
    def set_orientation_set(self,orientation_set):
        self.orientation_set = orientation_set
        print("orientation set set:", self.orientation_set)
    def get_ff_calculator(self,start_time,end_time,damage_output_handle,parent_dir_path):
        ff_calculator = Plotter(damage_output_handle,parent_dir_path,out_prefix_text = "Calculating form factors...")
        plt.close()
        ff_calculator.initialise_form_factor_params(start_time,end_time,self.max_q,self.photon_energy,t_fineness=self.t_fineness) # q_fineness isn't used for our purposes.   
        return ff_calculator
    
    def spooky_laser(self, start_time, end_time, sim_data_handle, sim_parent_dir_path, target, SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False):
        """ 
        end_time: The end time of the photon capture in femtoseconds. Not a real thing experimentally, but useful for choosing 
        a level of damage. Explicitly, it is used to determine the upper time limit for the integration of the form factor.
        pdb_fpath: The pdb file's path. Changes the variable self.atoms.
        y_orientations: 
            Number of unique y axis rotations to sample crystal. x_axis_rotations not implemented (yet?).
        random_orientation overrides the XFEL class's orientation_set, replacing each with a random orientation. (get same number of orientations though at present TODO.) 
        """
        ff_calculator = self.get_ff_calculator(start_time,end_time,sim_data_handle,sim_parent_dir_path)     
        target.set_ff_calculator(ff_calculator)    
        self.target = target

        if random_orientation == True and self.orientation_set != None:
            raise Exception("Ambiguity: random orientations set to True, but set orientations were provided.")
        if random_orientation == False and self.orientation_set == None:
            raise Exception("Providing an axis of orientation (e.g. format: [0,0,1] for z axis) or enabling random orientations is required (even with a pixel sampling simulation because I coded it bad sorry)") #TODO!
        

        if pixels_across == None and circle_grid == False and SPI:
            raise Exception("Require pixels_across argument for rectangular screen")
        
        # Create output folder for results
        directory = path.abspath(path.join(__file__ ,"../")) + results_parent_dir + self.experiment_name + "/"
        print("creating folder:",directory)
        exist_ok = True
        if (os.path.exists(directory)):
            if clear_output:
                exist_ok = False
                for filename in os.listdir(directory):
                    fpath = os.path.join(directory, filename)
                    if os.path.isfile(fpath): 
                        os.remove(fpath)
                os.rmdir(directory)          
        os.makedirs(directory, exist_ok=exist_ok)      

        
        if SPI and circle_grid:
            ring = np.empty(self.num_rings,dtype="object")
            result = Results_SPI()
            result.I = 0            
            q_sep = (self.max_q-self.min_q)/(self.num_rings)
            q_samples = np.linspace(self.min_q,self.max_q,self.num_rings) + q_sep/2
            for rot_x in range(self.x_orientations):
                self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))      
                self.y_rotation = 0                 
                for rot_y in range(self.y_orientations):              
                    print("Imaging at x, y, rotations:",self.x_rotation,self.y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))      
                    self.y_rotation += 2*np.pi/self.y_orientations              
                    for i, q in enumerate(q_samples):
                        # approximate angle subtended
                        angle_subtended = self.q_to_theta(q+q_sep/2) - self.q_to_theta(q-q_sep/2)
                        ring[i] = self.generate_ring(q,self.phi_array,angle_subtended)
                        #print("q:",q, "x:",ring[i].R,"I[alph=0]",ring[i].I[0])


                    # Initialise stuff that is constant between images (done here to access ring radii.)
                    if rot_y  == 0 and rot_x == 0:
                        phi = self.phi_array
                        #radii = np.zeros(self.num_rings)
                        #for i in range(len(ring)):
                            #radii[i] = ring[i].r           
                        #r, phi = np.meshgrid(radii, phi)     
                        q, phi = np.meshgrid(q_samples,phi)

                    I = np.zeros(q.shape)
                    for ang in range(len(I)):
                        for pos in range(len(I[ang])):
                            I[ang][pos] = ring[pos].I[ang]                                         
                
                    result.I += I/(self.y_orientations*self.x_orientations)
                self.x_rotation += 2*np.pi/self.x_orientations
                result.phi = phi
                result.q = q
            #result.package_up()
            return result
        elif SPI:
            N = pixels_across#100  NxN cells
            cell = np.zeros((N,N),dtype="object")
            result = Results_Grid()      
            
            ### Geometry
            # In crystallography, for a resolution d we have q = 2*pi/d [atomic units] as lambda=2dsin(theta), q = (4pi/lambda)sin(theta). Possibly a questionable definition for SPI without periodicity, but it is used for consistency, and Neutze 2000 makes no distinction.
                  
            # max q represents the actual limit of the change in momentum (at least in a 180 degree arc)
            # rim_q is the q we want to have at the largest unbroken ring
            rim_q = self.max_q
            if SPI_resolution!= None:
                d = SPI_resolution/ang_per_bohr # resolution
                rim_q = res_to_q(d)
                if self.max_q < rim_q:
                    print ("WARNING: resolution of " + str(SPI_resolution) + " angstroms requires q to go beyond its maximum. Using max_q="+str(self.max_q/ang_per_bohr)+" instead.")
                    rim_q = self.max_q
            print("rim q:",rim_q)

            max_theta = self.q_to_theta(rim_q) # maximum theta of full ring.
            #screen_width = resolution_to_X(d) * 2
            screen_distance = self.detector_distance # screen-target separation [a0]
            
            # X is the distance from the centre of the screen to the point of incidence (flat screen)
            # We know where the cells are, they are equally spaced. We also know the distance from the screen to the detector.
            # This gives us theta. 
            def X_to_theta(x):
                '''
                Assumes screen distance in same units as x (a0)
                '''
                return 0.5*np.arctan2(x,screen_distance)
            def theta_to_X(theta):
                '''
                Assumes screen distance in same units as x (a0)
                '''
                return np.abs(screen_distance*np.tan(2*theta))
            
            # We must determine q at each point for finding I
            def X_to_q(x):
                '''
                Returns q [1/a0]
                '''
                theta = X_to_theta(x)
                k0 = self.photon_momentum #2*np.pi/E_to_lamb(photon_energy)
                return 2*k0*np.sin(theta)   
            
            # res. at edge corner or centre? Surely at centre, for a full ring of information
            screen_width = 2*theta_to_X(max_theta)  # edge centre  # We have placed our screen to have the desired resolution at the "rim".
            #screen_width=  2 * (np.sin(2*max_theta)/np.sqrt(2)) * screen_distance # corner
            print("Screen width:",round(screen_width*ang_per_bohr/1e7,2),"mm")            
            # Trig consistency check
            assert round(X_to_q(screen_width/2),2) == round(rim_q,2)

            ## Calculate q for each cell.
            corner_q = X_to_q(np.sqrt(2*(screen_width/2)**2))
            result.cell_width = screen_width/len(cell)
            q_grid = np.empty(cell.shape)
            result.xy = np.empty(cell.shape+(2,))
            result.I = np.zeros(cell.shape) 
            result.zero_angle_I = 0 # In case dim of axis has even number of pixels           
            for i, x in enumerate(cell):
                for j, y in enumerate(x):
                    x = result.cell_width*(i-(len(cell)-1)/2)
                    y = result.cell_width*(j-(len(cell)-1)/2)
                    if self.hemisphere_screen:
                        print ("hemisphere not supported for SPI yet")
                    
                    dist = np.sqrt(x**2+y**2)
                    q_grid[i,j] = X_to_q(dist)
                    result.xy[i,j] = np.array([x,y])
            #result.q_scr_xy = X_to_q_scr(result.xy)
            result.q_xy = X_to_q(result.xy)
            # Store a mask for the values that we should ignore.
            mask = copy.deepcopy(q_grid)
            mask[mask < self.min_q] = 0; mask[mask > self.max_q] = 0
            mask[mask != 0] = 1

            #print("MASK",mask)
            #result.I[mask == 0] = None
            result.full_ring_mask = mask
            #result.I*=mask
            # Calculate I for each cell from q (need to vectorise q)
            for rot_x in range(self.x_orientations):
                self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))      
                #self.y_rotation = 0                 
                for rot_y in range(self.y_orientations):                  
                    print("Imaging at x, y, rotations:",self.x_rotation,self.y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))      
                    self.y_rotation += 2*np.pi/self.y_orientations              
                    result.I += self.generate_cell_intensity(q_grid,result.xy)
                    result.zero_angle_I += self.generate_cell_intensity(np.array([[0]]),np.array([[[0,0]]]))

                self.x_rotation += 2*np.pi/self.x_orientations             
            # convert to angstrom
            result.xy *= ang_per_bohr
            result.cell_width *= ang_per_bohr
            result.q_xy /= ang_per_bohr
            # Average out intensity # NOTE intensities aren't aligned. Shouldn't be using R factor directly on result from multiple orientations for SPI.
            # TODO Intensity should really be a list of results.
            result.I /= (self.y_orientations*self.x_orientations)
            #result.package_up()
            return result        

        else:
            # Iterate through each orientation of crystal, picklin' up a file for each orientation
            used_orientations = []
            if random_orientation:
                self.orientation_set = [(0,0,0)]*self.num_orientations  # Dummy orientations
            for j, cardan_angles in enumerate(self.orientation_set):             
                print("Imaging orientation",j)
                bragg_points, miller_indices,cardan_angles = self.bragg_points(target,cell_packing = target.cell_packing, cardan_angles = cardan_angles,random_orientation=random_orientation)
                used_orientations.append(cardan_angles)
                num_points = int(len(bragg_points))
                result = Results(num_points,j)
                # Get the q vectors where non-zero
                i = 0
                #TODO vectorise
                # (Assume pixel adjacent to bragg point does not capture remnants of sinc function)
                point = self.generate_point(bragg_points,cardan_angles)
                # fix this filling stuff
                result.phi = point.phi
                result.phi_aligned = point.phi_crystal_aligned
                result.X = point.X
                result.q = point.q
                result.I += point.I
                
                result.package_up(miller_indices)
                #Save the result object into its own file within the output folder for the experiment
                fpath = directory + str(cardan_angles) +".pickle"
                with open(fpath,"wb") as pickle_out:
                    pickle.dump(result,pickle_out)
            return used_orientations

    class Feature:
        def __init__(self,q,X,theta):
            self.q = q
            self.X = X # radial distance from centre of screen
            self.theta = theta          
    class Ring(Feature):
        def __init__(self,*args):
            super().__init__(*args)
    class Spot(Feature):
        def __init__(self,*args):
            super().__init__(*args)
    class Cell(Feature):
        def __init__(self,*args):
            super().__init__(*args)
            
    
    def generate_cell_intensity(self,q,xy_grid):
        solid_angle = 1  # approximate all as same for now.
        theta = self.q_to_theta(q)
        X = self.q_to_X(q)     
        if self.hemisphere_screen:
            print("hemisphere/spherical screen not supported for SPI yet.")           
        cell = self.Cell(q,X,theta)
        cell.q_parr_screen = self.q_to_q_scr(q)
        phis = np.arctan2(xy_grid[:,:,1],xy_grid[:,:,0])
        #print("phis",phis)
        return self.illuminate(cell,phis=phis,SPI_proj_solid_angle = solid_angle)

    def generate_ring(self,q,phi_array,angle_subtended):
        '''Returns the intensity(phi) array and the radius for given q.'''
        #print("q=",q)
        #r = self.q_to_r(q)# TODO
        theta = self.q_to_theta(q)
        X = self.q_to_X(q)
        solid_angle = (2*np.pi/len(phi_array)) * (np.pi/angle_subtended)   # solid angle in fractions [sp]. Solved for sphere but will be same when project to flat detector. 
        proj_solid_angle = solid_angle  # laser source.
        if self.hemisphere_screen:
            print("hemisphere/spherical screen not supported for SPI yet.")
        ring = self.Ring(q,X,theta)
        ring.I = self.illuminate(ring,phis=phi_array,SPI_proj_solid_angle=proj_solid_angle)
        return ring 
    def generate_point(self,G,cardan_angles): # G = vector
        ''' We store variables like G since they correspond to the special bragg points, unlike points from SPI that are just samples of the continuous pattern. 
        
        '''
        if len(G.shape) > 1:
            G = np.moveaxis(G,0,len(G.shape)-1)  # moves the axis corresponding to the individual momentum to the back, giving us dim = [3,num_G].
        q = np.sqrt(np.power(G[0],2)+np.power(G[1],2)+np.power(G[2],2))
        #r = self.q_to_r(q) # TODO
        X = self.q_to_X(q)
        if self.hemisphere_screen:
            X = self.q_to_q_scr_curved(G)/self.photon_momentum*self.detector_distance  # not tested...
        theta = self.q_to_theta(q)
        point = self.Spot(q,X,theta)
        point.phi = np.arctan2(G[1],G[0])
        point.G = G
        #Crystal-aligned phi. i.e. we realign all the images to be in the 0 0 0 orientation.
        G = self.rotate_G_to_orientation(G,*cardan_angles,inverse=True)[0]
        point.phi_crystal_aligned = np.arctan2(G[1],G[0])
                
        if DEBUG:
            print("X",X,"phi",point.phi,"G",G[0],G[1],G[2])


        #print("phi",point.phi,"q_parr_screen",point.q_parr_screen)
        
        point.I = self.illuminate(point,cardan_angles=cardan_angles)
        # # Trig check      
        # check = smth  # check = np.sqrt(G[0]**2+G[1]**2)
        # if q_parr_screen != check:
        #     print("Error, q_parr_screen =",q_parr_screen,"but expected",check)  
        return point   
    

    # Returns the relative intensity at point q for the target's unit cell, i.e. ignoring crystalline effects.
    # If the feature is a bragg spot, this gives its relative intensity, but due to photon conservation won't be the same as the intensity without crystallinity - additionally different factors for non-zero form factors occur across different crystal patterns.
    def illuminate(self,feature, phis = None,cardan_angles = None, SPI_proj_solid_angle=None):  # Feature = ring or spot.
        """Returns the intensity at q. Not crystalline yet.
        phis is an ndarray containing the angle of each point to calculate for (SPI) rings, but should contain only 1 element for points.
        """
        # if phis == None:
        #     phis = self.phi_array
        
        if type(feature) == self.Spot:
            SPI = False
        elif type(feature) in [self.Cell,self.Ring]:
            SPI = True

        F_shape = tuple()
        if type(feature) is self.Spot:
            F_shape += (self.t_fineness,)           
            if len(feature.G.shape) > 1:
                F_shape += (len(feature.G[2]),) #[times,num_G]   
        else:
            if type(feature) is self.Ring:
                F_shape += phis.shape  
            F_shape += (self.t_fineness,)# [?phis?,times]
            if type(feature.q) is np.ndarray:
                F_shape += feature.q.shape          # [?phis?,times,feature.q.shape]
        F_supercells = np.zeros(self.target.supercell_simulations,dtype="object")
        for S in range(self.target.supercell_simulations):   
            print("Simulating supercell", S)
            times_used = None
            for species in self.target.species_dict.values():
                species.set_stochastic_states()   
                species.set_coord_deviation()             
                species.set_scalar_form_factor()
                times_used = species.times_used
            if (times_used[-1] == times_used[0]):   
                raise Exception("Intensity array's final time equals its initial time")                  
            # Technically sum of F(t)*sqrt(J(t)), where F = sum(f(q,t)*T(q)), and J(t) is the incident intensity, thus accounting for the pulse profile. (J(t) is accounted for in get_stochastic_f)
            F_sum = np.zeros(F_shape,dtype="complex_")  
            for species in self.target.species_dict.values():
                print("------------------------------------------------------------")
                print("Getting contribution to integrand from species",species.name)
                if not np.array_equal(times_used,species.times_used):
                    raise Exception("Times used don't match between species.")        
                # iterate through every atom including in each symmetry of unit cell (each asymmetric unit)
                max_atoms_per_loop = 1000 # Restrict array size to prevent computer explosions. 
                for s in range(len(self.target.sym_rotations)):
                    #print("Working through symmetry",s)
                    num_atom_batches = int(len(species.coords)/max_atoms_per_loop)+1
                    for a_batch in range(num_atom_batches):
                        atm_idx = np.arange(len(species.coords)*s+max_atoms_per_loop*a_batch, len(species.coords)*s + min(len(species.coords), max_atoms_per_loop*(a_batch+1)))
                        if len(atm_idx) == 0:
                            break
                        if self.target.supercell_simulations < 10:
                            if len(self.target.sym_rotations) < 50 or (len(self.target.sym_rotations) < 2000 and s%100 == 0)  or s%1000 == 0:
                                print("symmetry:",s,"atoms:",species.name,atm_idx[0],"-",atm_idx[-1])
                        relative_atm_idx = np.arange(max_atoms_per_loop*a_batch, min(len(species.coords), max_atoms_per_loop*(a_batch+1)))
                        # Rotate to target's current orientation 
                        # rot matrices are from bio python and are LEFT multiplying. TODO should be consistent replace this with right mult. 
                        R = np.array(species.coords[relative_atm_idx[0]:relative_atm_idx[-1]+1]) 
                        coord = self.target.get_sym_xfmed_point(R,s)  + species.error[relative_atm_idx] # dim = [ N, 3], where N is number of coords.
                        coord = coord @ self.x_rot_matrix  
                        coord = coord @ self.y_rot_matrix
                        coord = coord @ self.z_rot_matrix
                        # Get spatial factor T
                        if SPI:
                            T = self.SPI_interference_factor(phis,coord,feature)  # if grid: [phis,qx,qy]  if ring: [phis,q] (unimplemented)
                        else:
                            T= self.interference_factor(coord,feature,cardan_angles)  #[num_G] 
                        f = species.get_stochastic_f(atm_idx, feature.q)  / np.sqrt(self.target.num_cells) # Dividing by np.sqrt(self.num_cells) so that fluence is same regardless of num cells. 
                        same_each_sym = False # (debug)
                        if same_each_sym:
                            f = species.get_stochastic_f(relative_atm_idx, feature.q)  / np.sqrt(self.target.num_cells) # Dividing by np.sqrt(self.num_cells) so that fluence is same regardless of num cells. 
                        #print(F_sum.shape,T.shape,f.shape) 
                        if SPI: 
                            if type(feature) is self.Cell:
                                F_sum += np.sum(T[:,None,...]*f,axis=0)          #[num_atoms,None,qX,qY] X [num_atoms,times,qX,qY] -> [times,qX,qY]
                            if type(feature) is self.Ring:           
                                F_sum += np.sum(T[:,:,None,:] * f[:,None,:,:],axis=0)        # [num_atoms,phis,None,num_|q|]X[num_atoms,None,times, num_|q|]  -> [phis,times,num_|q|]
                        else:
                            F_sum += np.sum(T[:,None,:] * f,axis=0)                           # [num_atoms,None,num_G]X[num_atoms,times,num_G]  ->[times,num_G] 
                        #print("s,F_sum",s,F_sum)
                                                                            #I =  np.square(np.abs(F_sum[:,0]))  # 
            F_supercells[S] = F_sum
        # Add supercells
        F_cry = np.zeros(F_shape,dtype="complex_")
        length = np.ceil(self.target.num_supercells**(1/3)) 
        super_batch_size = 50
        supercells_remaining = self.target.num_supercells
            
        x, y, z= np.meshgrid(np.arange(0, length), np.arange(0, length), np.arange(0, length))
        super_cube_coords = np.stack([ y.flatten(), x.flatten(), z.flatten()], axis = -1) # sadly not dimensionally-transcendental enough to be prefixed "hyper"        
        curr_cube_idx = 0
        while supercells_remaining > 0:
            end_cube_idx = min(supercells_remaining,super_batch_size)
            super_coords = np.empty((self.target.num_supercells,3))  
            super_coords = super_cube_coords[curr_cube_idx:end_cube_idx]*self.target.supercell_dim
            F_supercell_copies = np.zeros(shape=(len(super_coords),)+(F_supercells[0].shape),dtype=complex)
            for p in range(len(F_supercell_copies)):
                F_supercell_copies[p] = np.random.choice(F_supercells)
            if SPI:
                T_supercell = self.SPI_interference_factor(phis,super_coords,feature)
            else:
                T_supercell = self.interference_factor(super_coords,feature,cardan_angles)
            F_cry += np.sum(F_supercell_copies*T_supercell[:,None],axis=0)
            supercells_remaining -= super_batch_size
            curr_cube_idx = end_cube_idx
        # Integrate over time to get the intensity
        time_axis = 0
        if type(feature) is self.Ring:
            time_axis = 1
        I = np.trapz(np.square(np.abs(F_cry)),times_used,axis = time_axis) / (times_used[-1]-times_used[0])       #[num_G] for points, or for SPI: [phis,feature.q.shape], corresponding to rings or square grid
            

        # For unpolarised light (as used by Neutze 2000). (1/2)r_e^2(1+cos^2(theta)) is thomson scattering - recovering the correct equation for a lone electron, where |f|^2 = 1 by definition.    
        # Generally not important due to small angles involved.
        thet = self.q_to_theta(feature.q)
        r_e_sqr = 2.83570628e-9
        I*= r_e_sqr*(1/2)*(1+np.square(np.cos(2*thet))) 

        if SPI:
            # For crystal we can approximate infinitely small pixels and just consider bragg points.
            # But for SPI need to take the pixel's size into account. Neutze 2000 makes the following approximation:
            I *=  SPI_proj_solid_angle    # equiv. to *= solid_angle
        
        I/=10**12  # Scale down intensity.
        print("Total screen-incident intensity = ","{:e}".format(np.sum(I)))
        return I

    # def illuminate_average(self,feature,phis = None,cardan_angles = None):  # Feature = ring or spot.
    #     """Returns the intensity at q. Not crystalline yet."""
    #     if phis == None:
    #         phis = self.phi_array
    #     F = np.zeros(phis.shape,dtype="complex_")
    #     for species in self.target.species_dict.values():
    #         species.set_scalar_form_factor(feature.q)
    #         # iterate through each symmetry of unit cell (each asymmetric unit)
    #         for s in range(len(self.target.sym_rotations)):
    #             for R in species.coords:
    #                     # Rotate to crystal's current orientation 
    #                     R = R.left_multiply(self.y_rot_matrix)  
    #                     R = R.left_multiply(self.x_rot_matrix)   
    #                     # PDB format note: if x axis has dim of X, we have a point between [- X, X].
    #                     coord = np.multiply(R.get_array(),self.target.sym_rotations[s]) + np.multiply(self.target.supercell_dim,self.target.sym_translations[s])
    #                     # Get spatial factor T
    #                     T = np.zeros(phis.shape,dtype="complex_")
    #                     if SPI:
    #                         T = self.SPI_interference_factor(phis,coord,feature)
    #                     else:
    #                         T= self.interference_factor(coord,feature,cardan_angles)
    #                     F += species.ff_average*T   
    #                     # Rotate atom for next sample            
    #     I = np.square(np.abs(F))        
        

    #     return I 

    def interference_factor(self,coord,feature,cardan_angles):
        """ theta = scattering angle relative to z-y plane """ 
        # Rotate our G vector BACK to the real laser orientation relative to the crystal.
        q_vect = self.rotate_G_to_orientation(feature.G.copy(),*cardan_angles,inverse=True)[0]       
        coord = np.moveaxis(coord,0,-1)  #  dim = [xyz,atoms]
        q_vect = np.moveaxis(q_vect,0,-1) # dim = [momenta,xyz]
        q_dot_r = np.apply_along_axis(np.dot,len(q_vect.shape)-1,q_vect,coord) # dim = [num_G]  
        q_dot_r = np.moveaxis(q_dot_r,-1,0)
        coord = np.moveaxis(coord,-1,0)                            
        T = np.exp(-1j*q_dot_r)
        return T

    def SPI_interference_factor(self,phi_array,coord,feature):  #TODO refactor SPI features so can just use above (allows for realignment)
        """ theta = scattering angle relative to z-y plane 
        
        """ 
        q_z = np.multiply(feature.q,np.sin(feature.theta))                          # dim = [qX,qY] for square screen with cell coords given by X,Y 
        # phi is angle of vector relative to x-y plane, pointing from screen centre to point of incidence. 
        # TODO need to get working with rings again since they were designed to have multi phi for each q. Here it is 1:1
        if type(feature) is not self.Cell:
            raise Exception("accidentally removed ring implementation sorry!")
        q_z = q_z       
        q_y = feature.q_parr_screen * np.sin(phi_array)
        q_x = feature.q_parr_screen * np.cos(phi_array)   # dim = [qX,qY]
        #print("====="); print(phi_array);print("-----");print(q_z.shape,q_y.shape,q_x.shape)   
        q_vect = np.array([q_x,q_y,q_z])
        q_vect = np.moveaxis(q_vect,0,len(q_vect.shape)-1)    # dim = [qX,qX,3]
        coord = np.moveaxis(coord,0,-1)
        q_dot_r = np.apply_along_axis(np.matmul,len(q_vect.shape)-1,q_vect,coord) # dim = [phis,qX,qY]       q_dot_r = np.tensordot(q_vect,coord,axes=len(q_vect.shape)-1)
        q_dot_r = np.moveaxis(q_dot_r,-1,0)
        coord = np.moveaxis(coord,-1,0)
        if (type(feature.q) == np.ndarray): 
            if len(coord.shape) == 2 and q_dot_r.shape != (coord.shape[0],) + feature.q.shape or len(coord.shape) == 1 and q_dot_r.shape != feature.q.shape:
                raise Exception("Unexpected q_dot_r shape.",q_dot_r.shape)
        T = np.exp(-1j*q_dot_r) 
        return T   
    def bragg_points(self,crystal, cell_packing, cardan_angles,random_orientation=False):
        ''' 
        Using the unit cell structure, find non-zero values of q for which bragg 
        points appear.
        lattice_vectors e.g. = [a,b,c] - length of each spatial vector in orthogonal basis.
        Currently only checked to work with cubics. 
        '''

        def get_G(miller_indices,cartesian=True):
            '''
            Get G in global cartesian coordinates (crystal/zone axes not implemented),
            where G = G[0]x + G[1]y + G[2]z
            '''       #TODO vectorise somehow      
            G = np.zeros(miller_indices.shape)
            if cartesian:
                for i in range(len(miller_indices)):
                    G[i] = np.array(np.dot(miller_indices[i],b))
                    if (i == 2 or i == 10) and DEBUG:
                        print("G Debug")
                        print(b)
                        print(miller_indices[i],G[i])
                #G = G.reshape(3,)
                G, used_angles = self.rotate_G_to_orientation(G,*cardan_angles,random = random_orientation,G_idx_first = True) 
            else: 
                G = "non-cartesian Unimplemented"
            return G, used_angles
        
        # (h,k,l) is G (subject to selection condition) in lattice vector basis (lattice vector length = 1 in each dimension):
        # a = lattice vectors, b = reciprocal lattice vector
        if cell_packing == "SC" or cell_packing == "FCC" or cell_packing == "BCC" or cell_packing == "FCC-D":
            if cell_packing == "SC":
                a = np.array([[1,0,0],[0,1,0],[0,0,1]])
            if cell_packing == "BCC":
                a = 0.5*np.array([[-1,1,1],[1,-1,1],[1,1,-1]])
            if cell_packing == "FCC":
                a = 0.5*np.array([[0,1,1],[1,0,1],[1,1,0]])
            if self.custom_cell_dims_for_miller_indices is None:
                a = np.multiply(a,crystal.cell_dim)
            else:
                a = np.multiply(a,self.custom_cell_dims_for_miller_indices) 
        else: 
            raise Exception("Unknown cell packing type")
        if not np.array_equal(crystal.supercell_dim,crystal.cell_dim):
            print("Warning: This code samples the intensity at the miller indices, and thus does not capture peak broadening - consider using the SPI imaging mode instead.")
        b1 = np.cross(a[1],a[2])
        b2 = np.cross(a[2],a[0])
        b3 = np.cross(a[0],a[1])
        b = 2*np.pi*np.array([b1,b2,b3])/(np.dot(a[0],np.cross(a[1],a[2])))
        
        # Cast a wide net, catching all possible permutations of miller indices.
        #   G = hb1 + kb2 + lb3. (bi = lattice vector)
        #   !Attention! Assuming vectors are orthogonal.
        # TODO double check not cutting off possible values.
        if not self.override_max_q:
            q_1_max = self.max_q
            q_2_max = self.max_q
            q_3_max = self.max_q        # q = (0,0,l) case.
            h_max = 0
            while np.sqrt(sum(pow(h_max*element, 2) for element in b[0])) <= q_1_max:
                h_max += 1
            k_max = 0
            while np.sqrt(sum(pow(k_max*element, 2) for element in b[1])) <= q_2_max:
                k_max += 1 
            l_max = 0
            while np.sqrt(sum(pow(l_max*element, 2) for element in b[2])) <= q_3_max:
                l_max += 1             
        else:
            h_max = k_max = l_max = self.max_miller_idx     

        h_set = np.arange(-h_max,h_max+1,1)
        k_set = np.arange(-k_max,k_max+1,1)
        l_set = np.arange(-l_max,l_max+1,1)
        indices = list(itertools.product(h_set,k_set,l_set))
        indices = np.array([*set(indices)])   
    
        # Selection rules
        if cell_packing == "SC":
            selection_rule = lambda f: True
        if cell_packing == "BCC":
            selection_rule = lambda f: (f[0]+f[1]+f[2])%2==0  # All even
        if cell_packing == "FCC":
            selection_rule = lambda f: (np.abs(f[0])%2+np.abs(f[1])%2+np.abs(f[2])%2) in [0,3]   # All odd or all even.
        if cell_packing == "FCC-D":
            selection_rule = lambda f: (np.abs(f[0])%2+np.abs(f[1])%2+np.abs(f[2])%2) == 3 or ((np.abs(f[0])%2+np.abs(f[1])%2+np.abs(f[2])%2) == 0 and (f[0] + f[1] + f[2])%4 == 0)  # All odd or all even.    

        # Set G, and retrieve the cardan angles used (in case of random orientations)
        G_temp, cardan_angles = get_G(indices) 
        # If we used a random orientation, we now lock in the orientations just generated and contained in cardan_angles.
        random_orientation = False 

        def miller_selection_rule():
            return None
        if self.max_miller_idx != None:
            m = self.max_miller_idx
            max_g_vect = get_G(np.full((1,3),m))[0][0]
            if not self.override_max_q:
                self.max_q = min(self.max_q,np.sqrt(((max_g_vect[0])**2+(max_g_vect[1])**2+(max_g_vect[2])**2)))
            else:
                self.max_q = np.sqrt(((max_g_vect[0])**2+(max_g_vect[1])**2+(max_g_vect[2])**2))
            def miller_selection_rule(indices): # Probably unnecessary I'm just making sure...
                return  (abs(indices[0]) <= m and abs(indices[1]) <= m and abs(indices[2]) <= m)
        print("max q (i.e. rim q):",self.max_q)
        
        print("using q range of ", self.min_q/ang_per_bohr,"-",self.max_q/ang_per_bohr," angstrom-1")
        def min_max_q_rule(g):
            return self.min_q <= np.sqrt(((g[0])**2+(g[1])**2+(g[2])**2)) <= self.max_q
        #max_q_rule = lambda f: np.sqrt(((f[0]*np.average(cell_dim))**2+(f[1]*np.average(cell_dim))**2+(f[2]*np.average(cell_dim))**2))<= self.max_q
        
        # Catch the miller indices with a boolean mask
        if self.all_miller_indices:
            mask=np.apply_along_axis(miller_selection_rule,1,indices)*np.apply_along_axis(selection_rule,1,indices)
            if not self.override_max_q:
                mask*=np.apply_along_axis(min_max_q_rule,1,G_temp)
        else:
            mask = np.apply_along_axis(selection_rule,1,indices)*np.apply_along_axis(min_max_q_rule,1,G_temp)*np.apply_along_axis(self.mosaic_elastic_condition,1,G_temp)
            if self.max_miller_idx != None:
                mask*=np.apply_along_axis(miller_selection_rule,1,indices)
        indices = indices[mask]

        print("Cardan angles:",cardan_angles)
        print("Number of points:", len(indices))   
        assert len(indices) < 10920, "Output size failsafe triggered, output size likely >~ 1 GiB."
        if len(indices) > 2000:
            print("WARNING: very high number of points!")        
        for elem in indices:
            if DEBUG:
                print(elem)
        

        G, cardan_angles = get_G(indices,cardan_angles)
        if DEBUG:
            print(G)
        return G, indices, cardan_angles

    def rotate_G_to_orientation(self,G,alpha=0,beta=0,gamma=0,random=False,inverse = False, G_idx_first = False):
        '''
        Assumes G.shape = (num_momenta,3)
        Rotates the G vectors to correspond to the crystal when oriented to the given cardan angles. 
        alpha: x-axis angle [radians].
        beta:  y-axis angle [radians].
        gamma: z-axis angle [radians].
        ## Returns:
        G: transformed G
        list: [alpha,beta,gamma] (the angles used)
        '''
        moved_axis = False
        if G_idx_first and G.shape[len(G.shape)-1] != 3:
            raise Exception("please give array of form [...,num_G,3]")
        elif not G_idx_first and G.shape[len(G.shape) -2] !=3:
            raise Exception("please give array of form [...,3,num_G]")
            
        if G_idx_first:
            G = np.swapaxes(G,len(G.shape)-2,len(G.shape)-1)

        if type(random) != bool:
            raise Exception("random_orientation must be boolean")
        if random == True:
            alpha = np.random.random_sample()*2*np.pi
            beta = np.random.random_sample()*2*np.pi
            gamma = np.random.random_sample()*2*np.pi

        R = Rotation.from_euler('xyz',[alpha,beta,gamma])
        rot_matrix = R.as_matrix()
        if inverse: 
            rot_matrix = rot_matrix.T # transpose == inverse
        # Apply rotation
        #print("Rotation matrix:")
        #print(R.as_matrix())
        G = np.around(rot_matrix,decimals=10) @ G     # [3,3]X[3,num_G]

        if G_idx_first:
            G = np.swapaxes(G,len(G.shape)-2,len(G.shape)-1)  
        return G, [alpha,beta,gamma]  # [num_G,3]
                
    def q_to_theta(self,q):
        '''
        q momentum transfer [1/bohr]
        '''
        lamb = E_to_lamb(self.photon_energy)
        theta = np.arcsin(lamb*q/(4*np.pi)) 
        #if not (theta >= 0 and theta < 45):
            #print("warning, theta out of bounds")        
        return theta


    # def q_to_r(self, q):
    #     """Returns the distance from centre of screen that is struck by photon which transferred q (magnitude)"""
    #     D = self.detector_distance #bohr
    #     theta = self.q_to_theta(q)
    #     return D*np.tan(2*theta)

    #TODO I'll have to check later if this is crazy or not
    '''   
    def q_to_scr_pos(self,q):
        """ 
        """
        theta = self.q_to_theta(q)
        # From here it will help to consider q as referring to the final momentum. /// what is this comment??????? 
        v_x = c_au*np.cos(2*theta)
        v_y = c_au*np.sin(2*theta)
        D = self.detector_distance
        t = D/v_x
        radius = v_y*t
        return radius
    '''
    def q_to_q_scr(self,q):
            """ Returns the screen-parallel component of q - useful for getting q_x and q_y
            """
            theta = self.q_to_theta(q)
            return q/(np.sin(np.pi/2-theta))      
    def q_to_X(self,q):
        '''
        Assumes screen distance in same units as x (a0)
        '''
        theta = self.q_to_theta(q)
        return np.abs(self.detector_distance*np.tan(2*theta))    
               
    def q_to_q_scr_curved(self,G):
        return np.sqrt(G[0]**2+G[1]**2)
    
    def mosaic_elastic_condition(self,q_vect):
        ''' determines whether vector q is allowed.'''
        # Rocking angle version
        q = np.sqrt(q_vect[0]**2 + q_vect[1]**2 + q_vect[2]**2)
        if q > self.max_q:
            return False
        theta = self.q_to_theta(q)
        k = self.photon_momentum
        #Assuming theta > 0:
        delt = self.target.rocking_angle
        min_err = k*sin(theta - delt/2)
        max_err = k*sin(theta+delt/2)
        err = np.sqrt(q_vect[0]**2 + q_vect[1]**2 + (q_vect[2]-k)**2) - k
        if min_err <= err <= max_err:
            return True
        return False
        

    def r_to_q(self,r):
        D = self.detector_distance
        lamb = E_to_lamb(self.photon_energy) 
        theta = np.arctan(r/D)/2
        q = 4*np.pi/lamb*np.sin(theta)
        return q 
    # def r_to_q_scr(self,r):
    #     q = self.r_to_q(r)
    #     return self.q_to_q_scr(q)

def E_to_lamb(photon_energy):
    """Energy (eV) to wavelength in A.U."""
    E = photon_energy  # eV
    return 2*np.pi*c_au/(E/eV_per_Ha)
        #
        # q = 4*pi*sin(theta)/lambda = 2pi*u, where q is the momentum in AU (a_0^-1), u is the spatial frequency.
        # Bragg's law: n*lambda = 2*d*sin(theta). d = gap between atoms n layers apart i.e. a measure of theoretical resolution.

def scatter_scatter_plot(get_R_only = False,neutze_R = True, crystal_aligned_frame = False ,SPI_result1 = None, SPI_result2 = None, full_range = True,num_arcs = 50,num_subdivisions = 40, result_handle = None, results_parent_dir = RESULTS_LOCAL_PATH, compare_handle = None, normalise_intensity_map = False, show_grid = False, cmap_power = 1, cmap = None, min_alpha = 0.05, max_alpha = 1, 
                         bg_colour = "grey",solid_colour = "white", show_labels = False, radial_lim = None, plot_against_q=False,log_I = True, log_dot = False,  fixed_dot_size = False, dot_size = 1, crystal_pattern_only = False, log_radial=False,cutoff_log_intensity = None,spi_full_rings_only=True,min_R_dmg_pixel = 0.1,cmap2=None,log_range=None,custom_fig_width=None,custom_fig_height=None,log_diff_vmin = -0.5, log_diff_vmax = 0.5, normalize_to_centre = True, dpi=100):
    ''' (Complete spaghetti at this point.)
    Plots the simulated scattering image.
    result_handle:

    compare_handle:
        results/compare_handle/ is the directory of results that will be subtracted from those in the results directory.
        Must have same orientations.
    
    '''
    # Pixel plots...
    if custom_fig_height is None:
        custom_fig_height = plt.rcParams['figure.figsize'][1]
    if custom_fig_width is None:
        custom_fig_width = plt.rcParams['figure.figsize'][0]
    LOG10 = False
    if LOG10:
        log_function = np.log10
        if log_range is None:
            log_range = 5
    else:
        log_function = np.log
        if log_range is None:
            log_range = 10
    print("=====================Plotting===========================")
    mpl.rcParams['text.color'] = "blue"
    mpl.rcParams['axes.labelcolor'] = "blue"
    mpl.rcParams['xtick.color'] = "black"
    mpl.rcParams['ytick.color'] = "blue"    
    #https://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
    if cmap != None:
        class MplColorHelper:

            def __init__(self, cmap_name, start_val, stop_val):
                self.cmap_name = cmap_name
                self.cmap = plt.get_cmap(cmap_name)
                self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
                self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

            def get_rgb(self, val):
                bad_val = False
                mod_val = val
                if compare_handle == None:
                    # Plot actual image!
                    for elem in (val,min_z,max_z,cmap_power):
                        if elem == -np.inf or elem == None or elem == np.inf:
                            bad_val = True   
                    if val < min_z or val > max_z:
                        bad_val = True
                    if bad_val:
                        raise Exception("invalid value within",val,min_z,max_z,cmap_power)
                    mod_val = ((val-min_z)/(max_z-min_z))**cmap_power
                
                if mod_val < 0 or mod_val > 1:
                    raise Exception("rgb val outside acceptable range")
                return self.scalarMap.to_rgba(mod_val)            
    else:
        # Broken TODO
        r,g,b = to_rgb(solid_colour)
        colours = [(r,g,b,a) for a in np.clip(z/max_z,min_alpha,max_alpha)]
    if cmap2 is None:
        cmap2 = cmap
    def add_screen_properties(fig_width=14,fig_height=8.4):
        if radial_lim:
            bottom,top = plt.ylim()
            plt.ylim(bottom,radial_lim)
        if log_radial:
            plt.yscale("log")  
        plt.gca().set_facecolor(bg_colour)      
        #plt.gcf().set_figwidth(20)         # if not using widget magic
        #plt.gcf().set_figheight(20)        #       
            
        plt.gcf().set_figwidth(fig_width)
        plt.gcf().set_figheight(fig_height)                  

    if result_handle != None: #TODO replace this atrocious way of distinguishing between inf. crystal and finite
        if plot_against_q:
            radial_lim /= ang_per_bohr
        else:
            radial_lim*= ang_per_bohr        
        
        results_dir = results_parent_dir+result_handle+"/"
        compare_dir = None
        if compare_handle!= None:
            compare_dir = results_parent_dir+compare_handle+"/"        

        ## Point-like (Crystalline)
        # Initialise R factor sector comparison plot. 
        sector_histogram = np.zeros((num_arcs,num_subdivisions)).T
        sector_num_histogram = np.zeros(sector_histogram.shape)
        sector_den_histogram = np.zeros(sector_histogram.shape)
        R_histogram = np.zeros((num_arcs,num_subdivisions)).T
        R_num_histogram = np.zeros(R_histogram.shape)
        R_den_histogram = np.zeros(R_histogram.shape)
        # Bin edges
        phi_edges = np.linspace(-np.pi,np.pi,num_arcs+1)
        radial_edges = np.linspace(0,radial_lim,num_subdivisions+1)        
        sector_phi,sector_radial = np.meshgrid(phi_edges,radial_edges)
        num_orientations = 0
        def get_old_histogram_contribution(I,phi,radial_axis):
                # We don't set density = True, because we don't want to normalise the weights.
                num_samples = np.histogram2d(phi, radial_axis, bins=(phi_edges, radial_edges))[0]
                num_samples[num_samples == 0] = 1 # avoid division by zero
                H = np.histogram2d(phi, radial_axis, weights=I, bins=(phi_edges, radial_edges))[0]
                H = H.T
                H = np.divide(H, num_samples.T)
                return H/num_orientations   
        
        def get_sector_histogram_contribution(I_ideal,I_real,phi,radial_axis,phi_edges = phi_edges):
            '''I_ideal = Intensities of undamaged target
               I_real = Intensities of damaged target
            '''#                                                                            _|-|_
            # R = Σ|sqrt(I_ideal) - sqrt(I_real|) / (Σ sqrt(I_ideal))  |   me -> ( ._.)ヽ(￣┏＿┓￣ R)  "thousands of lines of code just for you Señor R factor".
            N = np.abs(np.sqrt(I_ideal) - np.sqrt(I_real))
            D = np.sqrt(I_ideal)
            numerators = np.histogram2d(phi, radial_axis, weights=N, bins=(phi_edges, radial_edges))[0]
            denominators  = np.histogram2d(phi, radial_axis, weights=D, bins=(phi_edges, radial_edges))[0]
            numerators = numerators.T; denominators = denominators.T
            return numerators,denominators


        def plot_spots(I_ideal,I_real, miller_indices):
            plt.close()
            stringified = []
            for elem in miller_indices:
                stringified.append(str(elem[0])+str(elem[1])+str(elem[2])) 
            I_real *= np.sum(I_ideal)/np.sum(I_real) #normalise
            plt.bar(stringified,np.sqrt(I_ideal),alpha=1)
            plt.bar(stringified,np.sqrt(I_real),alpha=1,color='r',width=0.4)
            plt.ylim(0,np.sqrt(max(np.max(I_ideal),np.max(I_real))))
            plt.show()
        def plot_sectors(sector_histogram):            
            plt.close()
            fig = plt.figure()
            ax2 = fig.add_subplot(projection="polar",aspect="equal")
            sector_histogram = np.ma.masked_where(sector_histogram ==0, sector_histogram)  # masking for colour (replace with black)
            if full_range:
                pcolour = ax2.pcolormesh(sector_phi,sector_radial,sector_histogram,cmap=cmap,vmin=0,vmax=1)
            else: pcolour = ax2.pcolormesh(sector_phi,sector_radial,sector_histogram,cmap=cmap)
            fig.colorbar(pcolour)
            # print(sector_phi)
            # print()
            # print(sector_radial)
            # print()
            # print(sector_histogram)
            add_screen_properties()
            plt.show()

        if compare_handle != None:
            if log_I:
                print("Not plotting logarithmic I, not supported for comparisons.")
                log_I = False  
            
            all_I_ideal = []
            all_I_real = []
        
        # Iterate through each orientation (file) to get minimum/maximum for normalisation of plot of scattering image (not difference image).
        max_z = -np.inf
        min_z = np.inf

        fig = plt.figure()
        # fig.canvas.layout.width = '40%'
        # fig.canvas.layout.height = '40%'
        # fig.canvas       
        ax = fig.add_subplot(projection="polar")        
        for filename in os.listdir(results_dir):
            #result = get_result(filename)
            result1,result2 = get_result(filename,results_dir,compare_dir)
            if result1 == "__PASS__":
                continue            
            if result1 == None:
                break  
            num_orientations += 1
            if compare_dir != None:
                max_z = max(max_z,np.max(result1.R))
                min_z = min(min_z,np.min(result1.R))
                if max_z > 1:
                    print("error, max_R =",max_z)
                if min_z < 0:
                    print("error, min_R =",min_z)
            else:
                max_z = max(max_z,np.max(result1.I))
                min_z = min(min_z,np.min(result1.I))
        # Apply log (to non-diff image)
        if log_I:
            max_z = log_function(max_z)
            min_z = log_function(min_z)
            print("log(I) max,min ",max_z,min_z)
        else:
            print("I max,min",max_z,min_z)
        if max_z == min_z:
            print("Single intensity detected. Ignoring max/min z")
            min_z = max_z-1
        # Plot each orientation's scattering pattern/add each to sector histogram
        added_colorbar = False
        for filename in os.listdir(results_dir):
            result1,result2 = get_result(filename,results_dir,compare_dir)
            if result1 == "__PASS__":
                continue
            if result1 == None:
                break         
            # get points at same position on screen.
            radial_axis = result1.X*ang_per_bohr/1e7
            if plot_against_q:
                radial_axis = result1.q_scr/ang_per_bohr
            radial_axis = radial_axis[0]    

            # Plot spot intensities
            plot_all_the_spots = False
            if plot_all_the_spots:
                if result2 != None and crystal_aligned_frame:
                    plot_spots(result2.I.flatten(), result1.I.flatten(),result1.miller_indices)


            identical_count = np.zeros(result1.I.shape)
            if crystal_aligned_frame:
                phi = result1.phi_aligned
            else:
                phi  = result1.phi
            tmp_I1 = np.zeros(result1.I.shape)
            tmp_I2 = np.zeros(tmp_I1.shape)
            processed_copies = []
            unique_values_mask = np.zeros(result1.I.shape,dtype="bool")
            # Catch for multiple overlapping points (within same orientation only!!!) (does work, but would be unlikely.)
            # TODO need to do something about fact that overlapping points with many plots will hide points. Not critical rn thanks to sectors. 
            for i in range(len(result1.I)):
                if (radial_axis[i], phi[i],result1.image_index)  in processed_copies:
                    unique_values_mask[i] = False
                    continue
                else:
                    unique_values_mask[i] = True
                # if (radial_axis[i] > radial_lim):
                #     unique_values_mask[i] = False
                    
                # Average out spots (TODO need to check they are the same for no stochastic variation.)
                # TODO need to change this to be the bragg indices.
                matching_rad_idx = np.nonzero(radial_axis == radial_axis[i])  # numpy note: equiv. to np.where(condition). Non-zero part irrelevant.
                matching_phi_idx = np.nonzero(phi == phi[i])
                for elem in matching_rad_idx[0]:
                    if elem in matching_phi_idx[0]:
                        identical_count[i] += 1
                        matching1 = result1.I[(radial_axis == radial_axis[i])*(phi == phi[i])]
                        for value1 in matching1:
                            tmp_I1[i] += value1    
                        if result2 != None:
                            matching2 = result2.I[(radial_axis == radial_axis[i])*(phi == phi[i])]
                            for value2 in matching2:
                                tmp_I2[i] += value2      
            #print(identical_count)                                                
            tmp_I1 /= (identical_count)
            tmp_I2 /= (identical_count)

                #processed_copies.append((radial_axis[i],phi[i],result.image_index[i]))
            #print("processed copies:",processed_copies)

            # Get the dependent variable.
            result = copy.deepcopy(result1)
            result.I = tmp_I1 
            # get z used for colour of scatter plot.
            if compare_dir != None:
                result2.I = tmp_I2                    
                result.diff(result2)
                z = result.R[unique_values_mask]   # Making comparison, set z to be measure of difference 
            else:
                z = result.I[unique_values_mask]   # no comparison, z is intensity

            identical_count = identical_count[unique_values_mask]
            # Intensities used for histogram 
            I1 = tmp_I1[unique_values_mask]
            I2 = tmp_I2[unique_values_mask]
            radial_axis = radial_axis[unique_values_mask]
            phi = phi[unique_values_mask] + np.pi/2 # to align with the SPI pixel plot
            phi[phi>=np.pi] -= 2*np.pi

            #sector_histogram += get_histogram_contribution(z,result.phi,radial_axis)
            if result2 != None:
                ## Sectors/Fragmented rings ##
                numerators,denominators = get_sector_histogram_contribution(I1,I2,phi,radial_axis)
                # Add to sum of all orientations' results
                sector_num_histogram += numerators
                sector_den_histogram += denominators                

                ## Full rings ## Note: we are filling 2d arrays with same-valued arcs for each radius so that pcolormesh can plot circles. 
                # The numerators/denominators are effectively tiled, e.g.: numerators = np.tile(numerators,(1,num_arcs))
                numerators,denominators = get_sector_histogram_contribution(I1,I2,phi,radial_axis,phi_edges = np.array([-np.pi,np.pi]))
                R_num_histogram += numerators
                R_den_histogram += denominators            
                if neutze_R:
                    all_I_real.extend(I1)
                    all_I_ideal.extend(I2)


            if log_I: 
                z = log_function(z)      

            #debug_mask = (0.01 < radial_axis[0])*(radial_axis[0] < 100)
            #print(identical_count[debug_mask])
            #print(radial_axis[0][debug_mask])
            #print(phi[debug_mask ]*180/np.pi)
            if not get_R_only:
                # Dot size
                dot_param = z
                if crystal_pattern_only:
                    dot_param = identical_count
                if not log_dot:
                    dot_param = np.e**(dot_param)     
                norm = np.max(dot_param)
                s = [100*dot_size*x/norm for x in dot_param]
                if fixed_dot_size:
                    s = [100*dot_size for x in dot_param]

                # use cmap to get colours but with alpha following a specific rule
                COL = MplColorHelper(cmap, 0, 1) 
                alpha_modified_cmap_colours = np.empty((len(z),4)) 
                for i, K in enumerate(z):
                    try:
                        rgba = COL.get_rgb(K)
                    except Exception as e:
                        raise Exception("(max/min z:" +str(max_z) + "/" + str(min_z) + ") - val " + str(K) + " did not work for get_rgb().Original error: " + str(e))     
                    #Replace alpha
                    rgba=rgba[0:3] + (np.clip((K*(max_alpha-min_alpha))/(max_z) + min_alpha,min_alpha,None),)
                    alpha_modified_cmap_colours[i] = np.array(rgba)                           
                

                colours = [(r,g,b,a) for r,g,b,a in alpha_modified_cmap_colours] 
                colours = np.around(colours,10)   

                if not added_colorbar:
                    if compare_dir == None: 
                        print("Warning, colorbar not working with color power at present. Need to create cmap from COL")
                        if not normalise_intensity_map:
                            # Good for debugging
                            fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=np.min(z),vmax=np.max(z)),cmap=cmap),ax=ax)
                            print("Attention: Not normalising intensities, but still arbitrary units")
                            print("Warning: not yet taking into account combined dots!! Scale is off!")
                        else:
                            fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0,vmax=1),cmap=cmap),ax=ax) 
                    else:
                        #TODO get this truncation of bar working https://stackoverflow.com/questions/40982050/matplotlib-how-to-cut-the-unwanted-part-of-a-colorbar
                        # from matplotlib import colorbar
                        # colors = cmap(np.linspace(1.-(0.5-0.3)/float(0.5), 1, cmap.N))
                        # cbar_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(cmap, colors)
                        # cax,_ = colorbar.make_axes(ax)
                        # norm= mpl.colors.Normalize(vmin=0,vmax=1)
                        # cbar = colorbar.ColorbarBase(cax, cmap=cbar_cmap, norm=norm)
                        # cbar.set_ticks([0.3,0.4,0.5])
                        # cbar.set_ticklabels([0.3,0.4,0.5])
                        fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0,vmax=1),cmap=cmap),ax=ax) 
                            
                    #print(z)
                    added_colorbar = True

                sc = ax.scatter(phi,radial_axis,c=colours,s=s)
                plt.grid(alpha = min(show_grid,0.6),dashes=(5,10))   
                if show_labels:
                    rad_max = 0
                    q_max_idx = None
                    max_miller_index = 0
                    q_max = 0
                    for i, miller_indices in enumerate(result.miller_indices[unique_values_mask]):
                        #if rad_max < radial_axis[i]:
                        if q_max < result1.q[0][unique_values_mask][i]/ang_per_bohr:
                            rad_max = radial_axis[i]
                            q_max = result1.q[0][unique_values_mask][i]/ang_per_bohr
                            q_max_idx = i
                        max_miller_index = max(max_miller_index,np.max(miller_indices))
                        ax.annotate("%.0f" % miller_indices[0]+","+"%.0f" % miller_indices[1]+","+"%.0f" % miller_indices[2], (phi[i], radial_axis[i]),ha='center')
                    print("Num points", len(result.miller_indices[unique_values_mask]))
                    print("Max miller index",max_miller_index)
                    print("Miller indices of max q =","%.0f" % result.miller_indices[unique_values_mask][q_max_idx][0]+","+"%.0f" % result.miller_indices[unique_values_mask][q_max_idx][1]+","+"%.0f" % result.miller_indices[unique_values_mask][q_max_idx][2]) 
                    if plot_against_q:
                        print("max q_scr =",rad_max)
                    else:
                        print("max radius=",rad_max)
                    print("q max=", q_max)
        if not get_R_only:
            add_screen_properties()
            plt.show()  
        if compare_handle != None:
            if not get_R_only:
                #print("Plotting orientation-averaged R factor") # Doesn't work because when sector is empty it reduces the average.
                #plot_sectors(sector_histogram h= sector_histogram)
                non_zero_denon_histogram = sector_den_histogram.copy()
                non_zero_denon_histogram[non_zero_denon_histogram== 0] = 1        
                sector_histogram = np.divide(sector_num_histogram,non_zero_denon_histogram)
                print("Plotting total R factor (incorrect, need to normalise I's)")
                plot_sectors(sector_histogram = sector_histogram)

                #print("plotting full ring orientation-averaged R factor") # Doesn't work because when sector is empty it reduces the average.
                #plot_sectors(R_histogram)

                
                non_zero_denon_histogram = R_den_histogram.copy()
                non_zero_denon_histogram[non_zero_denon_histogram== 0] = 1        
                R_histogram = np.divide(R_num_histogram,non_zero_denon_histogram)       
                print("plotting full ring total R factor (incorrect, need to normalise I's)") 
                plot_sectors(R_histogram)   

            if neutze_R:
                sqrt_ideal = np.sqrt(all_I_ideal)
                sqrt_real = np.sqrt(all_I_real)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)
                R = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                # print("sqrt_ideal",sqrt_ideal)
                # print("-------")
                # print("-------")
                # print("sqrt_real normed",sqrt_real*inv_K) 

                print("R: ",R)
                # print("R:") (not normalised)
                # R = np.sum(np.abs((sqrt_real - sqrt_ideal)))/np.sum(sqrt_ideal)
                # print(R)
                if not get_R_only:
                    print("sum of (sqrt) real","{:e}".format(np.sum(sqrt_real)),"sum of (sqrt) ideal","{:e}".format(np.sum(sqrt_ideal)),"sum of abs difference","{:e}".format(np.sum(np.abs((sqrt_real - sqrt_ideal)))))
                    #neutze_histogram = np.histogram2d(phi, radial_axis, weights=R, bins=(np.array([-np.pi,np.pi]), radial_edges))[0]             
                    #plot_sectors(neutze_histogram)

                # Pearson Correlation Coefficient:
                x = all_I_ideal
                y = all_I_real
                x_bar = np.mean(x)
                y_bar = np.mean(y)   
                num = np.sum((x-x_bar)*(y-y_bar))
                den = np.sqrt(np.sum(x-x_bar)**2*np.sum(y-y_bar)**2)
                cc = num/den

                return R,cc,None # resolution lims not implemented
    
    ## Continuous (SPI)
    else:
        result1 = copy.deepcopy(SPI_result1)
        result2 = copy.deepcopy(SPI_result2)     
        if normalize_to_centre:
            # # Normalize I to 1 at centre 
            # # Average out centre pixels/select central pixels depending on axis dimension's integer parity.
            # A,B,C,D = ( int(np.floor((result1.I.shape[0]+1)/2)), int(np.ceil((result1.I.shape[0]+1)/2)), 
            #             int(np.floor((result1.I.shape[1]+1)/2)), int(np.ceil((result1.I.shape[1]+1)/2)) )
            # result1.I/=np.average(result1.I[A:B+1, C:D+1])
            # if result2 != None:
            #     result2.I/=np.average(result2.I[A:B+1, C:D+1])  
            result1.I/=result1.zero_angle_I
            result2.I/=result2.zero_angle_I
            result1.zero_angle_I = result2.zero_angle_I = 1
        else:
            # Normalize  total magnitude of intensity to 1.
            result1.I/=np.sum(result1.I)
            result1.zero_angle_I/=np.sum(result1.I)
            if result2 != None:
                result2.I/=np.sum(result2.I)
                result2.zero_angle_I/=np.sum(result2.I)
        ## Square grid
        if len(result1.I.shape) > 1:#if type(result1) == Results_Grid:
            if spi_full_rings_only:
                # for calculating R we remove the non-full rings - we represent this visually:
                result1.I *=  (result1.full_ring_mask + 0.1)/1.1         
                result2.I *= (result2.full_ring_mask+0.1)/1.1 

            if log_I: 
                z1 = log_function(result1.I)   
                result1.zero_angle_I = log_function(result1.zero_angle_I)
                z1[result1.I == 0] = None
                if result2 != None:
                    z2 = log_function(result2.I)
                    result2.zero_angle_I = log_function(result2.zero_angle_I)
                    z2[result2.I == 0] = None
            #else:
                #TODO fix up this masked stuff i didnt finish 
                # z1 -= cutoff_log_intensity
                # z1[z1<0] = 0
                #z1 = result1.I.copy()/np.sum(result1.I)    
                #if result2 != None:
                    #z2 -= cutoff_log_intensity
                    #z2[z2<0] = 0
                    #z2 = result2.I.copy()/np.sum(result2.I)        
            if log_I and cutoff_log_intensity != None:
                np.ma.array(z1, mask=(z1<cutoff_log_intensity)*(np.isnan(z2)))
                if result2 != None:
                    z2 = np.ma.array(z2, mask=(z2<cutoff_log_intensity)*(np.isnan(z2)))
            else:
                z1 = np.ma.array(z1,mask = np.isnan(z1))
                z2 = np.ma.array(z2,mask = np.isnan(z2)) 

            if spi_full_rings_only:
                # Now remove the non-full rings
                result1.I *=  result1.full_ring_mask          
                result2.I *= result2.full_ring_mask
                        
            if result2 != None:
                alpha2 = (result2.full_ring_mask + 1)/2             
            print("Result 1 (Damaged):")
            print("Total screen-incident intensity:","{:e}".format(np.sum(result1.I)))
            if result2 != None:
                combined_data = np.array([z1,z2])
            else: 
                combined_data = z1
            z_min, z_max = np.nanmin(combined_data), np.nanmax(combined_data)       
            if log_I:
                z_min = max(z_max-log_range,z_min)
            current_cmap = plt.colormaps.get_cmap("plasma")
            current_cmap.set_bad(color='black')
            if not get_R_only:
                z1_map = plt.imshow(z1,vmin=z_min,vmax=result1.zero_angle_I,cmap=current_cmap)
                plt.colorbar(z1_map)
                plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                plt.savefig("tmp_I_real.pdf",format="pdf",dpi=dpi)
                plt.show()
            if result2 != None:       
                if not get_R_only:
                    print("Result 2 (Undamaged):")
                    print("Total screen-incident intensity:","{:e}".format(np.sum(result2.I)))
                    z2_map = plt.imshow(z2,vmin=z_min,vmax=result2.zero_angle_I,cmap=current_cmap)
                    plt.colorbar(z2_map)
                    plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                    plt.savefig("tmp_I_ideal.pdf",format="pdf",dpi=dpi)
                    plt.show()
                    print("R:")
                    fig, ax = plt.subplots()
                    I_tmp = result2.I
                    if log_I:
                        I_tmp = z2.copy()
                    alpha = (I_tmp-np.nanmin(I_tmp))/(np.nanmax(I_tmp) - np.nanmin(I_tmp))
                    alpha[np.isnan(alpha)] = 0
                sqrt_real = np.sqrt(result1.I)
                sqrt_ideal = np.sqrt(result2.I)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # Scales I_real to I_ideal's tot intensity
                R_cells = np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal))
                R = np.sum(R_cells)
                print(R)               
                if not get_R_only:
                    #min_R_dmg_pixel = 0.1 # The minimum R factor contributed by a pixel (multiplied by the number of pixels), for it to be displayed.
                    alpha_prop_to_I = False
                    bg = np.full((*z1.shape, 3), 0, dtype=np.uint8) #bg = np.full((*z1.shape, 3), 70, dtype=np.uint8)

                    ax.imshow(bg)
                    R_cells *= len(R_cells)**2# multiply by num cells to give the 'weighted contribution' (such that R is now like the weighted average) 
                    alpha[R_cells < min_R_dmg_pixel] = 0  # hide insignificant pixels
                    if not alpha_prop_to_I:
                        alpha[alpha != 0] = 1  # override.
                    R_map = ax.imshow(R_cells,vmin=0,vmax=0.4,alpha=alpha,cmap=cmap)     
                    plt.colorbar(R_map)
                    ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                    ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                    ticks = ticks[len(ticks)//10::len(ticks)//5]
                    ticklabels = ticklabels[len(ticklabels)//10::len(ticklabels)//5]                    
                    plt.xticks(ticks,ticklabels)
                    plt.yticks(ticks,ticklabels)
                    plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                    plt.savefig("tmp_R_contributions.pdf",format="pdf",dpi=dpi)
                    plt.show()

                    R_num = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)))        

                    norm_root_diff_map = True
                    if norm_root_diff_map:  #TODO clean mess
                        print("Plotting normalised root difference map")
                        fig, ax = plt.subplots()
                        R_cells = np.abs((inv_K*sqrt_real - sqrt_ideal)/sqrt_ideal)
                        ax.imshow(bg)
                        R_map = ax.imshow(R_cells,vmin=0,vmax=2,alpha=alpha,cmap=cmap2)     
                        plt.colorbar(R_map)
                        ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                        ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                        ticks = ticks[len(ticks)//10::len(ticks)//5]
                        ticklabels = ticklabels[len(ticklabels)//10::len(ticklabels)//5]                        
                        plt.xticks(ticks,ticklabels)
                        plt.yticks(ticks,ticklabels)
                        plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                        plt.savefig("tmp_R_pixel.pdf",format="pdf",dpi=dpi)
                        plt.show()   
                        # print("Plotting change map") 
                        # fig, ax = plt.subplots()
                        # R_cells = inv_K*sqrt_real/sqrt_ideal
                        # ax.imshow(bg)
                        # R_map = ax.imshow(R_cells,vmin=None,vmax=2,alpha=alpha,cmap=cmap2)     
                        # plt.colorbar(R_map)
                        # ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                        # ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                        # plt.xticks(ticks,ticklabels)
                        # plt.yticks(ticks,ticklabels)
                        # plt.show()                          
                    print ("Plotting log ratio")
                    I1 = result1.I #real
                    I2 = result2.I #ideal
                    # Average out centre pixels/select central pixels depending on axis dimension's integer parity.
                    A,B,C,D = ( int(np.floor((I1.shape[0]+1)/2)), int(np.ceil((I1.shape[0]+1)/2)), 
                                int(np.floor((I1.shape[1]+1)/2)), int(np.ceil((I1.shape[1]+1)/2)) )
                    I1_divisor = np.average(I1[A:B+1, C:D+1])
                    I2_divisor = np.average(I2[A:B+1, C:D+1])                    
                        
                    log_ratio = log_function((I1/I1_divisor)/(I2/I2_divisor))
                    fig, ax = plt.subplots()
                    ax.imshow(bg)
                    norm = TwoSlopeNorm(vmin=log_diff_vmin, vcenter=0, vmax=log_diff_vmax)
                    I_map = ax.imshow(log_ratio,norm=norm,cmap = cmap2)#cmap="plasma")#"nipy_spectral_r")
                    cb = plt.colorbar(I_map)
                    cb.ax.set_yscale('linear')
                    ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                    ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                    ticks = ticks[len(ticks)//10::len(ticks)//5]
                    ticklabels = ticklabels[len(ticklabels)//10::len(ticklabels)//5]
                    plt.xticks(ticks,ticklabels)
                    plt.yticks(ticks,ticklabels)
                    plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                    plt.savefig("tmp_log_diff.pdf",format="pdf",dpi=dpi)
                    plt.show()                      


                    print("sum of real","{:e}".format(np.sum(sqrt_real)),"sum of ideal","{:e}".format(np.sum(sqrt_ideal)),"sum of abs difference","{:e}".format(R_num))       
                

                ### Calculate damage measures
                # Iterate through bins of resolutions
                resolutions = q_to_res(np.sqrt(np.apply_along_axis(np.sum,2, result1.q_xy**2))) 
                #min_res = np.min(resolutions)  # THIS IS NOT FROM THE RIM!!! It includes all.
                min_res = np.max([np.max(resolutions[0,:]),np.max(resolutions[:,0])])  # rim resolution
                max_res = np.min([6,np.max(resolutions)])
                res_lims = np.linspace(min_res,max_res,20)
                cc = np.zeros(len(res_lims)) # Pearson cc
                R_vect = cc.copy() # R_dmg (for full dataset up to the given resolution)
                binned_R_vect = cc.copy() # R_dmg for individual resolution bins/rings
                binned_q_R_dict = {}
                binned_q_I_ratio_dict  = {}
                # Generate dictionary of evenly spaced q, spanning the range of 0 to q corresponding to min_res. Bins centred on 1/2*delta_q, 3/2*delta_q and so on.
                delta_q = 0
                for i,val in enumerate(np.linspace(0,res_to_q(min_res),20,endpoint=False)):
                    binned_q_R_dict[val+0.5*delta_q] = np.nan
                    if i == 1:
                        delta_q = val
                        binned_q_R_dict.pop(0)
                        binned_q_R_dict.pop(val)
                        binned_q_R_dict[0.5*delta_q],binned_q_R_dict[1.5*delta_q] = [np.nan]*2
                binned_q_I_ratio_dict  = copy.deepcopy(binned_q_R_dict)                        
                
                for i in range(len(res_lims)):
                    lower = resolutions * 0 + res_lims[i]
                    upper = np.Infinity
                    # Select data up to some resolution limit
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    # cc
                    x_bar = np.mean(x)
                    y_bar = np.mean(y)   
                    num = np.sum((x-x_bar)*(y-y_bar))
                    den = np.sqrt(np.sum((x-x_bar)**2)*np.sum((y-y_bar)**2))
                    cc[i] = None
                    if den != 0:
                        cc[i] = num/den
                    # R factor
                    sqrt_real = np.sqrt(x)
                    sqrt_ideal = np.sqrt(y)
                    inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # normalises I_real to I_ideal's tot intensity
                    R_vect[i] = None
                    if np.sum(sqrt_ideal) != 0:
                        R_vect[i] = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))                    
                    # R factor but binned to resolution (note the resolutions still correspond to the maximum scattering angles for the bin, i.e. the bin is not centred on the corresponding q for the resolution.)
                    delta_res = res_lims[1]-res_lims[0]
                    lower = resolutions * 0 + res_lims[i]
                    upper = resolutions * 0 + res_lims[i] + delta_res/2
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    binned_sqrt_real = np.sqrt(x)
                    binned_sqrt_ideal = np.sqrt(y)
                    binned_inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # normalises I_real to I_ideal's tot intensity
                    binned_R_vect[i] = None                    
                    if np.sum(binned_sqrt_ideal) != 0:
                        binned_R_vect[i] = np.sum(np.abs((binned_inv_K*binned_sqrt_real - binned_sqrt_ideal)/np.sum(binned_sqrt_ideal)))
                # R factor but bins separated by constant delta q.
                for key in binned_q_R_dict.keys():
                    lower = resolutions*0 + q_to_res(key+delta_q/2)
                    upper = resolutions*0 + q_to_res(key-delta_q/2)
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    #tmp = copy.deepcopy(result1.I)
                    #tmp[tmp == 0] = np.nan
                    #plt.imshow(tmp)
                    sqrt_real = np.sqrt(x)
                    sqrt_ideal = np.sqrt(y)
                    inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # normalises I_real to I_ideal's tot intensity
                    binned_q_R_dict[key] = None   
                    if np.sum(sqrt_ideal) != 0:
                        binned_q_R_dict[key] = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))                         
                # Average intensity
                I1_centre = result1.I[int(len(result1.I)/2)][int(len(result1.I)/2)]
                I2_centre = result2.I[int(len(result2.I)/2)][int(len(result2.I)/2)]
                for key in binned_q_I_ratio_dict.keys():
                    lower = resolutions*0 + q_to_res(key+delta_q/2)
                    upper = resolutions*0 + q_to_res(key-delta_q/2)
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    if np.sum(y) != 0:
                        binned_q_I_ratio_dict[key] = log_function(np.mean( (x/I1_centre)/(y/I2_centre) ))                    

                #plt.plot((bin_lims[0:-1] + bin_lims[1:])/2,cc)
                if not get_R_only:
                    plt.rcParams["text.usetex"] = True
                    plt.plot(res_lims,R_vect,color="red")
                    plt.ylabel("R")
                    plt.xlabel("d ($\AA$)")                    
                    plt.ylim(0,1.04)
                    plt.gca().invert_xaxis()
                    plt.show()                    
                    plt.plot(res_lims,cc)
                    plt.ylabel("CC")
                    plt.xlabel("d ($\AA$)")
                    plt.ylim(0,0.2)
                    plt.gca().invert_xaxis()
                    plt.show()
                damage_dict  = dict(
                    R = R_vect,
                    bin_res_R = binned_R_vect,
                    cc = cc,
                    resolutions = res_lims,
                    bin_q_R = binned_q_R_dict,
                    bin_q_I = binned_q_I_ratio_dict
                )
                return damage_dict
            
        ## Circle grid
        else:
            if log_I: 
                z = log_function(result1.I)
            else:
                z = result1.I        
            if log_I and cutoff_log_intensity != None:
                #cutoff_log_intensity = -1
                z -= cutoff_log_intensity
                z[z<0] = 0
                pass
            radial_axis = result1.X
            if plot_against_q:
                radial_axis = result1.q         
            fig = plt.figure()
            ax = fig.add_subplot(projection="polar")
            # phi_mesh = result.phi_mesh
            # if crystal_aligned_frame:
            #     phi_mesh = result.phi_aligned_mesh   
            #     print("crystal aligned not implemented for SPI yet")    
            ax.pcolormesh(result1.phi, radial_axis, z,cmap=cmap)
            ax.plot(result1.phi, radial_axis, color = 'k',ls='none')
            plt.grid(alpha=min(show_grid,0.6),dashes=(5,10)) 

            if radial_lim:
                bottom,top = plt.ylim()
                plt.ylim(bottom,radial_lim)
            if log_radial:
                plt.yscale("log")

            if result2 != None:
                print("R:")
                sqrt_real = np.sqrt(result1.I)
                sqrt_ideal = np.sqrt(result2.I)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real) 
                R = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                print(R)   

                # Pearson Correlation Coefficient:
                x = result1.I
                y = result2.I
                x_bar = np.mean(x)
                y_bar = np.mean(y)   
                num = np.sum((x-x_bar)*(y-y_bar))
                den = np.sqrt(np.sum(x-x_bar)**2*np.sum(y-y_bar)**2)
                cc = num/den                             
                return R,cc,None  # resolution lims not implemented
                # print("R: (not normalised)")
                # R = np.sum(np.abs((sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                # print(R)

    def get_orientation_set_of_folder():
        '''Returns a list of orientations present in results subfolder'''

# https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def get_result(filename,results_dir,compare_dir = None):
    #Requires all orientations of result_handle in compare_handle, but not vice versa.
    fpath = os.path.join(results_dir, filename)
    if os.path.isfile(fpath):
        with open(fpath,'rb') as f:
            result1 = pickle.load(f)
    else: 
        return "__PASS__" "__PASS__"
    result2 = None
    if compare_dir != None:
        if filename in os.listdir(compare_dir):
            fpath2 = os.path.join(compare_dir, filename)
            with open(fpath2,'rb') as f:
                result2 = pickle.load(f)     
                result1.diff(result2)             
        else:
            print("ERROR, missing matching orientation in comparison directory")
            return None, None # No corresponding file found.          
    return result1, result2




##### https://scripts.iucr.org/cgi-bin/paper?S0021889807029238, http://superflip.fzu.cz/

#TODO figure out why we get zeros for reflection intensities sometimes.
import pandas as pd
import csv
def create_reflection_file(result_handle,results_parent_dir = RESULTS_LOCAL_PATH,overwrite=False):
    '''
    Generates a .rfl file, compatible with Superflip
    '''
    print("Creating reflection file for",result_handle)
    directory = "reflections/"
    results_dir = results_parent_dir+ result_handle+"/"
    os.makedirs(directory, exist_ok=True) 
    results_dir = results_parent_dir+result_handle+"/"
    init = False
    for filename in os.listdir(results_dir):
        result = get_result(filename,results_dir)[0]
        if result == "__PASS__":
            continue
        if result == None:
            break          
        miller_indices = result.miller_indices
        intensity = np.array([result.I]).T
        data = np.concatenate((miller_indices,intensity),axis=1)
        # Put in data frame
        columns = ["h","k","l","I"]
        new_df = pd.DataFrame(data=data, columns=columns)
        if init == False:
            df = new_df
            init = True
        else:
            df = pd.concat((df,new_df))
    # Save as file
    out_path = directory + result_handle + ".rfl"
    if os.path.isfile(out_path): 
        if not overwrite:
            print("Cannot write, file already present at",out_path)
            return
        os.remove(out_path)
    columns.reverse()
    df = df.sort_values(by=["l","k","h"],axis=0)
    for i in ("hkl"):
        df[i] = df[i].astype('int')
    df["I"] = df["I"].astype('float')
    df = df.round(6)
    df.drop_duplicates(subset = ["h","k","l"],inplace=True) # TODO should take average.
    df = df[df['I']>=0.01] # TODO temporary fix for appearance of low values that needs to be squashed.
    df.to_csv(out_path,header=False,index=False,float_format='%10f', sep=" ", quoting=csv.QUOTE_NONE, escapechar=" ")

def rfl_to_sca(result_handle,reflections_dir = "reflections/",overwrite=True):
    '''Converts .rfl file to scalepack .sca file
    See https://www.ccp4.ac.uk/html/scala.html#files
    '''
    directory = "scalepack/"
    os.makedirs(directory, exist_ok=True) 
    rfl_file_path = reflections_dir + result_handle + ".rfl"
    out_path = directory + result_handle + ".sca"
    save_action = "x"
    if overwrite:
        save_action = "w"    
    # First find max length of intensities, so can scale values down.
    max_length = 0
    with open(rfl_file_path, 'r') as a:
        for line in a:
            entries = line.split()
            if entries[0]+entries[1]+entries[2] == '0'*3: 
                continue # Ignore (0,0,0) reflection
            
            max_length = max(max_length,len(line.split()[3].split('.')[0]))
    #
    with open(rfl_file_path, 'r') as a, open(out_path, save_action) as b:
        # placeholder boilerplate 
        indent = ' '*3
        b.write(indent+' 1\n -987\n')
        # Cell geometry TODO integrate with actual input.
        b.write(indent+' 79.000    79.000    38.000    90.000    90.000    90.000 P43212\n')
        # Miller indices (hkl) | IMEAN_dataset | SIGIMEAN_dataset
        for line in a:
            # Initialise elements
            h=k=l = ' '*4
            I_mean=sigI_mean = ' '*8
            I_scaling_power = max_length - 7

            # convert data to usable format from .rfl file
            entries = line.split()
            #   hkl
            if entries[0]+entries[1]+entries[2] == '0'*3:
                continue # Ignore (0,0,0) reflection
            #   I
            entries[3] =  '%.0f'%(float(entries[3])/10**I_scaling_power)
            # Add in entry for sigmaI_mean
            entries.append('50.0')

            # Populate elements
            for i,q in enumerate([h,k,l,I_mean,sigI_mean]):
                assert len(entries[i]) < len(q), "Error, entry of '"+entries[i]+"' has length >= "+str(len(q))   # Use '<' not '<=' because need a space.
                q = q[:-len(entries[i])]+ entries[i]
                b.write(q)
            b.write('\n')


#create_reflection_file("hen_v7__eal",True)
#rfl_to_sca("hen_v7_real")

#####
# stylin' 
def stylin(exp_name1,exp_name2,radial_lim,get_R_only = False,SPI=False,SPI_max_q=None,SPI_result1=None,SPI_result2=None,results_parent_dir = RESULTS_LOCAL_PATH,show_labels=False,**kwargs):
    experiment1_name = exp_name1#"Lys_9.95_random"#exp_name1
    experiment2_name = exp_name2#"lys_9.80_random"#exp_name2 

    #####


    font = {'family': 'serif',
            'size'   : 10}

    plt.rc('font', **font)

    use_q = False # TODO remove this option
    log_radial = False
    log_I = True
    cutoff_log_intensity = -1#-1
    try:
        cmap = shiftedColorMap(matplotlib.cm.RdYlGn_r,midpoint=0.2,name="shiftedcmap")#"plasma"#"YlGnBu_r"#cc.m_fire#"inferno"#cmr.ghostlight#cmr.prinsenvlag_r#cmr.eclipse#cc.m_bjy#"viridis"#'Greys'#'binary'
    except: 
        cmap =  plt.get_cmap("shiftedcmap")
    cmap.set_bad(color='black')
    cmap_power = 1.6
    min_alpha = 0.3
    max_alpha = 1
    colour = "y"
    full_crange_sectors = False

    cmap_intensity = "inferno"


    # screen_radius = 150#55#165    #
    # q_scr_lim = experiment.r_to_q(screen_radius) #experiment.r_to_q_scr(screen_radius)#3.9  #NOTE this won't be the actual max q_parr_screen but ah well.
    # zoom_to_fit = True
    # ####### n'
    # # plottin'

    # as q = ksin(theta).  
    #experiment1.max_q*(1/np.sqrt(2)) # for flat screen. as q_scr = qcos(theta). (q_z = qsin(theta) =ksin^2(theta)), max theta is 45. (Though experimentally ~ 22 as of HR paper)

    zoom_to_fit = False
    if not zoom_to_fit:
        radial_lim*=1.02
        #els:
    #else:
    #     radial_lim = screen_radius#min(screen_radius,experiment.q_to_r(experiment.max_q))
    #     print(radial_lim)
    #     if use_q:
    #         #radial_lim = experiment.r_to_q_scr(radial_lim)
    #         radial_lim = q_scr_lim #radial_lim = min(q_scr_lim,experiment.q_to_q_scr(experiment.max_q))
    #         print(radial_lim)
    # else:
    #     radial_lim = None

    #TODO fix above to work with distance

    # R Sectors
    if not SPI:
        if not get_R_only:
            print("----R Sectors unaligned----")
            scatter_scatter_plot(crystal_aligned_frame = False,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True,results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs) 
            print("----Intensity of experiment 1----")
            scatter_scatter_plot(crystal_aligned_frame = False,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print("----Intensity of experiment 2----")
            scatter_scatter_plot(crystal_aligned_frame = False,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment2_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print("----R Sectors aligned----")
            damage_dict = scatter_scatter_plot(crystal_aligned_frame = True,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print("----Intensity of experiment 1 (damaged) aligned----") 
            scatter_scatter_plot(crystal_aligned_frame = True,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print("----Intensity of experiment 2 (undamaged) aligned----")
            scatter_scatter_plot(crystal_aligned_frame = True,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment2_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
        else:
           damage_dict = scatter_scatter_plot(get_R_only = True,crystal_aligned_frame = True,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)

    else:
        use_q = True
        log_radial = False
        log_I = True
        cutoff_log_intensity = None # -1 or None
        cmap = "PiYG_r"# "plasma"
        cmap2 = "seismic"
        radial_lim = SPI_max_q
        damage_dict = scatter_scatter_plot(get_R_only=get_R_only,log_I = log_I, cutoff_log_intensity = cutoff_log_intensity, SPI_result1=SPI_result1,SPI_result2=SPI_result2,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,cmap2=cmap2,**kwargs)


    return damage_dict

def res_to_q(d):
    '''
    Pass d = None to use default resolution (when put into exp_args)
    '''
    if d == None:
        return None
    return 2*np.pi/d
def q_to_res(q):
    return 2*np.pi/q

#TODO 
# log doesnt work atm.
# Get the rings to correspond to actual rings

#%% vvvvv Scattering Playground vvvv
DEBUG = False

if __name__ == "__main__":
    fig_width = 3.49751 # 20
    fig_height = fig_width*3/4 # 20
    ### Simulate
    target_options = ["neutze","hen","tetra","glycine","fcc"]
    #============------------User params---------==========#

    target = "hen"#"glycine"  #target_options[2]
    best_resolution = 2 # 1.58 (abdullah) # 2   # resolution (determining max q)
    worst_resolution = None#30 # 'resolution' corresponding to min q

    #### Individual experiment arguments 
    tag = "D" # Non-SPI i.e. Crystal only, tag to add to folder name. Reflections saved in directory named version_number + target + tag named according to orientation .
    start_time = -18#-12#-6
    end_time = 18#12#6
    laser_firing_qwargs = dict(
        # pixel sampling method (Neutze) if True - Miller indices if False
        SPI = True,  # sampling method, if False, bragg spots. if True, detector pixels. TODO change name
        SPI_resolution = best_resolution,
        pixels_across = 300,  # for SPI, shld go on xfel params.
        # miller
        random_orientation = True, #bragg spot sampling only, TODO refactor to be in same place as other orients...# orientation is synced with second 
    )
    ##### Crystal params
    crystal_qwargs = dict(
        supercell_scale = 1,  # for SC: supercell_scale^3 "unit" cells per supercell # Bragg spots will be sampled based on the cell scale, not the supercell scale.
        num_supercells = 1,#100, # 35409
        supercell_simulations = 1, #150
        positional_stdv = 0,#0.2,  #Introduces disorder to positions. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if gauging serial crystallography R factor, as should average out. 0.2 neutze.
        include_symmetries = True,  # should unit cell contain symmetries?
        cell_packing = "SC",
        rocking_angle = 0.02,  #  (approximating mosaicity - use 0.02 for proper, use a high value, like 1-10, and set a low max triple miller indice to disallow seemingly impossible indices (due to rocking angle/our implementation of it via momentum conservation formulae) that mimic studies that use the first few miller indices )
        #CNO_to_N = True,   # whether the plasma simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    )

    #### XFEL params
    #TODO make it so reflections don't overwrite same orientation, as stochastic now.
    energy = 7112#7100 # eV
    exp_qwargs = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(worst_resolution),#None #angstrom
        q_cutoff = res_to_q(best_resolution), #(best_resolution),#2*np.pi/2
        t_fineness=25,   
        #####crystal stuff (miller)
        max_miller_idx = None, #None, # = m, [thus max q given by q with miller indices (m,m,m)]
        all_miller_indices = False, # False, whether to find all bragg points at or below the max miller index (and between min and max q)
        ####SPI stuff ( ab initio)
        num_rings = 20,
        pixels_per_ring = 20,
        # first image orientation cardan angles [degrees] 
        SPI_x_rotation = 0,
        SPI_y_rotation = 0,
        SPI_z_rotation = 0,
        #crystallographic orientations (not consistent with SPI yet)
        # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Overridden if random orientations.        
        num_orients_crys=5, # Miller indices orientations
        orientation_axis_crys = None,
        #orientation_axis_crys = [0,0,1],#None,#[1,1,0]
        
        # for debugging/comparison with other works
        custom_cell_dims_for_miller_indices = None,#[17.174,14.93,13.384], # None # Implemented for comparison with others that use supercells.  
        override_max_q = False # False # Also special, implemented for comparison purposes but should be left as False by default.
        ######
    )
    same_deviations = True # whether same position deviations between damaged and undamaged crystal (SPI only) 
    

    # Optional: Choose previous folder for crystal results
    chosen_root_handle = None # None for new. use e.g. "tetra_v1", if want to add images under same params to same results.
    #=========================-------------------------===========================#

    ## DEBUG
    # WARNING we often assume that first crystal is damaged and second is undamaged when plotting. 
    first_crystal_is_damaged = True # True  
    second_crystal_is_damaged = False  # False

    #---------------------------Result handle names---------------------------#
    exp1_qualifier = "real"
    exp2_qualifier = "ideal"
    if chosen_root_handle is None:
        version_number = 1
        count = 0
        if tag != "":
            tag = "_" + tag        
        while True:
            if count > 99:
                raise Exception("could not find valid file in " + str(count) + " loops")
            results_parent_folder = RESULTS_LOCAL_PATH # needs to be synced with other functions
            root_handle = str(target) + tag
            exp_name1 = root_handle + "_" + exp1_qualifier + "_v" + str(version_number)
            exp_name2 = root_handle + "_" + exp2_qualifier  + "_v" + str(version_number)
            if path.exists(path.dirname(results_parent_folder + exp_name1 + "/")) or path.exists(path.dirname(results_parent_folder + exp_name2 + "/")):
                version_number+=1
                count+=1
                continue 
            break
    else:
        exp_name1 = chosen_root_handle + "_" + exp1_qualifier
        exp_name2 = chosen_root_handle + "_" + exp2_qualifier

    #exp_name2 = None

    #---------------------------------#
    water_index = None # None TODO automate
    if target == "neutze": #T4 virus lys
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/2lzm.pdb"
        target_handle = "lys-1_2"  
        folder = "lys"
        allowed_atoms = ["N_fast","S_fast"]
        CNO_to_N = True
    elif target == "hen": # egg white lys
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/4et8.pdb"
        # Solvated targets
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/sol_4et8_full_struct_asym.xpdb"; water_index = 1089
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/sol_4et8_full_struct_unit_cell.pdb"; water_index = 8705
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/lys_asym_water.xpdb"; water_index = 
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/lys_8_cell.xpdb"; water_index = 69632      
        # target_handle = "lys_nass_2"
        # folder = "lys"
        #'''
        # Light + Heavy atoms handles.
        target_handle = "lys_nass_Gd_full_1"#"lys_nass_Gd_full_1"#"lys_nass_no_S_3" #"lys_all_light-typical"#"lys_full-94_1"#"lys_nass_15"#"lys-5_3"#" #12keV, 0.1/0.01 count, 10 fs
        # Light atoms only
        #target_handle = "lys_nass_no_S_3"
        #folder = "lys" 
        folder = "" 
        #background_targets = "lys_water"
        #''' 
        '''
        target_handle = "lys_no_S_1"#"lys_no_S_2" #12keV, 0.1/0.01 count, 10 fs
        folder = ""        
        '''
        #//
        allowed_atoms = ["C","N","O"]
        #allowed_atoms = ["C","N","O"]
        #allowed_atoms = ["N","S_fast"]
        #allowed_atoms = ["N_fast"]
        #allowed_atoms = ["S_fast"]
        CNO_to_N = False
        S_to_N = True
        #//
        CNO_to_N = False
    elif target == "tetra": 
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/5zck.pdb" 
        folder = ""#"tetra_CNO"
        target_handle = "lys_all_light-typical"#"6-5-2_tetra_CNO_3"
        #allowed_atoms = ["N_fast"]
        allowed_atoms = ["C","N","O"]
        CNO_to_N = False
        S_to_N = False
    elif target == "glycine":
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/glycine.pdb" 
        folder = ""
        target_handle = "glycine_abdullah_4"
        allowed_atoms = ["C","N","O"]
        CNO_to_N = False
        S_to_N = False
    else:
        raise Exception("'target' invalid")
    #-------------------------------#

    import inspect
    src_file_path = inspect.getfile(lambda: None)
    sim_data_dir = path.abspath(path.join(src_file_path ,"../../../output/__Molecular/"+folder)) + "/"


    # Set up experiments
    experiment1 = XFEL(exp_name1,energy,**exp_qwargs)
    experiment2 = XFEL(exp_name2,energy,**exp_qwargs)
    # Create Crystals

    crystal = Crystal(pdb_path,allowed_atoms,is_damaged=first_crystal_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, **crystal_qwargs)
    # The undamaged crystal uses the initial state but still performs the same integration step with the pulse profile weighting.
    if same_deviations:
        # we copy the other crystal so that it has the same deviations in coords
        crystal_undmged = copy.deepcopy(crystal)#Crystal(pdb_path,allowed_atoms,cell_dim,is_damaged=False,CNO_to_N = CNO_to_N, **crystal_qwargs)
        crystal_undmged.is_damaged = second_crystal_is_damaged
    else:
        crystal_undmged = Crystal(pdb_path,allowed_atoms,is_damaged=second_crystal_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, **crystal_qwargs)
    crystal.plot_me(300000,water_index = water_index,template="plotly_dark")
#%
    if laser_firing_qwargs["SPI"]:
        SPI_result1 = experiment1.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal,results_parent_dir=results_parent_folder, **laser_firing_qwargs)
        SPI_result2 = experiment2.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal_undmged,results_parent_dir=results_parent_folder,  **laser_firing_qwargs)
        stylin(exp_name1,exp_name2,experiment1.max_q,SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2,custom_fig_width=fig_width,custom_fig_height=fig_height)
    else:
        exp1_orientations = experiment1.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal, results_parent_dir=results_parent_folder, **laser_firing_qwargs)
        create_reflection_file(exp_name1,results_parent_dir=results_parent_folder)
        rfl_to_sca(exp_name1)
        if exp_name2 != None:
            laser_firing_qwargs["random_orientation"] = False
            experiment2.set_orientation_set(exp1_orientations)  # pass in orientations to next sim, random_orientation must be false!
            experiment2.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal_undmged, results_parent_dir=results_parent_folder, **laser_firing_qwargs)
            create_reflection_file(exp_name2,results_parent_dir=results_parent_folder)
            rfl_to_sca(exp_name2)

        stylin(exp_name1,exp_name2,experiment1.q_to_X(experiment1.max_q)/1e7,custom_fig_width=fig_width,custom_fig_height=fig_height) # Note we are passing the max q, not max q_scr.
#^^^^^^^
#%% Pixels/SPI
if __name__ == "__main__":
    #fig_width = 3.49751 # 20
    #fig_height = fig_width*3/4 # 20
    stylin(exp_name1,exp_name2,experiment1.max_q,SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2, custom_fig_height=fig_height,custom_fig_width=fig_width,
           min_R_dmg_pixel=0,spi_full_rings_only=False,log_range=None,
           log_diff_vmin = -0.4, 
           log_diff_vmax = 0.6, dpi=800,)
#%% Miller/macrocrystal
if interactive and __name__ == "__main__":
    stylin(exp_name1,exp_name2,experiment1.q_to_X(experiment1.max_q)/1e7,show_labels=False) # Note we are passing the max q, not max q_scr.
    #stylin("glycine_v36__real","glycine_v36__ideal",2.3)
    #stylin(exp_name1,exp_name2,3)
#%%--------STRUCTURE CONSTRUCTOR------

# Save full structures in pdb format for SOLVATE
# Using this structure is not amazing practice, it takes a lot of time and potentially memory!
# SOLVATE allows for generating just the water with the solute removed. So an alternative method might 
# be generating the water for an individual unit cell with different distributions (seems possible by using slightly different thickness), 
# stitching them together, and removing atoms that are outside the Wigner–Seitz cell.
# We then calculate the form factor for the crystal, followed by the form factor for each of N water cells by defining water_background = Crystal(water_background_N,allowed_atoms).
# Finally, the rest of the water drop could be calculated by generating a large distribution of water, then scaling its contribution to the form factor.
if interactive and __name__ == "__main__":
    ##### Crystal params
    pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/4et8.pdb"
    crystal_qwargs = dict(
        supercell_scale = 1,  # for SC: supercell_scale^3 unit cells
        positional_stdv = 0,  # Not used
        include_symmetries = True,  # should unit cell contain symmetries or just one asymmetric unit?
        cell_packing = "SC",
        rocking_angle = 0.1,  # (approximating mosaicity - use 0.02 for proper, use a high value, like 1, and set a low max triple miller indice to disallow seemingly impossible indices (due to rocking angle/our implementation of it via momentum conservation formulae) that mimic studies that use the first few miller indices )
        CNO_to_N = False,
    )
    allowed_atoms = ["C","N","O","S"]
    crystal = Crystal(pdb_path,allowed_atoms,is_damaged=False, **crystal_qwargs)
    crystal.save_structure()

# %%
def plot_recovered_atoms():
    # Plot atoms retrieved from Superflip/EDMA
    df = pd.read_csv('tet_plot_points.pl',delim_whitespace=True)
    x = df['x']
    y = df['y']
    z = df['z']


    xscale = 4.813 
    yscale = 17.151
    zscale = 29.564
    fig =plt.figure(figsize =(zscale,yscale))
    #ax = Axes3D(fig)
    #ax.scatter(x,y,z)


    x*= xscale
    y*= yscale
    plt.scatter(z,y,c=-x,cmap="Blues",s=200)
    plt.xlabel("z (Ang)")
    plt.ylabel("y (Ang)")
if interactive and __name__ == "__main__":
    plot_recovered_atoms()    