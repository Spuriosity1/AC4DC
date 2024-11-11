import sys
import numpy as np
import argparse
import os.path as path
import os
from types import SimpleNamespace 

#TODO need to adjust num time steps and grid update period based on fwhm

# parser = argparse.ArgumentParser(description = "Empty Description")
# parser.add_argument("infile_path")
# parser.add_argument("outfile_directory")
# parser.add_argument("-f", "--fast", help = "Uses shell approximation for atoms", required = False, default = False)
# parser.add_argument("-md", "--MD", help = "Generate MD file", required = False, default = False)
# parser.add_argument("-n", "--name", help = "Custom name for generated .mol file (no ext)", required = False, default = '')
# input_argument = parser.parse_args()

def main():
  # Choose set of parameters to vary
  PULSE_PARAMETERS = 0; SOURCE_PARAMETERS = 1
  ####
  MODE = 0
  ####
  GAUSSIAN = 0; SQUARE = 1

  if MODE is PULSE_PARAMETERS:
    GRID_TYPE = "M"
  # https://www.xfel.eu/sites/sites_custom/site_xfel/content/e35165/e46561/e46876/e179573/e179574/xfel_file179576/19042023_Parameters_factsheet_2024-01_Final_eng.pdf
    ENERGIES = [6000,7000,8000,9000,10000,11000,12000,13000,14000,15000]
    FWHMS = [5,10,20,25]#[15]
    PHOTON_COUNTS = [1]#[0.1,10**(0.5)/10,1,10**(0.5), 10]

    SOURCE_FRACTION = None
    SOURCE_ENERGY = None
    SOURCE_DURATION = None
    PULSE_SHAPE = GAUSSIAN


  if MODE is SOURCE_PARAMETERS:
     GRID_TYPE = "S"
     ELECTRON_SOURCE = True 
     ENERGY = 10000
     FWHM = 15
     PHOTON_COUNT = 1000

     SOURCE_FRACTIONS = [3]
     SOURCE_ENERGIES = [0,50, 500, 1000, 1500, 2000, 3000, 4000, 5000,6000, 7000, 8000,9000, 10000,11000,12000,13000, 14000, 15000, 16000,17000,18000,19000,20000]
     SOURCE_DURATIONS = [1]
     PULSE_SHAPE = GAUSSIAN

  ##########


  batch_folder_parent_dir = "input/_batches/"

  if len(sys.argv) < 2:
    print('Usage: python3 generate_batch.py /path/to/file/to/copy.mol')
    return 1
  
  if sys.argv[1].partition('.')[-1] != 'mol':
    print('Invalid file: argument does not end in ".mol".')
    return 1

  # User provides path of input file
  infile_path = sys.argv[1] 
  handle = ''.join(path.basename(infile_path).partition('.')[:-1])[:-1]
  # We create the batch folder
  batch_dir = batch_folder_parent_dir + "batch_" + handle + "/"
  if not path.exists(batch_dir):
    os.makedirs(batch_dir)

  atoms, secondary_ion_exclusions, volume = copy_params(infile_path,handle)
  i = 1

  pulse_shape = None
  if PULSE_SHAPE is GAUSSIAN:
     pulse_shape = "gaussian"
  if PULSE_SHAPE is SQUARE:
     pulse_shape = "square"
  if MODE is PULSE_PARAMETERS:
     
    param_space = np.stack(np.meshgrid(FWHMS,PHOTON_COUNTS,ENERGIES),-1).reshape(-1,3)  # Note that last variable is the "inner loop", iterated over first and the most.
    feedin_dicts = [{
        'grid_type':GRID_TYPE,
        'atoms' : atoms,
        'scndry_exclusions' : secondary_ion_exclusions,
        'volume' : volume,
        'pulse_shape':pulse_shape, 
        'energy' : params[2],
        'photon_count' : params[1],    
        'fwhm' : params[0],
        'source_energy':SOURCE_ENERGY,
        'source_fraction':SOURCE_FRACTION,           
        'source_duration':SOURCE_DURATION,
       } for params in param_space] 
  if MODE is SOURCE_PARAMETERS:
    param_space = np.stack(np.meshgrid(SOURCE_DURATIONS,SOURCE_FRACTIONS,SOURCE_ENERGIES),-1).reshape(-1,3)

    assert(ELECTRON_SOURCE is True)
    feedin_dicts = [{
          'grid_type':GRID_TYPE,
          'atoms' : atoms,
          'scndry_exclusions' : secondary_ion_exclusions,
          'volume' : volume,
          'pulse_shape':pulse_shape, 
          'energy' : ENERGY,
          'photon_count' : PHOTON_COUNT,    
          'fwhm' : FWHM,
          'source_energy':params[2],
          'source_fraction':params[1],    
          'source_duration':params[0],   
      } for params in param_space]     
     
  for feedin_dict in feedin_dicts:
      outfile_path = None
      # Get unoccupied path.
      while outfile_path == None or path.exists(outfile_path):
        outfile_name = handle + "-" + str(i)
        outfile_path = batch_dir + outfile_name + '.mol'
        i+=1
      make_mol_file(infile_path, outfile_path,feedin_dict)

  return 0

def copy_params(mol_path,handle):
    '''
    Reads the control file and returns the relevant parameters within
    '''    
    input_path = "input/"
    molfile = find_mol_file_from_directory(input_path,handle)
    atom_lines = []
    exclusion_lines = []
    volume = None    
    with open(molfile, 'r') as f:
        reading_atoms = False
        reading_volume = False
        reading_secondary_ion_exclusions = False
        n = 0
        for line in f:
            if line.startswith("#ATOMS"):
                reading_atoms = True
                continue             
            if line.startswith("#BOUND_FREE_EXCLUSIONS"):
                reading_secondary_ion_exclusions = True
                continue                     
            elif line.startswith("#VOLUME"):
                reading_volume=True
                continue
            elif line.startswith("#") or line.startswith("//") or len(line.strip()) == 0:
                reading_atoms = False
                reading_volume = False
                reading_secondary_ion_exclusions = False
                n = 0
                continue
            if line.startswith("####END####"):
                break
            if reading_atoms:
                atom_lines.append(line)
            if reading_secondary_ion_exclusions:
               exclusion_lines.append(line)
            if reading_volume:
                if n == 0:
                    volume = float(line.split(' ')[0])
                n += 1
    print("Atoms:",atom_lines)
    if len(exclusion_lines) > 0:
      print("Exclusions:",exclusion_lines)
    print("Volume:", volume)
    return atom_lines, exclusion_lines, volume

def make_mol_file(fname, outfile, param_dictionary):
  
  P = SimpleNamespace(**param_dictionary)

  # Check validity of inputs
  if P.source_fraction in [None,0]:
     assert all(var in [None,0] for var in [P.source_duration,P.source_energy]),"Source fraction must not be None!"

  print("Outputting to " + outfile)

  update_period = min(10,P.fwhm/2.91)
  num_guess_steps = max(1000*(1+2*(P.pulse_shape == "square")),round(P.fwhm*50))

  # Creare plasma input file
  plasma_file = open(outfile, 'w')
  plasma_file.write("""// AC4DC input file generated from: """ + fname + """\n\n""")
  plasma_file.write("""#ATOMS\n""")
  for line in P.atoms:
    plasma_file.write(line)

  plasma_file.write("""\n#BOUND_FREE_EXCLUSIONS\n""")
  for line in P.scndry_exclusions:
    plasma_file.write(line)  
  
  plasma_file.write("""\n#VOLUME\n""")
  plasma_file.write("""%.2f      // Volume per molecule in Angstrom^3.\n""" % P.volume)
  plasma_file.write("""2000         // Radius of a sample in Angstrom. Used for effective escape rate of photo-electrons.\n""")
  plasma_file.write("""none         // Spatial boundary shape - options are none (confined system), spherical, cylindrical, planar.\n""")

  plasma_file.write("""\n#PULSE\n""")
  plasma_file.write("""%.0f         // Photon energy in eV.\n""" % P.energy)
  plasma_file.write("""%.1f         // Pulse-dependent definition of duration in femtoseconds: FWHM for Gaussian pulse, pulse width for square pulse.\n""" % P.fwhm)
  plasma_file.write(P.pulse_shape+"""     // Pulse shape.\n""")

  plasma_file.write("""\n#USE_COUNT\n""")
  plasma_file.write("""true         // enabled? true/false.\n""")
  plasma_file.write("""%.2f            // Photon density (x10^12 ph.µm-2)\n""" %P.photon_count)

  plasma_file.write("""\n#NUMERICAL\n""")
  plasma_file.write("""%.0f         // Initial guess for number of time step points for rate equation (Note e-e scattering part takes num_stiff_ministeps = 500 substeps for each step).\n""" %num_guess_steps)
  plasma_file.write("""26           // Number of threads in OpenMP.\n""")

  plasma_file.write("""\n#DYNAMIC_GRID\n""")
  plasma_file.write(P.grid_type+"""            // Grid regions preset, options are 'low', 'medium', 'high', (accuracy) among others (see README).\n""")
  plasma_file.write("""%.2f         // Grid update period in fs, (dynamic grid only).\n""" %update_period)

  if P.source_fraction not in [None,0] and P.source_energy not in [None,0]:
    plasma_file.write("""\n#ELECTRON_SOURCE\n""")
    plasma_file.write("""%.2f           // Additional electrons as fraction of other photoion. events.\n""" %P.source_fraction)
    plasma_file.write("""%.0f          // Electron energy in eV.\n""" %P.source_energy)
    plasma_file.write("""%.2f           // Fraction of pulse to emit over.\n""" %P.source_duration)
     

  plasma_file.write("""\n#OUTPUT\n""")
  plasma_file.write("""800          // Number of time steps in the output files.\n""")
  plasma_file.write("""4000         // Number of free-electron grid points in output file.\n""")

  plasma_file.write("""\n#DEBUG\n""")
  plasma_file.write("""999           // Simulation cutoff time in fs (can end before but not after).\n""")
  plasma_file.write("""0.01          // Interval to display current timestep [fs].\n""")
  plasma_file.write("""5             // Interval to update live plot (steps).\n""")

  plasma_file.write("####END####\n\n")
  
  plasma_file.write("Notes:")
  
  plasma_file.close()


# copied from elsewhere. TODO give this its own file or something and import only in the function to stop HPC form being upset. 
def find_mol_file_from_directory(input_directory, mol):
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

if __name__ == "__main__":
  main()