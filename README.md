# AC4DC Suite

A suite for simulating XFEL damage and its effect on image recovery in the context of ultrafast serial imaging (SFX/SPI) techniques. 

## AC4DC

The main code of the suite, AC4DC (the literal acronym is no longer apt), simulates the changes in electronic structure of atoms in the bulk of a micro-scale target or amorphous media as it is illuminated by an X-ray free-electron laser (XFEL) pulse. An independent atom approximation is adopted to represent a molecule as a collection of independent atoms, submerged into a bath of electrons that are formed as a result of ionization.

+ Computes cross-sections for atomic processes: photoionisation, fluorescence, Auger decay and electron-impact ionisation
+ Asynchronously solves for the time evolution of the non-equilibrium free electrons and the bound states of the atomic population, via the Boltzmann transport equation.
+ Allows for an arbitrary form of the free-electron distribution f(E) interpolated over a basis of B-splines. 
  + The basis is adaptive; it is transformed periodically to automatically assign more splines where sharp, non-polynomial peaks are present in the distribution. 

## Scatter

An auxilliary simulation, Scatter, generates scattering patterns off realistic targets constructed from PDB structure files, with the atoms’ states selected entirely based off the probability distributions produced by AC4DC, without regard for selections of prior snapshots. These are then compared with the scattering pattern produced by a structure in the 'ideal', undamaged case where no ionisation occurs. 

### Installing AC4DC

This project is built using the meson build system.

```
git clone https://github.com/Spuriosity1/AC4DC/tree/master && cd AC4DC
meson setup build
ninja -C build
```

Compatibility is only promised for Linux variants and macOS. (lack of Windows support is mainly due to the UNIX-style path assumptions used throughout. In principle, this may be corrected for by rewriting with `boost::filesystem`.) Previous versions of AC4DC have been tested in the following environments:
+ Debian 10 (buster), gcc8
+ macOS 10.14.6 (Mojave), gcc8, gcc9, gcc10 This project is built with `make`, but does not have a `configure` script. Some fiddling may be required depending on your OS.
The current version has been tested with:
+ WSl2 (Ubuntu), gcc10.
**Note for macOS users**  - Apple's standard `gcc` distribution does not support openmp, used for parallelisation. Code was compiled with `g++-10` from Homebrew. In principle, Homebrew llvm's `clang++` should also work, however its implementation of `openmp` seems to cause runtime segfaults. Attention is needed from a C++ expert to refactor this code to resolve this issue.

### Instructions for a fresh install for C++ beginners. 

Tested on an Ubuntu environment in WSL2. 

#### Compiling
0. Run `sudo apt-get install build-essential`
1. Install brew (e.g. linux: https://docs.brew.sh/Homebrew-on-Linux#requirements)
2. Run `brew install [formula]` on following formulae: `gcc@10`, `eigen`, `ncurses`
3. Ensure correct links for INC and LIB in this file. 
4. Run `make`

#### For live plotting (optional)
`brew install python3.9`  <br> 
Run 'pip3.9 install' for: pybind11, plotly, scipy.

#### Running and testing

`mv bin/ac4dc ac4dc`  <br>
`./ac4dc input/lys_example.mol`  <br>   
If the file cannot execute it may be a linking error. Run ‘ldd ac4dc’ to check dependencies.

### Dependencies

Required:
+ Eigen 3 (vector algebra library, http://eigen.tuxfamily.org/index.php?title=Main_Page)
+ Wigner Symbols (3-j and 6-j symbols, https://github.com/joeydumont/wignerSymbols)
+ OpenMP (multi-threading, https://www.openmp.org/)
+ ncurses (text-based UI, https://invisible-island.net/ncurses/)  

The free-electron distribution may be plotted `live’ via regular writes to AC4DC/_live_plot.png. This requires python3.9, with the packages:
- Pybind11
- Plotly
- Scipy
This and the use of python may be disabled by commenting out `#define PYBIND` in include/config.h 



### Configuration

AC4DC reads the composition of the target, the pulse parameters, and various hyperparameters (e.g pertaining to the spline knot grid) from the molecular (.mol) file it is provided. See AC4DC/input/mol_template.mol for the style of these files and the parameters available.     

In AC4DC/include/config.h various features of the simulation can be disabled (e.g. plasma processes, live plotting, backing up of data). Most features are enabled by default, with the exception of tracking the electron cascades seeded by each element ("TRACK_SINGLE_CONTINUUM"), as this is computationally costly. 

### Grid regions preset

Relevant files:
  include/GridSpacing.hpp
  src/DynamicRegions.cpp

        case 'd':
            preset.selected = DynamicGridPreset::dismal_acc;
            break;
        case 'l':
            preset.selected = DynamicGridPreset::low_acc;
            break;
        case 'm':
            preset.selected = DynamicGridPreset::medium_acc;
            break;
        case 'h':
            preset.selected = DynamicGridPreset::high_acc;
            break;    
        case 'n':
            preset.selected = DynamicGridPreset::no_dirac;
            break;         
        case 't':
            preset.selected = DynamicGridPreset::training_wheels;
            break;                      
        case 'A':
            preset.selected = DynamicGridPreset::heavy_support;
            break;  
        case 'B':
            preset.selected = DynamicGridPreset::Zr_support;
            break;              
        case 'D':
            preset.selected = DynamicGridPreset::lower_dirac_support;
            break;            

### Grid regions preset

Relevant files:
  include/GridSpacing.hpp
  src/DynamicRegions.cpp

        case 'd':
            preset.selected = DynamicGridPreset::dismal_acc;
            break;
        case 'l':
            preset.selected = DynamicGridPreset::low_acc;
            break;
        case 'm':
            preset.selected = DynamicGridPreset::medium_acc;
            break;
        case 'h':
            preset.selected = DynamicGridPreset::high_acc;
            break;    
        case 'n':
            preset.selected = DynamicGridPreset::no_dirac;
            break;         
        case 't':
            preset.selected = DynamicGridPreset::training_wheels;
            break;                      
        case 'A':
            preset.selected = DynamicGridPreset::heavy_support;
            break;  
        case 'B':
            preset.selected = DynamicGridPreset::Zr_support;
            break;              
        case 'D':
            preset.selected = DynamicGridPreset::lower_dirac_support;
            break;            

### Running AC4DC and workflow

0. Compile  <br>
`make`  <br>
`mv bin/ac4dc ac4dc`  <br>
1. Simulate the damage  <br>
Run with `./ac4dc input/lys_example.mol`  <br>
Once completed, the simulation produces a folder containing the data in AC4DC/output/__Molecular/lys_example_#, where # is the number for the simulation (this will be 1 if it is the first time the program has been run with this .mol file). 
The simulation saves the data to AC4DC/output/backup_data hourly. This can be loaded if specified in the input file. Note that old files are never removed, but they are overwritten.  <br> 
2. Analysis  <br>
The plotting/scattering programs always take in a folder handle, which is searched for within the AC4DC/output/__Molecular directory and subdirectories (ambiguities are flagged and handled).  <br>
View an interactive plot for the electron distribution with  <br>
`python3.9 scripts/_generate_interactive lys_example_#`  <br>
and generate various plots with  <br>
`python3.9 scripts/generate_plots lys_example_#`  <br>
These will be contained in AC4DC/output/graphs/

Automated generation of batches of simulations is also supported.  <br>
Configure generate_batch.py and run:  <br>
`python3.9 generate_batch.py input/templates/lys` will generate a batch of simulations with each permutation of parameters, located in the folder AC4DC/input/batch_lys/.  <br>
run_sims.sh can be configured to run these batches.  <br>
If these results are manually moved to the folder AC4DC/output/__Molecular/my_batch/, we can run  <br>
`python3.9 scripts/_generate_interactive my_batch`  <br>
and  <br>
`python3.9 scripts/generate_plots my_batch`  <br>
which will generate figures for each simulation, saved in AC4DC/output/graphs/my_batch  <br>



### Key Approximations in AC4DC:

1. All atoms interact with X-rays independent of each other.
2. At any time, an atom is described by an average-over-configuration state.
3. Time evolution of atoms is described by a system of coupled rate equations, in which perturbation theory calculations are accurate.
4. There are enough free electrons per Debye sphere for a Boltzmann treatment to be appropriate.
5. Electron plasma is isotropic and homogeneous in position space; isotropic in momentum space.
6. Electrons generated by atomic processes are instantaneously spread in energy space to a width of order 10 eV.* 
Rate coefficients can be approximated with a monochromatic source.
The photon density is homogeneous at each point in time.

*At present, low-energy photoelectrons are emitted with line profile widths of order 10 eV, which is a rough approximation for the true profile, but this is more of an untreated edge case than a key approximation.


### Output format

Running the program generates an `output` folder in the root directory.
In the Xsections folder, the program stores a series of data files for Auger, Fluorescence and Photoionisation cross sections in the as `.txt`, in the format
`[ rate value ] [ index to ] [ index from ][ energy of transition ]`
Corresponding lines in the `index.txt` file explain which configuration the to/from indices correspond to.

EII parameters are stored in "sort-of-json" format - please note that the program does NOT use a robust json parser, and is sensitive to line formatting.



### Bibliography:

+ A. Kozlov and H. M. Quiney, _Comparison of Hartree-Fock and local density exchange approximations for calculation of radiation damage dynamics of light and heavy atoms in the field of x-ray free electron laser_, Phys. Scripta **94**, 075404 (2019). DOI: [10.1088/1402-4896/ab097c](https://doi.org/10.1088/1402-4896/ab097c)
+ C. P. Bhalla, N. O. Folland, and M. A. Hein, _Theoretical K-Shell Auger Rates, Transition Energies, and Fluorescence Yields for Multiply Ionized Neon_, Phys. Rev. A **8**, 649 (1973). DOI: [10.1103/PhysRevA.8.649](https://doi.org/10.1103/PhysRevA.8.649)
+ O. Yu. Gorobtsov, U. Lorenz, N. M. Kabachnik, and I. A. Vartanyants, _Theoretical study of electronic damage in single-particle imaging experiments at x-ray free-electron lasers for pulse durations from 0.1 to 10 fs_, Phys. Rev. E **91**, 062712 (2015). DOI: [10.1103/PhysRevE.91.062712](https://doi.org/10.1103/PhysRevE.91.062712)
+ W. R. Johnson, _Atomic Structure Theory: Lectures on Atomic Physics_, (Springer, Berlin, Heidelberg, 2007). DOI: [10.1007/978-3-540-68013-0](https://doi.org/10.1007/978-3-540-68013-0)
+ Leonov, A., D. Ksenzov, A. Benediktovitch, I. Feranchuk, and U. Pietsch, _Time Dependence of X-Ray Polarizability of a Crystal Induced by an Intense Femtosecond X-Ray Pulse._ IUCrJ **1**, no. 6 (2014): 402–17. DOI: [10.1107/S2052252514018156](https://doi.org/10.1107/S2052252514018156)
