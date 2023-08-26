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
`brew install python3.9` 
Run 'pip3.9 install' for: pybind11, plotly, scipy.

#### Running and testing

`mv bin/ac4dc ac4dc`
`./ac4dc input/lys_example.mol`   
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

AC4DC reads the composition of the target, the pulse parameters, and various hyperparameters (e.g pertaining to the spline knot grid)  from the molecular (.mol) file it is provided. See AC4DC/input/mol_template.mol for the style of these files and the parameters available.     

In AC4DC/include/config.h various features of the simulation can be disabled (e.g. plasma processes, live plotting, backing up of data).

### Running AC4DC and workflow

0. Compile
`make`
`mv bin/ac4dc ac4dc`
1. Simulate the damage
Run with `./ac4dc input/lys_example.mol`
Once completed, the simulation produces a folder containing the data in AC4DC/output/__Molecular/lys_example_#, where # is the number for the simulation (this will be 1 if it is the first time the program has been run with this .mol file). 
The simulation saves the data to AC4DC/output/backup_data hourly. This can be loaded if specified in the input file. Note that old files are never removed, but they are overwritten. 
2. Analysis
The plotting/scattering programs always take in a folder handle, which is searched for within the AC4DC/output/__Molecular directory and subdirectories (ambiguities are flagged and handled).
View an interactive plot for the electron distribution with
`python3.9 scripts/_generate_interactive lys_example_#`
and generate various plots with
`python3.9 scripts/generate_plots lys_example_#`
These will be contained in AC4DC/output/graphs/

Automated generation of batches of simulations is also supported.
Configure generate_batch.py and run:
`python3.9 generate_batch.py input/templates/lys` will generate a batch of simulations with each permutation of parameters, located in the folder AC4DC/input/batch_lys/.
run_sims.sh can be configured to run these batches. 
If these results are manually moved to the folder AC4DC/output/__Molecular/my_batch/, we can run
`python3.9 scripts/_generate_interactive my_batch`
and
`python3.9 scripts/generate_plots my_batch`
which will generate figures for each simulation, saved in AC4DC/output/graphs/my_batch



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


### Known Bugs

- When compiled under gcc8 on Debian, a malloc error is thrown when trying to run lysozyme. The error has something to do with the allocation of the Q tensors, but the source is not clear. Unknown whether this is still an issue.
- The fit for the free electron distribution gets 'stuck' when attempting to simulate lysozyme.Gd, where the photon energy is ~100 eV above the L-shell edge calculated by AC4DC (under the orbital-shell approximation), e.g. input/galli/lys_galli_HF_L_edge.

## TODO

1.  Refactor to decouple atomic physics from rate equation solving. Make Hartree-Fock atomic code interpolate the rates for a given beam energy and store them in a file, to be read by the dynamics solver. Very important for capability in simulating heavier atoms (i.e. with more electron configurations).
2. Fermi-sea collision kernel
3. Finish implementation of electron filtration by the water background.
4. Finish implementation of correction for bound transport.
5. Improve pulse shape code quality, add capability to input arbitrary temporal profiles, or mimic the stochastic temporal profile of a SASE pulse.
6. Implement stochastic spectral profile for low-energy (< 1000 eV) photoelectrons (current approximation is sufficient for high energies), or at least widen low-energy photoelectron emission profiles to a more realistic width.
7. Smoothing to stabilize fitting and grid updating when running with coarser grids.
8. The cutoff/transition energy (for an electron to be considered thermalised) currently only updates when the grid updates. It should instead update independently and frequently. 
9. Refactor to get rid of compiler warnings wherever possible (usually about signed comparisons)
10. Fix bug where `-s` flag causes crash
11. Move to a proper database system to store input/output data
12. Implement 'output version control' for atomic parameters in storage: avoid unnecessary recalculation, guarantee recalculation if new input parameters are incompatible
13. Add methods to `Input.cpp` to enable reading/writing salient parameters to file, e.g. `output/C/run_2021-04-11/input.txt`
14. Add linear search implementation to input logic
16. Incorporate minimum and maximum energy into GridSpacing (perhaps rename it to GridParams)
17. Cmake build system
18. Restructure parameter input and rate output files to use JSON format
19. GUI (Curent candidate framework: Qt)
20. Optimise with static arrays - promote state_type to a N_FREE-dimensioned template for faster reads.

### Bibliography:

+ A. Kozlov and H. M. Quiney, _Comparison of Hartree-Fock and local density exchange approximations for calculation of radiation damage dynamics of light and heavy atoms in the field of x-ray free electron laser_, Phys. Scripta **94**, 075404 (2019). DOI: [10.1088/1402-4896/ab097c](https://doi.org/10.1088/1402-4896/ab097c)
+ C. P. Bhalla, N. O. Folland, and M. A. Hein, _Theoretical K-Shell Auger Rates, Transition Energies, and Fluorescence Yields for Multiply Ionized Neon_, Phys. Rev. A **8**, 649 (1973). DOI: [10.1103/PhysRevA.8.649](https://doi.org/10.1103/PhysRevA.8.649)
+ O. Yu. Gorobtsov, U. Lorenz, N. M. Kabachnik, and I. A. Vartanyants, _Theoretical study of electronic damage in single-particle imaging experiments at x-ray free-electron lasers for pulse durations from 0.1 to 10 fs_, Phys. Rev. E **91**, 062712 (2015). DOI: [10.1103/PhysRevE.91.062712](https://doi.org/10.1103/PhysRevE.91.062712)
+ W. R. Johnson, _Atomic Structure Theory: Lectures on Atomic Physics_, (Springer, Berlin, Heidelberg, 2007). DOI: [10.1007/978-3-540-68013-0](https://doi.org/10.1007/978-3-540-68013-0)
+ Leonov, A., D. Ksenzov, A. Benediktovitch, I. Feranchuk, and U. Pietsch, _Time Dependence of X-Ray Polarizability of a Crystal Induced by an Intense Femtosecond X-Ray Pulse._ IUCrJ **1**, no. 6 (2014): 402–17. DOI: [10.1107/S2052252514018156](https://doi.org/10.1107/S2052252514018156)
