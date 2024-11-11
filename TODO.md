# Cleanup for publication

2. Write some automated unit tests
3. Change config -- stop using hacky #defines
4. Get it compiling headless (appropriate use of #define)
1. Refactor to decouple atomic physics from rate equation solving. Make Hartree-Fock atomic code interpolate the rates for a given beam energy and store them in a file, to be read by the dynamics solver. Very important for capability in simulating heavier atoms (i.e. with more electron configurations).
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
16. Incorporate minimum and maximum energy into GridSpacing (perhaps rename it to GridParams)
18. Restructure parameter input and rate output files to use JSON format
19. GUI (Current candidate framework: Qt)
20. Optimise with static arrays - promote state_type to a N_FREE-dimensioned template for faster reads.
<<<<<<< HEAD
21. ODE integration routines: borrow from [rodent](https://www.github.com/jeanluct/rodent)'s ideas, make the function a template parameter rather than a virtual member. (Probably not limiting, but it's a fairly glaring misuse of virtual functions)
22. Upgrade dynamic grid algorithm to handle low-energy photoelectron peaks - currently mistakes them for MB peaks and so the solver fails.
=======
21. Upgrade dynamic grid algorithm to handle low-energy photoelectron peaks - currently mistakes them for MB peaks and so the solver fails.
>>>>>>> spencer/working_branch
