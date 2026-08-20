
# Sources shared by both binaries (atomic input parsing, radial/potential/grid
# machinery, numerics, and the HDF5 rate serialization seam).
common_sources = files(
'Constant.cpp',
'Dipole.cpp',
'EigenSolver.cpp',
'Grid.cpp',
'Input.cpp',
'Numerics.cpp',
'Potential.cpp',
'RadialWF.cpp',
'RateHDF5.cpp',
)

# Stage 1 (Hartree-Fock atomic physics) — linked only into atomic_rate_data.
hf_sources = files(
'ComputeRateParam.cpp',
'DecayRates.cpp',
'HartreeFock.cpp',
)

# Stage 2 (plasma dynamics / Boltzmann solver) — linked only into ac4dc.
solver_sources = files(
'BoundTransport.cpp',
'Display.cpp',
'DynamicRegions.cpp',
'ElectronRateSolver.cpp',
'ElectronSolverGridDynamics.cpp',
'ElectronSolverIO.cpp',
'FreeDistribution.cpp',
'MolInp.cpp',
'Plotting.cpp',
'Pulse.cpp',
'RateSystem.cpp',
'SplineBasis.cpp',
'SplineIntegral.cpp',
)

# Entry points.
atomic_rate_data_main = files('atomic_rate_data.cpp')
ac4dc_main = files('main.cpp')

project_sources += common_sources
project_sources += hf_sources
project_sources += solver_sources

subdir('Wigner')
