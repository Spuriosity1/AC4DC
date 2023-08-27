
main_sources = files(
'BoundTransport.cpp',
'ComputeRateParam.cpp',
'Constant.cpp',
'DecayRates.cpp',
'Dipole.cpp',
'DynamicRegions.cpp',
'EigenSolver.cpp',
'ElectronRateSolver.cpp',
'ElectronSolverGridDynamics.cpp',
'ElectronSolverIO.cpp',
'FreeDistribution.cpp',
'Grid.cpp',
'HartreeFock.cpp',
'Input.cpp',
'MolInp.cpp',
'Numerics.cpp',
'Potential.cpp',
'Pulse.cpp',
'RadialWF.cpp',
'RateSystem.cpp',
'SplineBasis.cpp',
'SplineIntegral.cpp'
)

# sources requiring linking with Python
py_sources = files( 'Display.cpp','Plotting.cpp')

ac4dc_main_source = files('main.cpp')


project_sources += main_sources

subdir('Wigner')
