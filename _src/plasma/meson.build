plasma_sources = files([
    'DynamicRegions.cpp',
    'ElectronRateSolver.cpp',
    'ElectronSolverGridDynamics.cpp',
    'ElectronSolverIO.cpp',
    'SplineIntegral.cpp',
    'struct/SplineBasis.cpp',
    'struct/FreeDistribution.cpp',
    'struct/Pulse.cpp',
    'struct/RateSystem.cpp',
    'struct/MolInp.cpp'
])

plasma_lib = static_library('plasma', plasma_sources, 
    dependencies: global_deps,
    include_directories: [plasma_include, common_include, atomic_include],
    link_with: [common_lib, atomic_lib],
    cpp_args: main_args)
