atomic_sources = files( [
    'BoundTransport.cpp',
    'ComputeRateParam.cpp',
    'DecayRates.cpp',
    'HartreeFock.cpp',
    'Potential.cpp',
    'numerics/EigenSolver.cpp',
    'numerics/NumericalIntegral.cpp',
    'struct/Grid.cpp',
    'struct/RadialWF.cpp',
    'struct/RateData.cpp',
    'struct/HFInputParam.cpp'
    ])

atomic_lib = static_library('atomic', atomic_sources, 
    dependencies: global_deps, 
    include_directories: [atomic_include, wigner_include],
    link_with: wigner_lib,
    cpp_args: main_args)
