common_sources = files( [
  'Dipole.cpp',
  'RateData.cpp'
  ])

common_lib = static_library('common', common_sources, 
  dependencies: json_dep,
  include_directories: common_include,
  cpp_args: main_args)

