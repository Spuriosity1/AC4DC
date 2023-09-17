project('ac3dc','cpp',
  version:'0.2.0',
  license: ['GPLv3'],
  meson_version: '>=1.1.0',
  default_options: ['cpp_std=c++17', 'warning_level=1']
  )

main_args = []

if get_option('disable_tbr')
  main_args += '-DNO_TBR'
endif
if get_option('disable_ee')
  main_args += '-DNO_EE'
endif
if get_option('disable_eii')
  main_args += 'NO_EII'
endif
if not get_option('ministep_updating')
  main_args += '-DNO_MINISTEP_UPDATING'
endif

library_path = []
project_sources = []
project_header_files = []
project_test_sources = []
project_benchmark_sources = []


# get the includes
subdir('_include')
project_include = [
  include_directories('_include'),
  atomic_include, plasma_include, display_include, wigner_include, 
  common_include] 


if get_option('interactive')
  # pybind11
  py3_inst = import('python').find_installation('python3')

  pybind11_config = find_program('pybind11-config')
  pybind11_config_ret = run_command(pybind11_config, ['--includes'],check: true)
  pybind11 = declare_dependency(
      include_directories: [pybind11_config_ret.stdout().split('-I')[-1].strip()],
  )

  python3 = dependency('python3')
  py3_inst = import('python').find_installation('python3')
endif


# Dependency management
json_dep = dependency('nlohmann_json')
# json_schema_dep = dependency('nlohmann_json_schema_validator')
omp_dep = dependency('openmp')
eigen3_dep = dependency('eigen3')

global_deps = [
  json_dep,
  # json_schema_dep,
  omp_dep,
  eigen3_dep
  ]

ac4dc_bin_deps = global_deps

if get_option('interactive')
  ac4dc_bin_deps += [dependency('ncurses')]
  main_args += '-DNCURSES'
endif


subdir('_src')
subdir('_tests')

ac4dc_bin_libdeps = [
  atomic_lib,
  # plasma_lib,
  # display_lib,
  wigner_lib,
  common_lib
]

# ac4dc_bin = executable('ac4dc', 
#   main_sources,
#   include_directories: project_include,
#   dependencies: ac4dc_bin_deps,
#   link_with: ac4dc_bin_libdeps,
#   cpp_args: main_args
# )

ratecalc_bin = executable('ratecalc', 
  ratecalc_sources,
  include_directories: [common_include, atomic_include],
  dependencies: global_deps,
  link_with: [common_lib, atomic_lib],
  cpp_args: main_args
)


