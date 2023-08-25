project('ac4dc','cpp',
  version:'0.2.0',
  license: ['GPLv3'],
  meson_version: '>=0.63.0',
  default_options: ['warning_level=3', 'cpp_std=c++17']
  )

library_path = []
project_sources = []
project_header_files = []
project_test_sources = []
project_benchmark_sources = []

inc = [include_directories('include', 'include/Wigner')]

# check if we are on homebrew
#OSX_ADHOC =  (run_command('brew','-v', check:true).returncode() == 0)I
  
subdir('include')
# add more include directories here...

subdir('src')

# Dependencies
# pybind11
py3_inst = import('python').find_installation('python3')

pybind11_config = find_program('pybind11-config')
pybind11_config_ret = run_command(pybind11_config, ['--includes'],check: true)
pybind11 = declare_dependency(
    include_directories: [pybind11_config_ret.stdout().split('-I')[-1].strip()],
)

python3 = dependency('python3')
py3_inst = import('python').find_installation('python3')

ac4dc_bin_deps = [
  dependency('ncurses'),
  dependency('openmp'),
  dependency('eigen3')
  ]


ac4dc_bin_libdeps = [
  wigner_lib
]


ac4dc_lib = py3_inst.extension_module('ac4dc', 
main_sources,
include_directories: inc,
dependencies: [ac4dc_bin_deps, pybind11, python3],
link_with: [ac4dc_bin_libdeps],
link_language: 'cpp',
override_options: [
        'cpp_rtti=true',
    ],
cpp_args: '-DNCURSES'
)
ac4dc_bin = executable('ac4dc_headless', 
  main_sources,
  include_directories: inc,
  dependencies: ac4dc_bin_deps,
  link_with: ac4dc_bin_libdeps,
)
