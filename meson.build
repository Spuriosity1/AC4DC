project('ac4dc','cpp',
  version:'0.2.0',
  license: ['GPLv3'],
  meson_version: '>=1.1.0',
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

ac4dc_bin_deps = [
  dependency('openmp',include_type:'system'),
  dependency('eigen3',include_type:'system')
  ]

ac4dc_bin_libdeps = [
  wigner_lib
]
# the main executable
ac4dc_bin = executable('ac4dc', 
  [main_sources, ac4dc_main_source],
  include_directories: inc,
  dependencies: ac4dc_bin_deps,
  link_with: ac4dc_bin_libdeps,
)

### TESTS
if get_option('build_tests')
  subdir('tests')
  ode_solver = executable('ode_solver', ode_solver_src,include_directories: inc)
  test('Adams-Basforth-Moulton ODE', ode_solver, args: ['abm', 'sho'])
  test(' ODE', ode_solver, args: ['abm', 'sho'])


endif


# Python linking (to fix later)
if get_option('link_python')

  py3_inst = import('python').find_installation('python3')
  
  pybind11_config = find_program('pybind11-config')
  pybind11_config_ret = run_command(pybind11_config, ['--includes'],check: true)
  pybind11 = declare_dependency(
      include_directories: [pybind11_config_ret.stdout().split('-I')[-1].strip()],
  )
  
  python3 = dependency('python3')
  py3_inst = import('python').find_installation('python3')

  ac4dc_lib = py3_inst.extension_module('ac4dc_py', 
  main_sources,
  include_directories: inc,
  dependencies: [ac4dc_bin_deps, dependency('ncurses'),pybind11, python3],
  link_with: [ac4dc_bin_libdeps],
  link_language: 'cpp',
  override_options: [
        'cpp_rtti=true',
    ],
  cpp_args: '-DNCURSES'
  )
endif
