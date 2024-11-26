project('riemann-problem-solver', 'cpp', default_options : ['warning_level=3', 'cpp_std=c++23'])

riemann_solver_lib_include = include_directories('.')
riemann_solver_lib_deps = declare_dependency(include_directories : riemann_solver_lib_include)

riemann_solver_lib = static_library('riemann_solver', sources : ['riemann.cpp', 'riemann.hpp'], dependencies : [riemann_solver_lib_deps])
executable('riemann_solver_tests', 'tests.cpp', link_with : [riemann_solver_lib])
