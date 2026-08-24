cmake_minimum_required(VERSION 3.15...3.29)

project(wickd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build'. Running
  it directly will almost certainly not produce the desired result. If
  you are a user trying to install this package, please use the command
  below, which will install all necessary build dependencies, compile
  the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to re-run the above
  after editing C++ files.")
endif()

# Try to import all Python components potentially needed by nanobind
find_package(Python 3.10
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

# Look for the Boost libraries
cmake_policy(SET CMP0167 NEW)
find_package(Boost)

set(SRC_LIST
  src/algebra/equation.cc
  src/algebra/expression.cc
  src/algebra/expression_normal_order.cc
  src/algebra/index.cc
  src/algebra/sqoperator.cc
  src/algebra/symbolic_term.cc
  src/algebra/symbolic_term_normal_order.cc
  src/algebra/tensor.cc
  src/algebra/term.cc
  src/api/api_wicked.cc
  src/api/combinatorics_api.cc
  src/api/equation_api.cc
  src/api/expression_api.cc
  src/api/index_api.cc
  src/api/operator_api.cc
  src/api/orbital_space_api.cc
  src/api/rational_api.cc
  src/api/sqoperator_api.cc
  src/api/tensor_api.cc
  src/api/term_api.cc
  src/api/wick_theorem.cc
  src/diagrams/contraction.cc
  src/diagrams/graph_matrix.cc
  src/diagrams/operator.cc
  src/diagrams/operator_expression.cc
  src/diagrams/operator_product.cc
  src/diagrams/wick_theorem.cc
  src/diagrams/wick_theorem_canonicalize_graph.cc
  src/diagrams/wick_theorem_composite_contractions.cc
  src/diagrams/wick_theorem_elementary_contractions.cc
  src/diagrams/wick_theorem_helpers.cc
  src/diagrams/wick_theorem_process_contractions.cc
  src/helpers/combinatorics.cc
  src/helpers/helpers.cc
  src/helpers/orbital_space.cc
  src/helpers/rational.cc
)

# We are now ready to compile the actual extension module
nanobind_add_module(
  # Name of the extension
  _wickd

  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  STABLE_ABI

  # Build libnanobind statically and merge it into the
  # extension (which itself remains a shared library)
  #
  # If your project builds multiple extensions, you can
  # replace this flag by NB_SHARED to conserve space by
  # reusing a shared libnanobind across libraries
  NB_STATIC

  # Source code goes here
  ${SRC_LIST}
)

target_include_directories(_wickd PRIVATE src src/algebra src/diagrams)

# Install directive for scikit-build-core
install(TARGETS _wickd LIBRARY DESTINATION wickd)

# Check if Boost was found
if(Boost_FOUND)
  option(USE_BOOST_1024_INT "Use Boost multiprecision 1024-bit integers" ON)
else()
  option(USE_BOOST_1024_INT "Use Boost multiprecision 1024-bit integers" OFF)
endif()

if(USE_BOOST_1024_INT AND NOT Boost_FOUND)
  message(FATAL_ERROR "USE_BOOST_1024_INT=ON but Boost was not found")
endif()

if(USE_BOOST_1024_INT)
  message(STATUS "Using Boost 1024-bit integers (disable with -DUSE_BOOST_1024_INT=OFF)")
  target_compile_definitions(_wickd PRIVATE USE_BOOST_1024_INT)
  target_include_directories(_wickd SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
else()
  message(STATUS "Not using Boost 1024-bit integers (enable with -DUSE_BOOST_1024_INT=ON)")
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)

if(CODE_COVERAGE)
  message(STATUS "Code coverage enabled")
  target_compile_options(_wickd PRIVATE --coverage)
  target_link_options(_wickd PRIVATE --coverage)
endif()

message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")
