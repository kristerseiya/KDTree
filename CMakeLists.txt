cmake_minimum_required(VERSION 3.17)
project(kdtree)

set(CMAKE_CXX_STANDARD 11)

# OPTIONS
############################
option(BUILD_KDTREE_PYTHON "build kdtree python interface" ON)
option(BUILD_KDTREE_TEST "build kdtree library test app" ON)
############################

set(CMAKE_CXX_FLAGS -O3)

# EDIT PATHS ACCORDINGLY
############################
if (NOT DEFINED PCDIO_PYTHON_INSTALL_PATH)
set(KDTREE_PYTHON_INSTALL_PATH
  /Users/Krister/anaconda3/envs/python39/lib/python3.9/site-packages)
endif()

if (NOT DEFINED pybind11_FOUND)
find_package(pybind11 REQUIRED HINTS /usr/local/opt/pybind11)
endif()
############################


file(GLOB KDTREE_SOURCES source/kdtree.cpp)
add_library(kdtree STATIC ${KDTREE_SOURCES})
# set_target_properties(kdtree PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(kdtree PRIVATE include)

if (BUILD_KDTREE_TEST)
    add_subdirectory(test)
endif()

if (BUILD_KDTREE_PYTHON)
pybind11_add_module(pykdtree ${KDTREE_SOURCES} python/kdtree_pybind.cpp)
set_target_properties(pykdtree PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(pykdtree PRIVATE include)
install(TARGETS pykdtree DESTINATION ${KDTREE_PYTHON_INSTALL_PATH})
endif()
