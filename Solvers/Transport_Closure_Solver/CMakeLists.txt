#lL.lH.
#
# Copyright (c) 2025, Lawrence Livermore National Security, LLC
# and other MoFA project developers. All Rights reserved.
# See files LICENSE and NOTICE for details. LLNL-CODE-2006961.
#
# This file is part of the MoFA Project. For more information
# and source code availability visit https://github.com/pietrzyk1/MoFA.
#
# SPDX-License-Identifier: BSD-3-Clause
#
#lL.lH.

# Include the path configuration cmake file for the various source directory paths
include(${CMAKE_SOURCE_DIR}/../../cmake/Configuration_Paths.cmake)

# Tell CMake the minimum version of CMake that it's allowed to use for building the project
cmake_minimum_required(VERSION ${CMAKE_MIN_VERSION})


# Tell CMake the project name and version to be built
project(Transport_Closure_Solver)


# Tell CMake the directory containing the project's header files (i.e., the .h files). This is often in the /include directory.
include_directories(${MFEM_SOURCE_DIR})
include_directories(${UTILITY_SOURCE_DIR}/include)

# Tell CMake the directory containing the project's library files (i.e., the .a, .so, or .lib files). This is often the ./lib directory.
link_directories(${MFEM_SOURCE_DIR})


# Set the executable to be created in the bin directory of the project source directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Add the executable we want to build. The first argument will be the name of the produced executable..
add_executable(transport_solver ${PROJECT_SOURCE_DIR}/src/Transport_Closure_Solver.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${UTILITY_SOURCE_DIR}/src/mfem_util.cpp)

# Tell CMake to link the project's libraries to the executable. If you told CMake the directories containing the project's libraries in "link_directories", you can simply put the library names in the second argument. Otherwise, put the full path to the library file (i.e., the .a, .so, or .lib file) in the second argument.
target_link_libraries(transport_solver libmfem.a)


# If the OpenMP library can be found, create an executable for solving the closure problems as an ensemble using multi-threading parallelism.
find_package(OpenMP)
if(OpenMP_FOUND)
    # Define USE_OPENMP, which can be used to trigger OpenMP tasks in the executable
    add_definitions(-DUSE_OPENMP) # Tells C++ preprocessor to define the macro "USE_OPENMP"

    # Add the executable for solving closure problems as an ensemble
    add_executable(transport_solver_ensemble ${PROJECT_SOURCE_DIR}/src/Transport_Closure_Solver_Ensemble.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${UTILITY_SOURCE_DIR}/src/mfem_util.cpp)

    # Tell CMake to link the project's libraries to the executable (this includes the link to the OpenMP library).
    target_link_libraries(transport_solver_ensemble libmfem.a) 
    target_link_libraries(transport_solver_ensemble OpenMP::OpenMP_CXX) # For the OpenMP library
endif()
