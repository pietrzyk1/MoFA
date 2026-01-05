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
project(Porescale_Avger)

# Set the executable to be created in the bin directory of the project source directory (this line must be placed after defining the project)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Tell CMake the directory containing the project's header files (i.e., the .h files). This is often in the /include directory.
include_directories(${UTILITY_SOURCE_DIR}/include)


# Add the executable we want to build. The first argument will be the name of the produced executable..
add_executable(porescale_avger ${PROJECT_SOURCE_DIR}/src/Porescale_Avger.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${UTILITY_SOURCE_DIR}/src/mfem_util.cpp)

# Tell CMake any directories containing the project's library files (i.e., the .a, .so, or .lib files). This is often the ./lib directories.
target_link_directories(porescale_avger PRIVATE ${MFEM_SOURCE_DIR})

# Tell CMake any directories containing the project's header files (i.e., the .h files). This is often the ./include directories.
target_include_directories(porescale_avger PRIVATE ${MFEM_SOURCE_DIR})

# Tell CMake any libraries to link to the project's executable. If you told CMake the directories containing the project's libraries in "link_directories", you can simply put the library names in the second argument. Otherwise, put the full path to the library file (i.e., the .a, .so, or .lib file) in the second argument.
target_link_libraries(porescale_avger PRIVATE libmfem.a)
