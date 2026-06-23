# lL.lH.
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
# lL.lH.

# Include the path configuration cmake file for the various source directory paths
include(${CMAKE_SOURCE_DIR}/../cmake/Configuration_Paths.cmake)

# Tell CMake the minimum version of CMake that it's allowed to use for building the project
cmake_minimum_required(VERSION ${CMAKE_MIN_VERSION})


# Tell CMake the project name and version to be built
project(Mesh_Generation)


# Include the ./include directory for the project
include_directories(${GMSH_SOURCE_DIR}/include)
include_directories(${UTILITY_SOURCE_DIR}/include)
include_directories(${MESH_TOOLS_SOURCE_DIR}/include)

# Link the libraries in the ./lib directory
link_directories(${GMSH_SOURCE_DIR}/lib)


# Set the executable to be created in the bin directory of the project source directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)


#######################################
# Make the generate mesh executable
#######################################

# Add the executable we want to build
add_executable(generate_mesh ${PROJECT_SOURCE_DIR}/src/generate_mesh.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${MESH_TOOLS_SOURCE_DIR}/src/gmsh_utils.cpp ${MESH_TOOLS_SOURCE_DIR}/src/gmsh_import_STL.cpp)

# Now we actually link the libraries to the executable
target_link_libraries(generate_mesh libgmsh.so)


#######################################
# Make the scale mesh executable
#######################################

# Add the executable we want to build
add_executable(scale_mesh ${PROJECT_SOURCE_DIR}/src/scale_mesh.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${MESH_TOOLS_SOURCE_DIR}/src/gmsh_utils.cpp)

# Now we actually link the libraries to the executable
target_link_libraries(scale_mesh libgmsh.so)


#######################################
# Make the juxtapose unit cell executable
#######################################

# Add the executable we want to build
add_executable(juxtapose_unitcell ${PROJECT_SOURCE_DIR}/src/juxtapose_unitcell.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${MESH_TOOLS_SOURCE_DIR}/src/gmsh_utils.cpp)

# Now we actually link the libraries to the executable
target_link_libraries(juxtapose_unitcell libgmsh.so)


#######################################
# Make the light-weight generate mesh from STL executable
#######################################
if(EXISTS ${PROJECT_SOURCE_DIR}/src/generate_STL_mesh.cpp)

    # Add the executable we want to build. The first argument will be the name of the produced executable.
    add_executable(generate_STL_mesh ${PROJECT_SOURCE_DIR}/src/generate_STL_mesh.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${MESH_TOOLS_SOURCE_DIR}/src/gmsh_utils.cpp ${MESH_TOOLS_SOURCE_DIR}/src/gmsh_import_STL.cpp)

    # Tell CMake any libraries to link to the project's executable. If you told CMake the directories containing the project's libraries in "link_directories", you can simply put the library names in the second argument. Otherwise, put the full path to the library file (i.e., the .a, .so, or .lib file) in the second argument.
    target_link_libraries(generate_STL_mesh libgmsh.so)

endif()