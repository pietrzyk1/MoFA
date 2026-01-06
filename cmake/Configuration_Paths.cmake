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

# Set the minimum version of CMake that is allowed to be used for building MoFA
set(CMAKE_MIN_VERSION "3.22.1") # EDIT AS NECESSARY


# For HPC to use c++17 (uncomment if building MoFA on HPC)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)


# Define a cmake variable for the MFEM source directory path. This will be used for all cmake files in the MoFA package.
set(MFEM_SOURCE_DIR "$ENV{HOME}/MFEM/mfem-4.8") # EDIT AS NECESSARY

# Define a cmake variable for the Gmsh source directory path. This will be used for all cmake files in the MoFA package.
set(GMSH_SOURCE_DIR "$ENV{HOME}/gmsh/gmsh-4.13.1-Linux64-sdk") # EDIT AS NECESSARY


########################################################
# For Parallel Build #
########################################################

# Define which type of executables should be built.
option(BUILD_FOR_SERIAL "Build the executables for serial. If only the parallel build of MFEM is available (e.g., on an HPC environment), turn to OFF)" ON) # EDIT AS NECESSARY
option(BUILD_FOR_PARALLEL "Build the executables for MPI parallel." ON) # EDIT AS NECESSARY

# Define a cmake variable for the Parallel-built MFEM source directory path. This will be used for all cmake files of MPI builds in the MoFA package.
set(PAR_MFEM_SOURCE_DIR "$ENV{HOME}/MFEM/mfem-4.8_Par") # EDIT AS NECESSARY

# Define cmake variables for the Hypre include and lib directory paths. This will be used for all cmake files of MPI builds in the MoFA package.
set(HYPRE_INCLUDE_DIR "$ENV{HOME}/MFEM/hypre/src/hypre/include") # EDIT AS NECESSARY
set(HYPRE_LIB_DIR "$ENV{HOME}/MFEM/hypre/src/hypre/lib") # EDIT AS NECESSARY

# Define cmake variables for the Metis include and lib directory paths. This will be used for all cmake files of MPI builds in the MoFA package.
set(METIS_INCLUDE_DIR "$ENV{HOME}/MFEM/metis-5.1.0/include") # EDIT AS NECESSARY
set(METIS_LIB_DIR "$ENV{HOME}/MFEM/metis-5.1.0/lib") # EDIT AS NECESSARY


########################################################
# For Internal MoFA Structure #
########################################################

# Define a cmake variable for the top directory of MoFA.
get_filename_component(cmake_SOURCE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(MOFA_SOURCE_DIR "${cmake_SOURCE_DIR}/..")

# Define a cmake variable for the Utility source directory path. This will be used for all cmake files in the MoFA package.
set(UTILITY_SOURCE_DIR "${MOFA_SOURCE_DIR}/Utility")

# Define a cmake variable for the Mesh_tools source directory path. This will be used for all cmake files in the MoFA package.
set(MESH_TOOLS_SOURCE_DIR "${MOFA_SOURCE_DIR}/Mesh_Tools")
