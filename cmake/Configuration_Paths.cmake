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

# Define a cmake variable for the MFEM source directory path. This will be used for all cmake files in the MoFA package.
set(MFEM_SOURCE_DIR "$ENV{HOME}/MFEM/mfem-4.8") # EDIT AS NECESSARY

# Define a cmake variable for the Gmsh source directory path. This will be used for all cmake files in the MoFA package.
set(GMSH_SOURCE_DIR "$ENV{HOME}/gmsh/gmsh-4.13.1-Linux64-sdk") # EDIT AS NECESSARY



# Define a cmake variable for the top directory of MoFA.
get_filename_component(cmake_SOURCE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(MOFA_SOURCE_DIR "${cmake_SOURCE_DIR}/..")

# Define a cmake variable for the Utility source directory path. This will be used for all cmake files in the MoFA package.
set(UTILITY_SOURCE_DIR "${MOFA_SOURCE_DIR}/Utility")

# Define a cmake variable for the Mesh_tools source directory path. This will be used for all cmake files in the MoFA package.
set(MESH_TOOLS_SOURCE_DIR "${MOFA_SOURCE_DIR}/Mesh_Tools")
