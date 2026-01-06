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
project(Transport_Solver)

# Set the executable to be created in the bin directory of the project source directory (this line must be placed after defining the project)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Tell CMake any global directories containing the project's header files (i.e., the .h files). This is often in the /include directory.
include_directories(${UTILITY_SOURCE_DIR}/include)


################
# Serial build #
################

# If BUILD_FOR_SERIAL is 'ON' (check cmake/Configuration_Paths.cmake), then build the solver assuming the serial build of MFEM.
if(BUILD_FOR_SERIAL)

    # Add the executable we want to build. The first argument will be the name of the produced executable.
    add_executable(transport_solver ${PROJECT_SOURCE_DIR}/src/Transport_Solver.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${UTILITY_SOURCE_DIR}/src/mfem_util.cpp)

    # Tell CMake any directories containing the project's library files (i.e., the .a, .so, or .lib files). This is often the ./lib directories.
    target_link_directories(transport_solver PRIVATE ${MFEM_SOURCE_DIR})

    # Tell CMake any directories containing the project's header files (i.e., the .h files). This is often the ./include directories.
    target_include_directories(transport_solver PRIVATE ${MFEM_SOURCE_DIR})

    # Tell CMake any libraries to link to the project's executable. If you told CMake the directories containing the project's libraries in "link_directories", you can simply put the library names in the second argument. Otherwise, put the full path to the library file (i.e., the .a, .so, or .lib file) in the second argument.
    target_link_libraries(transport_solver PRIVATE libmfem.a)

endif()


######################
# MPI-Parallel build #
######################

# If BUILD_FOR_PARALLEL is ON, and an MPI library can be found, create an executable that uses the MPI-parallel build of MFEM with the path to the MPI-parallel build of MFEM provided in cmake/Configuration_Paths.cmake.
if(BUILD_FOR_PARALLEL)
    find_package(MPI) # Find the MPI library, and only build if it is found
    if(MPI_CXX_FOUND)
        # Add the executable for the parallel implementation.
        add_executable(transport_solver_p ${PROJECT_SOURCE_DIR}/src/Transport_Solver_p.cpp ${UTILITY_SOURCE_DIR}/src/JSON_IO.cpp ${UTILITY_SOURCE_DIR}/src/mfem_util.cpp)

        # Tell CMake any directories containing the project's library files (i.e., the .a, .so, or .lib files). This is often the ./lib directories.
        target_link_directories(transport_solver_p PRIVATE ${PAR_MFEM_SOURCE_DIR} ${HYPRE_LIB_DIR} ${METIS_LIB_DIR})
        
        # Tell CMake any directories containing the project's header files (i.e., the .h files). This is often the ./include directories.
        target_include_directories(transport_solver_p PRIVATE ${PAR_MFEM_SOURCE_DIR} ${HYPRE_INCLUDE_DIR} ${METIS_INCLUDE_DIR})
        
        # Define MPI_BUILD, which can be used to trigger MPI_BUILD tasks/functions/classes in the executable
        target_compile_definitions(transport_solver_p PRIVATE MPI_BUILD)

        # Tell CMake any libraries to link to the project's executable (this includes the link to the MPI library). If you told CMake the directories containing the project's libraries in "link_directories", you can simply put the library names in the second argument. Otherwise, put the full path to the library file (i.e., the .a, .so, or .lib file) in the second argument.
        target_link_libraries(transport_solver_p PRIVATE MPI::MPI_CXX libmfem.a libHYPRE.a libmetis.a)
    else()
        message(WARNING "MPI not found: parallel transport solver, transport_solver_p, will not be created.")
    endif()
endif()
