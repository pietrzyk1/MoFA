#!/bin/bash

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

#############################
#### Shell script set up ####
#############################

# Tell shell script to exit immediately if a command fails
set -e

# Get this shell script's name
SH_NAME=$(basename "$0")

# Get the directory that this shell script is in
SH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"



#####################################################
#### Define the paths to the "build" directories ####
#####################################################

BUILD_DIRS=("Input_File_Generator/build/" "Mesh_Tools/build/" "Solvers/Stokes_Solver/build/" "Solvers/Transport_Closure_Solver/build/" "Solvers/Transport_Solver/build/" "Solvers/Upscaled_Solver_Scalar/build/" "Post_Processing/Error_Calc/build/")

N_BUILD_DIRS=${#BUILD_DIRS[@]}



############################################################################
#### Change directory to provided directories and execute cmake command ####
############################################################################

# Navigate to the top directory of MoFA
TOP_DIR="$SH_DIR/../"
cd $TOP_DIR

# Navigate to the provided directories, execute the cmake command, and return to the top directory of MoFA
for ((i = 0; i < $N_BUILD_DIRS; i++))
do
    cd ${BUILD_DIRS[i]}
    rm -R *
    cmake ../ -DCMAKE_EXPORT_COMPILE_COMMANDS=1 #-DCMAKE_CXX_COMPILER=mpicxx # Uncomment for HPC
    cd $TOP_DIR
done

# Navigate back to the shell script directory
cd $SH_DIR



######################
#### Exit message ####
######################

# Check status
if [ $? -eq 0 ]; then
    echo "$SH_NAME: All tasks completed successfully!"
else
    echo "$SH_NAME: Execution failed with error code $?"
fi
