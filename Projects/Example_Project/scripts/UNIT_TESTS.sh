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



#############################################################################
#### Source the bash config file (which provides all the relevant paths) ####
#############################################################################

source ./config.sh
source $BASH_UTILITY_FUNCTIONS_PATH



###############################
#### Define parser options ####
###############################

# Define default values
N_THREADS=-1

# Parse Options
while getopts "j:" opt; do
    case $opt in
        j) N_THREADS=$OPTARG ;;
        /?) echo "Invalid option -$OPTARG" >&2 ;;
    esac
done



##################################################################
#### Perform full development and execution of the MoFA model ####
##################################################################

# Start the timer
SECONDS=0

# Create the MoFA model input file
./make_and_run_input_file_generator.sh -t "uu" -u "1"
# Create the mesh for the closure problems/MoFA model
./make_and_run_mesh_generator.sh -t "u"
# Solve the stokes equations for the flow field for the MoFA model
./make_and_run_stokes_solver.sh -t "u"
# Solve the closure problems
if (( N_THREADS == -1 )); then
    ./full_make_run_scalar_closure_solver.sh
else
    ./transport_closure_solver_ensemble.sh -j "$N_THREADS"
fi
# Solve the MoFA model
./make_and_run_upscaled_model.sh



########################################################################
#### Perform full development and execution of the pore-scale model ####
########################################################################

# Create the pore-scale model input file
./make_and_run_input_file_generator.sh -t "pu" -u "1"
# Create the mesh for the pore-scale model
./make_and_run_mesh_generator.sh -t "p"
# Solve the stokes equations for the flow field for the pore-scale model
./make_and_run_stokes_solver.sh -t "p"
# Solve the pore-scale model
./make_and_run_porescale_solver.sh



###############################
#### Compare model results ####
###############################

# Get the error threshold from the MoFA upscaled configuration file 
EPSILON_KEYPATH=("mesh" "AR" "epsilon")
read -r EPSILON < <(extract_number $MOFA_SIM_CONFIG_PATH "${EPSILON_KEYPATH[@]}")

# Compare model results by computing the absolute error between average solutions
./make_and_run_error_calc.sh -E "$EPSILON"

# Unit 0: Error is ~0.003318
# Unit 1: Error is ~0.068486 or 0.068392



######################
#### Exit message ####
######################

# Stop the timer
echo "$SH_NAME: Tutorial unit test completed in $SECONDS seconds"

# Check status
if [ $? -eq 0 ]; then
  echo "$SH_NAME: All tasks completed successfully!"
else
  echo "$SH_NAME: Execution failed with error code $?"
fi
