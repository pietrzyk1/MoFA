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
#source $BASH_UTILITY_FUNCTIONS_PATH



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


# Check the number of threads provided and determine how many to export for use in the ensemble solve
if (( N_THREADS == -1 )); then
    export OMP_NUM_THREADS=1
else
    export OMP_NUM_THREADS=$N_THREADS
fi



##############################################################
#### Define closure solver file paths and other variables ####
##############################################################

EXE_DIR="$MOFA_DIR/Solvers/Transport_Closure_Solver"
EXE_NAME_ENSEMBLE="transport_solver_ensemble"
RUN_DIR="$PROJECT_DIR/output"

# Define the executable and run directory options
EXE_DIR_OPTION="-E $EXE_DIR"
RUN_DIR_OPTION="-R $RUN_DIR"

# Define the file path to the "run_executable.sh" file
RUN_EXE_PATH="$SH_DIR/subscripts/run_executable.sh"
RUN_EXE_PATH_OPTION="-r $RUN_EXE_PATH"

# Define the configuration file path and corresponding executable argument 
CONFIG_PATH="$MOFA_SIM_CONFIG_PATH"
CONFIG_PATH_OPTION="-C $CONFIG_PATH"

# Compile all options to pass to run_executable.sh
EXE_ARGS="$EXE_DIR_OPTION $RUN_DIR_OPTION $RUN_EXE_PATH_OPTION $CONFIG_PATH_OPTION"



############################################
#### Make the closure solver executable ####
############################################

./subscripts/make_executable.sh -d $EXE_DIR



#############################################
#### Call the "run_executable.sh" script ####
#############################################

# Start the timer
SECONDS=0

# Call the "run_executable.sh" script
./subscripts/run_executable.sh -d $EXE_DIR -e $EXE_NAME_ENSEMBLE -r $RUN_DIR -a "$EXE_ARGS"



######################
#### Exit message ####
######################

# Stop the timer
echo "$SH_NAME: Closure problems solved in $SECONDS seconds"

# Check status
if [ $? -eq 0 ]; then
  echo "$SH_NAME: All tasks completed successfully!"
else
  echo "$SH_NAME: Execution failed with error code $?"
fi
