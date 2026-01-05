#!/bin/bash

#lL.lH.
#
# Copyright (c) 2025, Lawrence Livermore National Security, LLC
# and other MoFA project developers. All Rights reserved.
# See files LICENSE and NOTICE for details. LLNL-CODE-2006961.
#
# This file is part of the MoFA Project. For more information
# and source code availability visit https://github.com/.
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



####################################
#### Record the necessary paths ####
####################################

EXE_DIR="$MOFA_DIR/Solvers/Transport_Solver"
EXE_NAME="transport_solver_p"
RUN_DIR="$PROJECT_DIR/output"

RUN_MODE="serial"
N_CORES="1"

# Parse Options
while getopts "m:n:" opt; do
  case $opt in
    m) RUN_MODE=$OPTARG ;;
    n) N_CORES=$OPTARG ;;
    /?) echo "Invalid option -$OPTARG" >&2 ;;
  esac
done

# Define solver arguments
EXE_ARGS="-C $PORE_SIM_CONFIG_PATH"



##############################################
#### Call the "make_executable.sh" script ####
##############################################

./subscripts/make_executable.sh -d $EXE_DIR



#############################################
#### Call the "run_executable.sh" script ####
#############################################

# Start the timer
SECONDS=0

./subscripts/run_executable.sh -d $EXE_DIR -e $EXE_NAME -r $RUN_DIR -a "$EXE_ARGS" -m "$RUN_MODE" -n "$N_CORES"



######################
#### Exit message ####
######################

# Stop the timer
echo "$SH_NAME: Porescale transport problem solved in $SECONDS seconds"

# Check status
if [ $? -eq 0 ]; then
  echo "$SH_NAME: All tasks completed successfully!"
else
  echo "$SH_NAME: Execution failed with error code $?"
fi
