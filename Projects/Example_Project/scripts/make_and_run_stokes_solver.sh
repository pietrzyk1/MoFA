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



####################################
#### Record the necessary paths ####
####################################

EXE_DIR="$MOFA_DIR/Solvers/Stokes_Solver"
EXE_NAME="stokes_solver"
RUN_DIR="$PROJECT_DIR/output"
SIM_TYPE="none"

# Parse Options
while getopts "t:" opt; do
  case $opt in
    t) SIM_TYPE=$OPTARG ;;
    /?) echo "Invalid option -$OPTARG" >&2 ;;
  esac
done

# Check the simulation type provided to determine the configuration file
if [ "$SIM_TYPE" == "none" ]; then
  echo "$SH_NAME: Error: Expected input -t u or -t p"
  exit 1
elif [ "$SIM_TYPE" == "u" ]; then
  CONFIG_PATH="$MOFA_SIM_CONFIG_PATH"
elif [ "$SIM_TYPE" == "p" ]; then
  CONFIG_PATH="$PORE_SIM_CONFIG_PATH"
else
  echo "$SH_NAME: Error: Unrecognized input for -t (expected u or p)"
  exit 1
fi

EXE_ARGS="-C $CONFIG_PATH"



##############################################
#### Call the "make_executable.sh" script ####
##############################################

./subscripts/make_executable.sh -d $EXE_DIR



#############################################
#### Call the "run_executable.sh" script ####
#############################################

./subscripts/run_executable.sh -d $EXE_DIR -e $EXE_NAME -r $RUN_DIR -a "$EXE_ARGS"



######################
#### Exit message ####
######################

# Check status
if [ $? -eq 0 ]; then
  echo "$SH_NAME: All tasks completed successfully!"
else
  echo "$SH_NAME: Execution failed with error code $?"
fi
