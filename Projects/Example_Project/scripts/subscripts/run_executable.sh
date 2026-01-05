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



###############################
#### Define parser options ####
###############################

# Define default values
EXE_MAIN_DIR="none"
EXE_NAME="none"
RUN_DIR="none"
ARGS=""
RUN_MODE="serial"
N_CORES="1"

# Parse options
while getopts "d:e:r:a:m:n:" opt; do
  case $opt in
    d) EXE_MAIN_DIR=$OPTARG ;;
    e) EXE_NAME=$OPTARG ;;
    r) RUN_DIR=$OPTARG ;;
    a) ARGS=$OPTARG ;;
    m) RUN_MODE=$OPTARG ;;
    n) N_CORES=$OPTARG ;;
    /?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

# Check that a directory was provided
if [ "$EXE_MAIN_DIR" == "none" ]; then
  echo "$SH_NAME: Error: Expected directory input (use -d)"
  exit 1
fi

# Check that an executable name was provided
if [ "$EXE_NAME" == "none" ]; then
  echo "$SH_NAME: Error: Expected executable name input (use -e)"
  exit 1
fi

# Check that a run directory was provided
if [ "$RUN_DIR" == "none" ]; then
  echo "$SH_NAME: Error: Expected run directory input (use -r)"
  exit 1
fi



####################################
#### Record the necessary paths ####
####################################

EXE_MAIN_BUILD_DIR="$EXE_MAIN_DIR/bin"



############################
#### Run the executable ####
############################

# Go to the run directory for this executable/project
cd $RUN_DIR

# Decide how to call the executable
CHOSEN=0
echo "$SH_NAME: Recognized run mode as '$RUN_MODE'."
if [ "$RUN_MODE" == "srun" ]; then
    # Run the executable in parallel using srun
    srun -n $N_CORES $EXE_MAIN_BUILD_DIR/$EXE_NAME $ARGS
    CHOSEN=1
fi
if [ "$RUN_MODE" == "mpirun" ]; then
    # Run the executable in parallel using mpirun
    mpirun -n $N_CORES $EXE_MAIN_BUILD_DIR/$EXE_NAME $ARGS
    CHOSEN=1
fi
if [ "$RUN_MODE" == "serial" ]; then
    # Run the executable in serial
    $EXE_MAIN_BUILD_DIR/$EXE_NAME $ARGS
    CHOSEN=1
fi
if [ $CHOSEN -eq 0 ]; then
    echo "$SH_NAME: Error: Unrecognized RUN_MODE = $RUN_MODE provided in input (via -m)"
    exit 1
fi

# Go back to this script's directory
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

