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

# Parse options
while getopts "d:" opt; do
  case $opt in
    d) EXE_MAIN_DIR=$OPTARG ;;
    /?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

# Check that a directory was provided
if [ "$EXE_MAIN_DIR" == "none" ]; then
  echo "$SH_NAME: Error: Expected directory input (use -d)"
  exit 1
fi



####################################
#### Define the necessary paths ####
####################################

EXE_MAIN_BUILD_DIR="$EXE_MAIN_DIR/build"



#############################
#### Make the executable ####
#############################

# Go to the executable project build directory
cd $EXE_MAIN_BUILD_DIR

# Make the executable
make

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

