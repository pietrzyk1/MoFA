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



##########################################
#### Define closure solver file paths ####
##########################################

EXE_DIR="$MOFA_DIR/Solvers/Transport_Closure_Solver"
EXE_NAME="transport_solver"
RUN_DIR="$PROJECT_DIR/output"



############################################
#### Make the closure solver executable ####
############################################

./subscripts/make_executable.sh -d $EXE_DIR



#######################################################################
#### Define flags, iteration parameters, and simulation parameters ####
#######################################################################

# Solver config file
CONFIG_FILE_PATH="$MOFA_SIM_CONFIG_PATH"

# Flags
CONFIG_PATH_FLAG="-C"
AR_FORCING_FLAG="-a"
BC_FORCING_FLAG="-i"
CLOSURE_PATH_FLAG="-N"
FORCING_FUNCTION_FILE_NAME_FLAG="-F"
RECURSIVE_ITER_FLAG="-I"

# File names
CLOSURE_FILE_NAME_PREFIX="chi"
CLOSURE_FILE_NAME_PREFIX_AR="_AR"
CLOSURE_FILE_NAME_PREFIX_BC="_BC"
CLOSURE_FILE_ITER_ID="I"
CLOSURE_FILE_NAME_SUFFIX=".gf"

# Simulation parameters
N_AR=10 # TODO Need to automate getting these...
N_ITER=2

# Initialize
CONFIG_FILE_PATH_OPTION="$CONFIG_PATH_FLAG $CONFIG_FILE_PATH"
PREVIOUS_CLOSURE_FILE_NAME=""

for ((i_AR = 0; i_AR < $N_AR; i_AR++))
do
    for ((i_ITER = 0; i_ITER < $N_ITER; i_ITER++))
    do
        # Ensure that the BC forcing is off. We do the BC forcing separately; after the AR forcing
        BC_FORCING_OPTION="$BC_FORCING_FLAG 0"
        
        # Determine whether the system will be forced by the AR forcing (if i_ITER=0) or the previous closure variable solution
        FORCING_OPTION="$AR_FORCING_FLAG $i_AR"
        if [ $i_ITER -gt 0 ]; then
            FORCING_OPTION="$FORCING_OPTION $FORCING_FUNCTION_FILE_NAME_FLAG $PREVIOUS_CLOSURE_FILE_NAME"
        fi

        # Create the recursive iteration option
        RECURSIVE_ITER_OPTION="$RECURSIVE_ITER_FLAG $i_ITER"
        
        # Piece together the closure solution file name and option
        CLOSURE_FILE_NAME_ITER_AND_END="$CLOSURE_FILE_ITER_ID$i_ITER$CLOSURE_FILE_NAME_SUFFIX"
        CLOSURE_FILE_NAME="$CLOSURE_FILE_NAME_PREFIX$CLOSURE_FILE_NAME_PREFIX_AR$i_AR$CLOSURE_FILE_NAME_ITER_AND_END"
        CLOSURE_FILE_NAME_OPTION="$CLOSURE_PATH_FLAG $CLOSURE_FILE_NAME"
        

        # Put together the options to give the solver
        RUN_ARGS="$CONFIG_FILE_PATH_OPTION $FORCING_OPTION $BC_FORCING_OPTION $CLOSURE_FILE_NAME_OPTION $RECURSIVE_ITER_OPTION"
        
        # Call the "run_executable.sh" script
        ./subscripts/run_executable.sh -d $EXE_DIR -e $EXE_NAME -r $RUN_DIR -a "$RUN_ARGS"
        
        # Collect the closure file name before moving to the next iteration
        PREVIOUS_CLOSURE_FILE_NAME="$CLOSURE_FILE_NAME"
    done
done



#######################################################################
#### Define flags, iteration parameters, and simulation parameters ####
#######################################################################

i_BC=1

for ((i_ITER = 0; i_ITER < $N_ITER; i_ITER++))
do
    # Determine whether the system will be forced by the AR forcing (if i_ITER=0) or the previous closure variable solution
    FORCING_OPTION="$BC_FORCING_FLAG 1"
    if [ $i_ITER -gt 0 ]; then
        FORCING_OPTION="$FORCING_OPTION $FORCING_FUNCTION_FILE_NAME_FLAG $PREVIOUS_CLOSURE_FILE_NAME"
    fi
    
    # Create the recursive iteration option
    RECURSIVE_ITER_OPTION="$RECURSIVE_ITER_FLAG $i_ITER"
    
    # Piece together the closure solution file name and option
    CLOSURE_FILE_NAME_ITER_AND_END="$CLOSURE_FILE_ITER_ID$i_ITER$CLOSURE_FILE_NAME_SUFFIX"
    CLOSURE_FILE_NAME="$CLOSURE_FILE_NAME_PREFIX$CLOSURE_FILE_NAME_PREFIX_BC$i_BC$CLOSURE_FILE_NAME_ITER_AND_END"
    CLOSURE_FILE_NAME_OPTION="$CLOSURE_PATH_FLAG $CLOSURE_FILE_NAME"
    

    # Put together the options to give the solver
    RUN_ARGS="$CONFIG_FILE_PATH_OPTION $FORCING_OPTION $CLOSURE_FILE_NAME_OPTION $RECURSIVE_ITER_OPTION"
    
    # Call the "run_executable.sh" script
    ./subscripts/run_executable.sh -d $EXE_DIR -e $EXE_NAME -r $RUN_DIR -a "$RUN_ARGS"

    # Collect the closure file name before moving to the next iteration
    PREVIOUS_CLOSURE_FILE_NAME="$CLOSURE_FILE_NAME"
done



######################
#### Exit message ####
######################

# Check status
if [ $? -eq 0 ]; then
  echo "$SH_NAME: All tasks completed successfully!"
else
  echo "$SH_NAME: Execution failed with error code $?"
fi
