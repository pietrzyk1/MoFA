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



####################################
#### Record the necessary paths ####
####################################

PROJECT_DIR="$SH_DIR/.."
MOFA_DIR="$PROJECT_DIR/../.."


# Define the directory containing the simulation configuration files 
SIM_CONFIG_DIR="$PROJECT_DIR/config"

# Define the paths to the simulation configuration files
MOFA_SIM_CONFIG_FILE_NAME="MoFA_config.txt"
MOFA_SIM_CONFIG_PATH="$SIM_CONFIG_DIR/$MOFA_SIM_CONFIG_FILE_NAME"

PORE_SIM_CONFIG_FILE_NAME="porescale_config.txt"
PORE_SIM_CONFIG_PATH="$SIM_CONFIG_DIR/$PORE_SIM_CONFIG_FILE_NAME"
