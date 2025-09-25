/*lL.lH.
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
lL.lH.*/

#pragma once
#include "JSON_IO.h"
#include "generate_input_file_TUTORIAL.h"
#include "generate_input_file_TOOLS.h"



using namespace std;



// ===============================================================
//   Define default function for generating a config file for the
//   MoFA model formulation process. It is suggested that users
//   copy and paste the contents of the createTutorialUpscaledConfig
//   function here and edit the entries accordingly to create their
//   config files.
// ===============================================================

void createUpscaledConfig(string &project_dir, string &config_file_name)
{
    createTutorialUpscaledConfig(project_dir, config_file_name);
}



// ===============================================================
//   Define default function for generating a config file for the
//   porescale model formulation process. It is suggested that users
//   copy and paste the contents of the createTutorialPorescaleConfig
//   function here and edit the entries accordingly to create their
//   config files.
// ===============================================================

void createPorescaleConfig(string &project_dir, string &config_file_name)
{
    createTutorialPorescaleConfig(project_dir, config_file_name);
}
