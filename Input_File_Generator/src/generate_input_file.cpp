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

//                                      Input/Config File Generator
//
// Description: This code and its header can be used to generate project input/config files
//              quickly and reliably. The current options create input/config files for 1.)
//              the quick-start tutorial, and 2.) a general project. The header will be the
//              primary file where changes are made.
//

#include "JSON_IO.h"
#include "generate_input_file_TUTORIAL.h"
#if __has_include( "generate_input_file_STAGED.h" )
#   include "generate_input_file_STAGED.h"
#else
#   include "generate_input_file.h"
#endif



using namespace std;



int main(int argc, char *argv[])
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string project_dir = "./../";
    string sim_type = "upscaled";
    string config_file_name = "MoFA_config.txt";
    

    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++)
    {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--project_dir") && i + 1 < argc)
        {
            project_dir = argv[i + 1];
            cout << "generate_input_file.cpp: Project directory obtained from parser options: " << project_dir << endl;
        }
        else if ((string(argv[i]) == "-t" || string(argv[i]) == "--sim_type") && i + 1 < argc)
        {
            sim_type = argv[i + 1];
            cout << "generate_input_file.cpp: Simulation type obtained from parser options: " << sim_type << endl;
        }
        else if ((string(argv[i]) == "-F" || string(argv[i]) == "--config_file_name") && i + 1 < argc)
        {
            config_file_name = argv[i + 1];
            cout << "generate_input_file.cpp: Config file name obtained from parser options: " << config_file_name << endl;
        }
    }
    
    string project_config_dir = project_dir + "config/";
    

    // ===============================================================
    //   Decide which config file to create based on the input.
    // ===============================================================
    if (sim_type == "upscaled")
    {
        createUpscaledConfig(project_dir, config_file_name);
    }
    else if (sim_type == "upscaled_tutorial")
    {
        createTutorialUpscaledConfig(project_dir, config_file_name);
    }
    else if (sim_type == "porescale")
    {
        createPorescaleConfig(project_dir, config_file_name);
    }
    else if (sim_type == "porescale_tutorial")
    {
        createTutorialPorescaleConfig(project_dir, config_file_name);
    }
    else
    {
        cout << "CRITICAL ERROR: generate_input_file.cpp: main: Provided 'sim_type' (i.e., the argument for '-t') not recognized. Please use '-t u' or '-t p'." << endl;
    }

    return 0;
}
