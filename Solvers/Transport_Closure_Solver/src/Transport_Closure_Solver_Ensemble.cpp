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

//                                      Transport Closer Solver
//
// Description: This code solves the closure problem associated to scalar transport. This
//              problem is written as the following:
//
//              Advective-Diffusive Transport:
//                                  u . grad(c) - div(grad(c)) = 0
//
//              Average Conditions:
//                                    < c >_{B^n} = 0 or 1
//
//              where u is assumed to be an incompressible, steady flow field and the
//              average condition is 1 in one of the averaging regions, and 0 in others.
//              We assign a Dirichlet condition on the inlet boundary, which can have a
//              value of 0 or 1. All other boundries are no-flux. 2D/3D domains are handled.

#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include "JSON_IO.h"
#include "mfem_util.h"

#include <cstdlib> // for system()
#include <cstdio> // for remove()

// OpenMP
#ifdef USE_OPENMP
    #include <omp.h>
#else
    cerr << "Transport_Closure_Solver_Ensemble.cpp: CRITICAL ERROR: CMake was unable to find OpenMP during build. As such, this code cannot run.";
    exit(1);
#endif




using namespace std;
using namespace mfem;




// Define a global "parameters struct"
struct Params
{
    string FILENAME = "Transport_Closure_Solver_Ensemble.cpp";
};
static Params globalVars; 


int main(int argc, char *argv[])
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options 
    // ===============================================================
    string config_path = "./";
    string config_path_ensemble = "./";
    string run_exe_file_path = "./";
    
    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";

    int max_recursive_iters = 2;
    
    string residual_output_dir = "./";
    string residual_output_file_name_prefix = "a_sol";
    string residual_output_file_name_suffix = ".txt";
    
    int useInlet = 0;
    
    string log_file_file_name_prefix = "log_thread";
    string log_file_file_name_suffix = ".txt";





    // TODO: Get from the bash code
    string EXE_DIR = "./";
    string EXE_NAME = "transport_solver";
    string RUN_DIR = "./";

    // Solver flags
    string CONFIG_PATH_FLAG = "-C";
    string AR_FORCING_FLAG = "-a";
    string BC_FORCING_FLAG = "-i";
    string CLOSURE_PATH_FLAG = "-N";
    string FORCING_FUNCTION_FILE_NAME_FLAG = "-F";
    string RECURSIVE_ITER_FLAG = "-I";
    string SAVE_MESH_FLAG = "-s";
    string RESIDUAL_OUTPUT_FILE_NAME_FLAG = "-r";

    string SOLVE_MODE_FLAG = "-S";
    string SOLVE_MODE = "ensemble";

    // File names
    string CLOSURE_FILE_NAME_PREFIX = "chi";
    string CLOSURE_FILE_NAME_SUFFIX = ".gf";
    string CLOSURE_FILE_NAME_PREFIX_AR = "_AR";
    string CLOSURE_FILE_NAME_PREFIX_BC = "_BC";
    string CLOSURE_FILE_ITER_ID = "I";
    

    


    // TODO: Get the config file made in bash for defining the previous variables, and load them



    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++)
    {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc)
        {
            config_path = argv[i + 1];
            cout << globalVars.FILENAME << ": Configuration path obtained from parser options: " << config_path << endl;
            break;
        }
    }

    // Search for the ensemble configuration path
    //for (int i = 1; i < argc; i++)
    //{
    //    if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc)
    //    {
    //        ensemble_config_path = argv[i + 1];
    //        cout << "Transport_Closure_Solver_Ensemble.cpp: Configuration path obtained from parser options: " << config_path << endl;
    //        break;
    //    }
    //}


    // ===============================================================
    //   Define config file struct and load file. Use data to initialize variables
    // ===============================================================
    JSONDict configData;
    int config_output = configData.loadFromFile(config_path);

    // If loading the config file was successful, adjust the values of the predefined variables
    if (config_output == 0)
    {
        JSONDict mesh_dict = *configData["mesh"];
        
        JSONDict sub_dict = *mesh_dict["info path"];
        sub_dict.getValue("directory", mesh_info_dir);
        sub_dict.getValue("file name", mesh_info_file_name);


        JSONDict closure_dict = *configData["scalar closure"];
        
        sub_dict = *closure_dict["simulation parameters"];
        sub_dict.getValue("use inlet", useInlet);

        sub_dict = *closure_dict["closure parameters"];
        sub_dict.getValue("max recursion iterations", max_recursive_iters);

        sub_dict = *closure_dict["residual path"];
        sub_dict.getValue("directory", residual_output_dir);
        sub_dict.getValue("file name prefix", residual_output_file_name_prefix);
        sub_dict.getValue("file name suffix", residual_output_file_name_suffix);

        sub_dict = *closure_dict["closure path"];
        sub_dict.getValue("file name prefix", CLOSURE_FILE_NAME_PREFIX);
        sub_dict.getValue("file name suffix", CLOSURE_FILE_NAME_SUFFIX);
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_info_file_path = mesh_info_dir + mesh_info_file_name;
    string CONFIG_FILE_PATH_OPTION = CONFIG_PATH_FLAG + " " + config_path;
    string SOLVE_MODE_OPTION = SOLVE_MODE_FLAG + " " + SOLVE_MODE;
        



    // TODO; Probably do not need to use the command line much, unless we always call solver from here (all-in-one)
    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    args.AddOption(&mesh_info_file_path, "-m", "--mesh_info_file_path", "Mesh info file path (best to define in the config file).");
    args.AddOption(&EXE_DIR, "-E", "--EXE_DIR", "Directory of the executable (i.e., transport solver).");
    args.AddOption(&RUN_DIR, "-R", "--RUN_DIR", "Directory of the simulation run output (i.e., the 'output' directory of the project).");
    args.AddOption(&run_exe_file_path, "-r", "--run_exe_file_path", "File path to the 'run_executable.sh' shell script.");
    //args.AddOption(&solve_mode, "-S", "--solve_mode", "Define how the code is being run. Options: 'serial' or 'ensemble' (i.e., closure problems are being solved as an ensamble in parallel)");
    args.ParseCheck();


    // ===============================================================
    //   Load variables from mesh info file
    // ===============================================================
    JSONDict mesh_info;
    mesh_info.loadFromFile(mesh_info_file_path);
    int N_AR = (int)(*mesh_info["AR"])["total_number"]; // Get the number of averaging regions defined in the mesh file
    

    // ===============================================================
    //   Solve the AR closure problems
    // ===============================================================
    // Declare a vector of strings to save the closure residual files' names
    vector<string> residual_files;
    bool completed_any_AR_problems = false;
    
    if (true)
    {
        #pragma omp parallel
        {
            // Get the current thread's ID
            int thread_ID = omp_get_thread_num();

            // Print "Hello!" from each thread to the consol
            #pragma omp critical
            {
                cout << globalVars.FILENAME << ": Thread " << thread_ID << ": Hello!" << std::endl;
            }

            // Have each thread create a log file for their solver outputs
            string log_file_file_path = residual_output_dir + log_file_file_name_prefix + to_string(thread_ID) + log_file_file_name_suffix;
            ofstream LOG(log_file_file_path);

            // Declare variables
            string msg;
            string BC_FORCING_OPTION;
            string FORCING_OPTION;
            string RECURSIVE_ITER_OPTION;
            string CLOSURE_FILE_NAME_ITER_AND_END;
            string CLOSURE_FILE_NAME;
            string CLOSURE_FILE_NAME_OPTION;
            string SAVE_MESH_OPTION;
            string RUN_ARGS;
            string string_command;
            string PREVIOUS_CLOSURE_FILE_NAME;

            // Define closure residual output file name and option. Add the name to the residual_files variable
            string RESIDUAL_OUTPUT_FILE_NAME = residual_output_file_name_prefix + to_string(thread_ID) + residual_output_file_name_suffix;
            string RESIDUAL_OUTPUT_FILE_OPTION = RESIDUAL_OUTPUT_FILE_NAME_FLAG + " " + RESIDUAL_OUTPUT_FILE_NAME;
            #pragma omp critical
            {
                residual_files.push_back( RESIDUAL_OUTPUT_FILE_NAME );
            }
            
            // Clear the closure residual file before starting to make sure only the thread's work is saved in there 
            remove((residual_output_dir + RESIDUAL_OUTPUT_FILE_NAME).c_str());
            
            
            #pragma omp for schedule(dynamic)
            for (int i_AR = 0; i_AR < N_AR; i_AR++)
            {
                for (int i_ITER = 0; i_ITER < max_recursive_iters; i_ITER++)
                {
                    // Print the current closure problem to the consol
                    msg = "Thread " + to_string(thread_ID) + ": Starting closure problem for AR " + to_string(i_AR) + ", iteration " + to_string(i_ITER) + "." ;
                    #pragma omp critical
                    {
                        cout << globalVars.FILENAME << ": " << msg << endl;
                    }
                    LOG << msg << "\n" << endl;
                    LOG.close();

                    // Ensure that the BC forcing is off. We do the BC forcing separately; after the AR forcing
                    FORCING_OPTION = BC_FORCING_FLAG + " " + "0";

                    // Determine whether the system will be forced by the AR forcing (if i_ITER=0) or the previous closure variable solution
                    FORCING_OPTION = FORCING_OPTION + " " + AR_FORCING_FLAG + " " + to_string(i_AR);
                    if (i_ITER > 0) { FORCING_OPTION = FORCING_OPTION + " " + FORCING_FUNCTION_FILE_NAME_FLAG + " " + PREVIOUS_CLOSURE_FILE_NAME; }
                    
                    // Create the recursive iteration option
                    RECURSIVE_ITER_OPTION = RECURSIVE_ITER_FLAG + " " +  to_string(i_ITER);
                    
                    // Piece together the closure solution file name and option
                    //CLOSURE_FILE_NAME_ITER_AND_END = CLOSURE_FILE_ITER_ID + to_string(i_ITER) + CLOSURE_FILE_NAME_SUFFIX;
                    //CLOSURE_FILE_NAME = CLOSURE_FILE_NAME_PREFIX + CLOSURE_FILE_NAME_PREFIX_AR + to_string(i_AR) + CLOSURE_FILE_NAME_ITER_AND_END;
                    CLOSURE_FILE_NAME = CLOSURE_FILE_NAME_PREFIX_AR + to_string(i_AR) + CLOSURE_FILE_ITER_ID + to_string(i_ITER);
                    CLOSURE_FILE_NAME_OPTION = CLOSURE_PATH_FLAG + " " + CLOSURE_FILE_NAME;
                    
                    // Determine if the simulation run will save the mesh (only needs to be done by one thread, one time)
                    if ( i_AR == 0 && i_ITER == 0 ) { SAVE_MESH_OPTION = SAVE_MESH_FLAG + " " + to_string(1); }
                    else { SAVE_MESH_OPTION = SAVE_MESH_FLAG + " " + to_string(0); }


                    // Put together the options to give the solver
                    RUN_ARGS = CONFIG_FILE_PATH_OPTION + " " + FORCING_OPTION + " " //+ BC_FORCING_OPTION + " "
                            + CLOSURE_FILE_NAME_OPTION + " " + RECURSIVE_ITER_OPTION + " " + SOLVE_MODE_OPTION
                            + " " + SAVE_MESH_OPTION + " " + RESIDUAL_OUTPUT_FILE_OPTION;
                    
                    // Call the "run_executable.sh" script
                    string_command = run_exe_file_path + " -d " + EXE_DIR + " -e " + EXE_NAME + " -r " + RUN_DIR + " -a \"" + RUN_ARGS + "\" >> " + log_file_file_path + " 2>&1";
                    int executable_result = system(string_command.c_str());
                    LOG.open(log_file_file_path, std::ios::app);
                    LOG << "\n";
                    

                    // Collect the closure file name before moving to the next iteration
                    PREVIOUS_CLOSURE_FILE_NAME = CLOSURE_FILE_NAME_PREFIX + CLOSURE_FILE_NAME + CLOSURE_FILE_NAME_SUFFIX;

                    // Print the current closure problem to the consol
                    msg = "Thread " + to_string(thread_ID) + ": Finished closure problem for AR " + to_string(i_AR) + ", iteration " + to_string(i_ITER) + ".";
                    #pragma omp critical
                    {
                        cout << globalVars.FILENAME << ": " << msg << endl;
                        completed_any_AR_problems = true;
                    }
                    LOG << msg << "\n" << endl;
                }
            }

            // Close the log file
            LOG.close();
        }
    }


    // ===============================================================
    //   Solve the BC closure problems
    // ===============================================================

    if (useInlet == 1)
    {
        //#pragma omp parallel
        int i_BC = 1;

        // Get the current thread's ID
        int thread_ID = 0; //omp_get_thread_num();

        // Have each thread create a log file for their solver outputs
        string log_file_file_path = residual_output_dir + log_file_file_name_prefix + to_string(thread_ID) + log_file_file_name_suffix;
        std::ios_base::openmode stream_mode;
        if (completed_any_AR_problems) { stream_mode = std::ios::app; }
        else { remove(log_file_file_path.c_str()); }
        ofstream LOG(log_file_file_path, stream_mode);
        
        // Declare variables
        string msg;
        string FORCING_OPTION;
        string RECURSIVE_ITER_OPTION;
        string CLOSURE_FILE_NAME_ITER_AND_END;
        string CLOSURE_FILE_NAME;
        string CLOSURE_FILE_NAME_OPTION;
        string RUN_ARGS;
        string string_command;
        string PREVIOUS_CLOSURE_FILE_NAME;

        // Define closure residual output file name and option. Add the name to the residual_files variable
        string RESIDUAL_OUTPUT_FILE_NAME = residual_output_file_name_prefix + to_string(thread_ID) + residual_output_file_name_suffix;
        string RESIDUAL_OUTPUT_FILE_OPTION = RESIDUAL_OUTPUT_FILE_NAME_FLAG + " " + RESIDUAL_OUTPUT_FILE_NAME;
        if (!completed_any_AR_problems) { residual_files.push_back( RESIDUAL_OUTPUT_FILE_NAME ); }
        

        for (int i_ITER = 0; i_ITER < max_recursive_iters; i_ITER++)
        {
            // Print the current closure problem to the consol
            msg = "Thread " + to_string(thread_ID) + ": Starting closure problem BC " + to_string(i_BC) + ", iteration " + to_string(i_ITER) + "." ;
            cout << globalVars.FILENAME << ": " << msg << endl;
            LOG << msg << "\n" << endl;
            LOG.close();

            // Determine whether the system will be forced by the AR forcing (if i_ITER=0) or the previous closure variable solution
            FORCING_OPTION = BC_FORCING_FLAG + " " + "1";
            if (i_ITER > 0 ) { FORCING_OPTION = FORCING_OPTION + " " + FORCING_FUNCTION_FILE_NAME_FLAG + " " + PREVIOUS_CLOSURE_FILE_NAME; }
            
            // Create the recursive iteration option
            RECURSIVE_ITER_OPTION = RECURSIVE_ITER_FLAG + " " +  to_string(i_ITER);
            
            // Piece together the closure solution file name and option
            //CLOSURE_FILE_NAME_ITER_AND_END = CLOSURE_FILE_ITER_ID + to_string(i_ITER) + CLOSURE_FILE_NAME_SUFFIX;
            //CLOSURE_FILE_NAME = CLOSURE_FILE_NAME_PREFIX + CLOSURE_FILE_NAME_PREFIX_BC + to_string(i_BC) + CLOSURE_FILE_NAME_ITER_AND_END;
            CLOSURE_FILE_NAME = CLOSURE_FILE_NAME_PREFIX_BC + to_string(i_BC) + CLOSURE_FILE_ITER_ID + to_string(i_ITER);
            CLOSURE_FILE_NAME_OPTION = CLOSURE_PATH_FLAG + " " + CLOSURE_FILE_NAME;

            // Piece together the thread's residual closure file name, and the associated option
            RESIDUAL_OUTPUT_FILE_NAME = residual_output_file_name_prefix + to_string(thread_ID) + residual_output_file_name_suffix;
            RESIDUAL_OUTPUT_FILE_OPTION = RESIDUAL_OUTPUT_FILE_NAME_FLAG + " " + RESIDUAL_OUTPUT_FILE_NAME;
            
            
            // Put together the options to give the solver
            RUN_ARGS = CONFIG_FILE_PATH_OPTION + " " + FORCING_OPTION + " "
                     + CLOSURE_FILE_NAME_OPTION + " " + RECURSIVE_ITER_OPTION
                     + " " + SOLVE_MODE_OPTION + " " + RESIDUAL_OUTPUT_FILE_OPTION;
            
            // Call the "run_executable.sh" script
            string_command = run_exe_file_path + " -d " + EXE_DIR + " -e " + EXE_NAME + " -r " + RUN_DIR + " -a \"" + RUN_ARGS + "\" >> " + log_file_file_path + " 2>&1";
            int executable_result = system(string_command.c_str());
            LOG.open(log_file_file_path, std::ios::app);
            LOG << "\n";
            
            
            // Collect the closure file name before moving to the next iteration
            PREVIOUS_CLOSURE_FILE_NAME = CLOSURE_FILE_NAME_PREFIX + CLOSURE_FILE_NAME + CLOSURE_FILE_NAME_SUFFIX;

            // Print the current closure problem to the consol
            msg = "Thread " + to_string(thread_ID) + ": Finished closure problem for BC " + to_string(i_BC) + ", iteration " + to_string(i_ITER) + ".";
            cout << globalVars.FILENAME << ": " << msg << endl;
            LOG << msg << "\n" << endl;
        }

        // Close the log file
        LOG.close();
    }

    cout << globalVars.FILENAME << ": All closure problems solved successfully!" << endl;        


    // ===============================================================
    //   Combine closure problem residual files
    // ===============================================================
    cout << globalVars.FILENAME << ": Compiling closure residuals to a single output file...";        

    // Load the first residual file to a compiled JSONDict
    JSONDict compiled_residual_file;
    int successful_load = compiled_residual_file.loadFromFile(residual_output_dir + residual_files[0]);
    assert (successful_load == 0);
    
    for (int i = 1; i < residual_files.size(); i++)
    {
        // Load the next residual file produced by a thread
        JSONDict thread_output_file;
        successful_load = thread_output_file.loadFromFile(residual_output_dir + residual_files[i]);
        assert (successful_load == 0);

        // Merge the information in the residual file to the compiled residual file
        compiled_residual_file.merge(thread_output_file);
    }

    // Save the combined residual file
    compiled_residual_file.saveToFile(residual_output_dir + residual_output_file_name_prefix + residual_output_file_name_suffix);

    cout << " Complete!" << endl;
    cout << globalVars.FILENAME << ": Deleting individual closure residual files (the ones from each thread, since they have been combined)...";

    // Delete the residual files from the threads (since we now have the compiled file)
    for (int i = 0; i < residual_files.size(); i++) { remove((residual_output_dir + residual_files[i]).c_str()); }

    cout << " Complete!" << endl;


    // ===============================================================
    //   Done
    // ===============================================================
    cout << globalVars.FILENAME << ": Ensemble solve complete!" << endl;
    
    return 0;
}
