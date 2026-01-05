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

//                                      Porescale Averager
//
// Description: This code averages the (transport) porescale solution according to the
//              placement of the averaging regions.
//

#include "mfem.hpp"
#include <vector>
#include <cassert>
#include "JSON_IO.h"
#include "mfem_util.h"



using namespace std;
using namespace mfem;



// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    string FILENAME = "Porescale_Avger.cpp";
    vector<double> L = {1., 1., 0.};
};
static Params globalVars; 



int main(int argc, char *argv[])
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string config_path = "./";

    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";
    
    int order;
    vector<int> isPeriodic = {0, 0, 0};
    int N_steps;
    int output_interval; // Number of time steps before saving the pore-scale and averaged pore-scale solutions
    
    string output_dir = "./";
    string output_file_name_prefix = "c_";
    string output_file_name_suffix = ".gf";
    
    string average_output_dir = "./";
    string average_output_file_name = "avg_sol.txt";
    vector<string> average_solution_keys = {"avg_c"};
    
    string mesh_output_dir = "./";
    string mesh_output_file_name = "mesh.mesh";


    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++) {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc) {
            config_path = argv[i + 1];
            cout << globalVars.FILENAME << ": Configuration path obtained from parser options: " << config_path << endl;
            break;
        }
    }


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
        

        JSONDict porescale_dict = *configData["scalar porescale"];
        
        sub_dict = *porescale_dict["simulation parameters"];
        sub_dict.getValue("order", order);
        sub_dict.getValue("isPeriodic", isPeriodic);
        sub_dict.getValue("N time steps", N_steps);
        sub_dict.getValue("output interval", output_interval);
        
        sub_dict = *porescale_dict["output path"];
        sub_dict.getValue("directory", output_dir);
        sub_dict.getValue("file name prefix", output_file_name_prefix);
        sub_dict.getValue("file name suffix", output_file_name_suffix);

        sub_dict = *porescale_dict["average output path"];
        sub_dict.getValue("directory", average_output_dir);
        sub_dict.getValue("file name", average_output_file_name);
        
        sub_dict = *porescale_dict["mesh output path"];
        sub_dict.getValue("directory", mesh_output_dir);
        sub_dict.getValue("file name", mesh_output_file_name);
    }
    else
    {
        cerr << globalVars.FILENAME << ": main(): Error in loading config file." << endl;
        exit(1);
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_file_path = mesh_output_dir + mesh_output_file_name; // + ".serial";
    string mesh_info_file_path = mesh_info_dir + mesh_info_file_name;
    
    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    args.ParseCheck();

    
    // ===============================================================
    //   Get the mesh, the number of averaging regions, and the defined boundary attributes.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Obtaining mesh and mesh info... ";

    // Get the boundary attributes for the mesh (NOTE: the provided attribute tags are 1 more than the actual Array index used to assign essential boundary conditions (SEE BELOW))
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_file_path);
    vector<int> AR_tags = (*mesh_info["AR"])["tags"];
    vector<double> AR_areas = (vector<double>)(*mesh_info["AR"])["total_areas"]; // Get the AR total areas from the mesh info
    int N_AR = (int)(*mesh_info["AR"])["total_number"];
    string avg_area_type = (*mesh_info["AR"])["area_type"];
    globalVars.L = (*mesh_info["geometry"])["L"];


    // Use the provided mesh file path to define the mesh 
    MeshManager mesh_manager(mesh_file_path);
    
    // Make the mesh periodic if required by isPeriodic
    mesh_manager.MakePeriodic(isPeriodic, globalVars.L);
    
    // Get a pointer to the mesh being used
    Mesh *mesh = mesh_manager.GetMesh();

    cout << "Complete." << endl;


    // ===============================================================
    //   Define a finite element space for the concentration.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Creating finite element spaces... ";

    // Finite element space for concentration
    FiniteElementCollection *fec_c = new H1_FECollection(order, mesh->Dimension());
    FiniteElementSpace *fespace_c = new FiniteElementSpace(mesh, fec_c, 1);

    cout << "Complete." << endl;


    // ===============================================================
    //   Create the averaging operator and obtain associated information.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Creating the averaging operator... ";

    // Initiate the averaging operator
    AveragingOperator avgOp(fespace_c, AR_tags);

    // Get the AR pore areas from the averaging operator
    vector<double> AR_pore_areas(avgOp.GetAR_areas_vector());
    
    // Assert that there is the same number of AR total areas as AR pore areas
    assert (AR_areas.size() == AR_pore_areas.size());

    // Compute the porosities
    vector<double> porosities;
    AveragingOperator::ComputePorosities(AR_pore_areas, AR_areas, porosities);

    cout << "Complete." << endl;


    // ===============================================================
    //   Averaging loop.
    // ===============================================================
    // Initialize the vectors for saving the average solutions
    Vector sol_avg(N_AR); // For interfacing with MFEM entities
    vector<vector<double>> avg_sol; for (int i = 0; i < N_AR; i++) { avg_sol.push_back({}); } // For saving the avg solution in time. Should be avg_sol[i_AR, i_timestep]
    
    // Compute the number of solution snapshots that will need to be averaged
    int N_sol = N_steps / output_interval + 1;
    
    // Begin the averaging loop
    for (int i_sol = 0; i_sol < N_sol; i_sol++)
    {
        // Print the snapshot number
        cout << globalVars.FILENAME << ":   Averaging snapshot " << i_sol << " out of " << N_sol - 1 << "." << endl;

        // Load the previous porescale solution through the GridFunctionCoefficientManager
        string output_file_path = output_dir + output_file_name_prefix + to_string(i_sol) + output_file_name_suffix; // + ".serial";
        GridFunctionManager sol_manager(output_file_path, mesh);
        GridFunction* gf = sol_manager.GetGridFunction();


        // --------------------------------------------------------
        // Creating the averaging operator/fe space on-the-fly was too expensive, despite being more general
        // --------------------------------------------------------

        // Get the finite element space from the manager
        //FiniteElementSpace *fespace = sol_manager.GetGridFunctionFES();
        
        // Initiate the averaging operator
        //AveragingOperator avgOp(fespace, AR_tags);

        // Get the AR pore areas from the averaging operator
        //vector<double> AR_pore_areas(avgOp.GetAR_areas_vector());
        
        // Assert that there is the same number of AR total areas as AR pore areas
        //assert (AR_areas.size() == AR_pore_areas.size());

        // Compute the porosities
        //vector<double> porosities;
        //AveragingOperator::ComputePorosities(AR_pore_areas, AR_areas, porosities);


        // Apply the averaging operator to the solution
        avgOp.ApplyAvgOperator(*gf, sol_avg, porosities);
        
        // Save the average solution in avg_sol
        for (int i = 0; i < sol_avg.Size(); i++) {
            if (std::abs(sol_avg.Elem(i)) <= 1.0E-300) { avg_sol[i].push_back( 0.0 ); }
            else { avg_sol[i].push_back(sol_avg.Elem(i)); }
        }
    }

    // Print the final average solutions (for debugging purposes...)
    cout << "The averaged concentrations :" << endl;
    for (int i = 0; i < sol_avg.Size(); i++) { cout << sol_avg.Elem(i) << endl; }


    // ===============================================================
    //   Save the solutions and mesh.
    // ===============================================================
    // Save the average solutions to the structure, and then the structure to a text file
    JSONDict avg_sol_struct, avg_c_struct, box_data, other_dict;
    for (int i = 0; i < avg_sol.size(); i++) { box_data[to_string(i)] = avg_sol[i]; }
    avg_c_struct[average_solution_keys[0]] = &box_data;
    avg_sol_struct["average solutions"] = &avg_c_struct;
    
    other_dict["N_AR"] = (int)avg_sol.size();
    other_dict["N_steps"] = (int)avg_sol[0].size();
    other_dict["N_sol"] = (int)average_solution_keys.size();
    other_dict["simulation keys"] = average_solution_keys;
    if (avg_area_type == "AR") { other_dict["average_type"] = "superficial"; }
    else { other_dict["average_type"] = "intrinsic"; }
    avg_sol_struct["other"] = &other_dict;
    
    avg_sol_struct.saveToFile(average_output_dir + average_output_file_name);
    

    // ===============================================================
    //   Free the used memory by deleting the pointers.
    // ===============================================================
    delete fespace_c;
    delete fec_c;
    
    return 0;
}
