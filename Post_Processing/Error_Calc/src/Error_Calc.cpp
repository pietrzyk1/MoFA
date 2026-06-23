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

//                                              Error Calculator
//
// Description:     This code allows the absolute error between averaged pore-scale and MoFA
//                  model results to be computed. Upon finding the absolute error below the
//                  modeling error threshold (provided by the MoFA methodologu), we deem the
//                  model valid.
//

//#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include "JSON_IO.h"
#include "mfem_util.h"



using namespace std;
using namespace mfem;



// Declare functions
double max_in_column(const std::vector<std::vector<double>> &matrix, size_t col);



// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    string FILENAME = "Error_Calc.cpp";
};
static Params globalVars; 



int main(int argc, char *argv[])
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string config_path = "./";

    string mesh_dir = "./";
    string mesh_file_name = "mesh.msh";
    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";

    // For the output error text file
    string output_dir = "./";
    string output_file_name = "abs_error.txt";
    string output_error_gf_file_name_prefix = "c_avg_error_";
    string output_error_gf_file_name_suffix = ".gf";
    string output_error_gf_mesh_file_name = "error_mesh.mesh";

    // For the average pore-scale solutions
    string avg_porescale_dir = "./";
    string avg_porescale_file_name = "avg_c.txt";
    
    // For the upscaled solutions
    string upscaled_dir = "./";
    string upscaled_file_name = "upscaled_c.txt";


    string error_threshold_str;
    double error_threshold = -1.0;

    int N_AR_x;
    int N_AR_y;
    

    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++)
    {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc) {
            config_path = argv[i + 1];
            cout << globalVars.FILENAME << ": Configuration path obtained from parser options: " << config_path << endl;
        }
        if ((string(argv[i]) == "-E" || string(argv[i]) == "--error_threshold") && i + 1 < argc) {
            error_threshold_str = argv[i + 1];
            error_threshold = std::stod(error_threshold_str);
            cout << globalVars.FILENAME << ":Error_Calc.cpp: Error threshold obtained from parser options: " << error_threshold << endl;
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
        
        JSONDict sub_dict = *mesh_dict["output path"];
        sub_dict.getValue("directory", mesh_dir);
        sub_dict.getValue("file name", mesh_file_name);

        sub_dict = *mesh_dict["info path"];
        sub_dict.getValue("directory", mesh_info_dir);
        sub_dict.getValue("file name", mesh_info_file_name);
        
        sub_dict = *mesh_dict["AR"];
        JSONDict subsub_dict = *sub_dict["N_AR"];
        subsub_dict.getValue("x", N_AR_x);
        subsub_dict.getValue("y", N_AR_y);
        

        JSONDict error_dict = *configData["error calc"];

        sub_dict = *error_dict["output path"];
        sub_dict.getValue("directory", output_dir);
        sub_dict.getValue("file name", output_file_name);
        sub_dict.getValue("file name prefix", output_error_gf_file_name_prefix);
        sub_dict.getValue("file name suffix", output_error_gf_file_name_suffix);
        sub_dict.getValue("mesh file name", output_error_gf_mesh_file_name);

        sub_dict = *error_dict["avg porescale path"];
        sub_dict.getValue("directory", avg_porescale_dir);
        sub_dict.getValue("file name", avg_porescale_file_name);

        sub_dict = *error_dict["upscaled path"];
        sub_dict.getValue("directory", upscaled_dir);
        sub_dict.getValue("file name", upscaled_file_name);
    }
    else
    {
        cerr << globalVars.FILENAME << ": main(): Error in loading config file." << endl;
        return 1;
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_file_path = mesh_dir + mesh_file_name;
    string mesh_info_file_path = mesh_info_dir + mesh_info_file_name;
    string output_file_path = output_dir + output_file_name;
    string avg_porescale_file_path = avg_porescale_dir + avg_porescale_file_name;
    string upscaled_file_path = upscaled_dir + upscaled_file_name;

    
    // ===============================================================
    //   Get the mesh, the number of averaging regions, and the defined boundary attributes.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Obtaining mesh and mesh info... ";

    // Get the boundary attributes for the mesh (NOTE: the provided attribute tags are 1 more than the actual Array index used to assign essential boundary conditions (SEE BELOW))
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_file_path);
    vector<int> AR_tags = (*mesh_info["AR"])["tags"];
    string avg_area_type = (*mesh_info["AR"])["area_type"];
    

    int print_error_contour = 1;
            
    
    // Create a mesh, function space, and gridfunction for plotting the average pore-scale solution.
    std::unique_ptr<Mesh> mesh_avg;
    std::unique_ptr<FiniteElementCollection> fec_c_avg;
    std::unique_ptr<FiniteElementSpace> fespace_c_avg;
    std::unique_ptr<GridFunction> sol_error_GF;
    if (print_error_contour == 1)
    {
        // Create a superficial mesh (i.e., mesh of ARs), a function space for that mesh, and a gridfunction on
        // that mesh in case the averaging area type calls for it. If the averaging area type is "AR", it means
        // the superficial average (i.e., average over the AR, not the pore-space) will be used. If it is "pore",
        // then the pore-space average will be used
        if (avg_area_type == "AR") { mesh_avg = std::make_unique<Mesh>(Mesh::MakeCartesian2D(N_AR_x, N_AR_y, Element::Type::QUADRILATERAL, false, N_AR_x, N_AR_y)); }
        else { mesh_avg = std::make_unique<Mesh>(mesh_file_path); } // Use the provided mesh file to define the mesh
        
        // Save the mesh
        mesh_avg->Save(output_dir + output_error_gf_mesh_file_name);

        // Create the finite element space for concentration (this is for plotting the average in the AR or pore space depending on the porosities avaiable)
        fec_c_avg = std::make_unique<DG_FECollection>(1, mesh_avg->Dimension());
        fespace_c_avg = std::make_unique<FiniteElementSpace>(mesh_avg.get(), fec_c_avg.get(), 1);
        
        // Define a grid function for plotting the average concentration within the AR or pore space depending on the porosities avaiable)
        sol_error_GF = std::make_unique<GridFunction>(fespace_c_avg.get());
    }

    cout << "Complete." << endl;


    // ===============================================================
    // 
    // ===============================================================
    // Load the average porescale and upscaled data
    JSONDict avg_porescale_dict, upscaled_dict;
    avg_porescale_dict.loadFromFile(avg_porescale_file_path);
    upscaled_dict.loadFromFile(upscaled_file_path);

    // Get the additional information
    JSONDict avg_porescale_data_info = *avg_porescale_dict["other"];
    JSONDict upscaled_data_info = *upscaled_dict["other"];
    
    // Get N_AR from the two data sources and make sure they are the same
    int N_AR_p = (int)(avg_porescale_data_info["N_AR"]);
    int N_AR_u = (int)(upscaled_data_info["N_AR"]);
    assert(N_AR_p == N_AR_u);

    // Get the number of saved time steps from the two data sources and make sure they are the same
    int N_steps_p = (int)(avg_porescale_data_info["N_steps"]);
    int N_steps_u = (int)(upscaled_data_info["N_steps"]);
    //assert(N_steps_p == N_steps_u);

    // Get N_sol, the number of solution variables, from the two data sources and make sure they are the same 
    int N_sol_p = (int)(avg_porescale_data_info["N_sol"]);
    int N_sol_u = (int)(upscaled_data_info["N_sol"]);
    assert(N_sol_p == N_sol_u);

    // Get sim_keys, the keys for the averages solutions, from the two data sources and make sure they are the same
    vector<string> sim_keys_p = avg_porescale_data_info["simulation keys"];
    vector<string> sim_keys_u = upscaled_data_info["simulation keys"];
    assert(sim_keys_p.size() == sim_keys_u.size());

    // Get the data
    vector<JSONDict*> avg_porescale_data;
    vector<JSONDict*> upscaled_data;

    for (int i = 0; i < N_sol_p; i++)
    {
        avg_porescale_data.push_back( new JSONDict() );
        upscaled_data.push_back( new JSONDict() );

        avg_porescale_data[i] = (*avg_porescale_dict["average solutions"])[sim_keys_p[i]];
        upscaled_data[i] = (*upscaled_dict["average solutions"])[sim_keys_u[i]];
    }
    
    //JSONDict *avg_porescale_data = (*avg_porescale_dict["average_solutions"])["avg_c"];
    //JSONDict *upscaled_data = (*upscaled_dict["average_solutions"])["avg_c"];
    
    

    // ===============================================================
    //   
    // ===============================================================
    double avg_p_c, up_c, abs_error;
    vector<double> avg_p_c_vec, up_c_vec;
    vector<vector<double>> max_error_sol_time;
    vector<vector<vector<double>>> error_data;
    bool validError = true;
    double max_error = -1.0;
    
    for (int i_sol = 0; i_sol < N_sol_p; i_sol++)
    {
        error_data.push_back( vector<vector<double>>() );
        for (int i_AR = 0; i_AR < N_AR_p; i_AR++)
        {
            error_data[i_sol].push_back( vector<double>() );
            avg_p_c_vec = (*avg_porescale_data[i_sol])[to_string(i_AR)];
            up_c_vec = (*upscaled_data[i_sol])[to_string(i_AR)];
            for (int i_step = 0; i_step < N_steps_p; i_step++)
            {
                error_data[i_sol][i_AR].push_back( abs(avg_p_c_vec[i_step] - up_c_vec[i_step]) );
            }
        }
    }


    //Vector AR_label(N_AR_p); for (int i_AR = 0; i_AR < N_AR_p; i_AR++) { AR_label[i_AR] = i_AR; }
    //CreateAvgSolGridFunction(AR_label, *sol_error_GF, avg_area_type, AR_tags);
    //sol_error_GF->Save((output_dir + "TEST_AR_LABEL" + output_error_gf_file_name_suffix).c_str());
    //exit(1);

    // Print the maximum error in space at each time step
    for (int i_sol = 0; i_sol < N_sol_p; i_sol++)
    {
        max_error_sol_time.push_back( vector<double>() );
        cout << endl;
        cout << "Maximum Absolute Error for Solution " << i_sol << " at each Time Step:" << endl;
        for (int i_step = 0; i_step < N_steps_p; i_step++)
        {
            max_error_sol_time[i_sol].push_back( max_in_column(error_data[i_sol], i_step) );
            cout << "Max Abs. Error in space at Time step " << i_step << ": " << max_error_sol_time[i_sol][i_step] << endl;
            
            // Print the average solution on an "AR mesh"
            if (print_error_contour == 1) {
                Vector sol_error(N_AR_p); for (int i_AR = 0; i_AR < N_AR_p; i_AR++) { sol_error[i_AR] = error_data[i_sol][i_AR][i_step]; }
                CreateAvgSolGridFunction(sol_error, *sol_error_GF, avg_area_type, AR_tags);
                sol_error_GF->Save((output_dir + output_error_gf_file_name_prefix + to_string(i_step) + output_error_gf_file_name_suffix).c_str());
            }
        }
    }
    
    // Print the maximum error over all time steps, and all of space, for each solution
    for (int i_sol = 0; i_sol < N_sol_p; i_sol++)
    {
        cout << endl;
        auto maxIterator = std::max_element(max_error_sol_time[i_sol].begin(), max_error_sol_time[i_sol].end());
        cout << "Max Abs. Error over Space and Time for Solution " << i_sol << ": " << *maxIterator << endl;

        // Save the maximum error
        if (max_error < *maxIterator) { max_error = *maxIterator; }
        
        // If error_threshold was provided in the inputs, compare the maximum error over space and time to it
        if (error_threshold != -1.0) {
            if (*maxIterator > error_threshold) {
                validError = false;
            }
        }
    }
    cout << endl;


    // ===============================================================
    //   Save the computed absolute errors.
    // ===============================================================
    // Save the average solutions to the structure, and then the structure to a text file
    JSONDict error_dict, other_dict, avg_c_error;
    vector<JSONDict*> box_data;

    for (int i_sol = 0; i_sol < N_sol_u; i_sol++)
    {
        box_data.push_back( new JSONDict() );
        for (int i_AR = 0; i_AR < N_AR_p; i_AR++)
        {
            (*box_data[i_sol])[to_string(i_AR)] = error_data[i_sol][i_AR];
        }
        avg_c_error[sim_keys_u[i_sol]] = box_data[i_sol];
    }
    error_dict["absolute error"] = &avg_c_error;

    other_dict["N_AR"] = N_AR_u;
    other_dict["N_steps"] = N_steps_p;
    other_dict["N_sol"] = N_sol_u;
    other_dict["simulation keys"] = sim_keys_u;
    error_dict["other"] = &other_dict;

    error_dict.saveToFile(output_file_path);


    // ===============================================================
    //   Free the used memory by deleting the pointers.
    // ===============================================================
    

    // ===============================================================
    //   Report whether the error threshold was satisfied, if it was defined.
    // ===============================================================
    if (error_threshold != -1.0) {
        if (validError) {
            cout << "Error_Calc.cpp: main(): RESULT: Error threshold satisfied! Maximum error was " << max_error << " and error threshold was " << error_threshold << "." << endl;
            cout << endl;
            return 0;
        }
        else {
            cerr << "Error_Calc.cpp: main(): RESULT: Error threshold not satisfied. Maximum error was " << max_error << " and error threshold was " << error_threshold << "." << endl;
            cout << endl;
            return 1;
        }
    }

    return 0;
}





// Function for finding the maximum in a given column of a matrix
double max_in_column(const std::vector<std::vector<double>> &matrix, size_t col)
{
    double max_val = -1;
    for (const auto &row : matrix)
    {
        if (col < row.size())  // safety check
            max_val = std::max(max_val, row[col]);
    }
    return max_val;
}
