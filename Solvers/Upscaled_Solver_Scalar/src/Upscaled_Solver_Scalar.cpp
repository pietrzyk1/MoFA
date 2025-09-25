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

//                                      Transport MoFA Model Solver
//
// Description: This code solves the scalar transport MoFA model.

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include "JSON_IO.h"
#include "mfem_util.h"



using namespace std;
using namespace mfem;



// Declare BC functions. They are defined after "main".
void inlet_BC_func_T(const double &t, double &F);

// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    double epsilon = 0.1;
};
static Params globalVars; 



int main(int argc, char *argv[])
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string config_path = "./";

    string mesh_dir = "./";
    //string mesh_file_name = "mesh.msh";
    unique_ptr<string> mesh_info_dir = make_unique<string>(mesh_dir);
    string mesh_info_file_name = "mesh_info.txt";

    int N_steps = 100;
    double dt = 1.0e-6;
    int output_interval = N_steps; // Number of time steps before saving the pore-scale and averaged pore-scale solutions
    
    string output_dir = "./";
    string output_file_name = "upscaled_c.txt";
    vector<string> average_solution_keys = {"avg_c"};
    
    unique_ptr<string> residual_dir = make_unique<string>(output_dir);
    string residual_file_name = "a_sol.txt";

    int N_iter = 2;
    double omega = 1.0;


    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++)
    {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc)
        {
            config_path = argv[i + 1];
            cout << "Upscaled_Solver_Scalar.cpp: Configuration path obtained from parser options: " << config_path << endl;
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

        sub_dict = *mesh_dict["AR"];
        sub_dict.getValue("epsilon", globalVars.epsilon);
        

        JSONDict upscaled_dict = *configData["upscaled"];
        
        sub_dict = *upscaled_dict["output path"];
        sub_dict.getValue("directory", output_dir);
        sub_dict.getValue("file name", output_file_name);
        
        sub_dict = *upscaled_dict["simulation parameters"];
        sub_dict.getValue("N time steps", N_steps);
        sub_dict.getValue("dt", dt);
        sub_dict.getValue("output interval", output_interval);


        JSONDict closure_dict = *configData["scalar closure"];

        sub_dict = *closure_dict["closure parameters"];
        sub_dict.getValue("max recursion iterations", N_iter);

        sub_dict = *closure_dict["physics parameters"];
        sub_dict.getValue("omega", omega);

        sub_dict = *closure_dict["residual path"];
        sub_dict.getValue("directory", residual_dir);
        sub_dict.getValue("file name", residual_file_name);
    }
    else
    {
        err << "Upscaled_Solver_Scalar.cpp: main: Error in loading config file." << endl;
        return 1;
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_info_file_path = *mesh_info_dir + mesh_info_file_name;
    string residual_file_path = *residual_dir + residual_file_name;
    string output_file_path = output_dir + output_file_name;
    

    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    args.AddOption(&output_file_path, "-O", "--output_file_path", "Output file path for the upscaled solution (use config file for more detailed file names, directories, etc.).");
    args.ParseCheck();
    
    
    // ===============================================================
    // 
    // ===============================================================
    // Load the residual map
    JSONDict closure_residuals_dict; closure_residuals_dict.loadFromFile(residual_file_path);

    // Load the mesh info file (it has N_AR, AR areas, pore areas in ARs, etc.)
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_file_path);
    
    // Load the number of averaging regions and number of iterations
    int N_AR = (int)(*mesh_info["AR"])["total_number"];
    
    
    // Calculate the porosities from the averaging region areas and pore areas in the averaging regions
    vector<double> AR_areas = (vector<double>)(*mesh_info["AR"])["total_areas"];
    vector<double> AR_pore_areas = (vector<double>)(*mesh_info["AR"])["pore_areas"];
    assert(AR_areas.size() == AR_pore_areas.size());
    vector<double> porosities;
    for (int i = 0; i < AR_areas.size(); i++) { porosities.push_back( AR_pore_areas[i]/AR_areas[i] ); }


    // Load the coefficients into the K and M tensors, and Kf and Mf vectors
    vector<double> K_col_major, M_col_major, temp;
    Vector Kf(N_AR), Mf(N_AR);
    
    int i_iter = 0;
    JSONDict res_dict_dummy; res_dict_dummy = *(*(*closure_residuals_dict["avg_c"])["iter_" + to_string(i_iter)])["residuals"];
    for (int i_AR = 0; i_AR < N_AR; i_AR++) // Go through the cases of active AR
    {
        temp = res_dict_dummy[to_string(i_AR)];
        for (int j_AR = 0; j_AR < N_AR; j_AR++) { temp[j_AR] *= 1/porosities[j_AR]; } // divide by the porosity
        K_col_major.insert(K_col_major.end(), temp.begin(), temp.end());
    }
    res_dict_dummy = *(*(*closure_residuals_dict["inlet"])["iter_" + to_string(i_iter)])["residuals"];
    temp = res_dict_dummy["0"];
    for (int i = 0; i < N_AR; i++) { Kf[i] = temp[i]; }


    i_iter = 1;
    res_dict_dummy = *(*(*closure_residuals_dict["avg_c"])["iter_" + to_string(i_iter)])["residuals"];
    for (int i_AR = 0; i_AR < N_AR; i_AR++) // Go through the cases of active AR
    {
        temp = res_dict_dummy[to_string(i_AR)];
        for (int j_AR = 0; j_AR < N_AR; j_AR++) { temp[j_AR] *= 1/porosities[j_AR]; } // divide by the porosity
        M_col_major.insert(M_col_major.end(), temp.begin(), temp.end());
    }
    res_dict_dummy = *(*(*closure_residuals_dict["inlet"])["iter_" + to_string(i_iter)])["residuals"];
    temp = res_dict_dummy["0"];
    for (int i = 0; i < N_AR; i++) { Mf[i] = temp[i]; }



    // Create the dense mass and stiffness matrices
    DenseMatrix M_dense(M_col_major.data(), N_AR, N_AR), K_dense(K_col_major.data(), N_AR, N_AR);    
    K_dense *= 1/globalVars.epsilon/globalVars.epsilon /omega/omega;
    Kf *= 1/globalVars.epsilon/globalVars.epsilon /omega/omega;
    SparseMatrix M(N_AR), K(N_AR);

    for (int i = 0; i < N_AR; i++)
    {
        for (int j = 0; j < N_AR; j++)
        {
            double val = M_dense(i, j);
            if (val != 0.0)
                M.Add(i, j, val);
            val = K_dense(i, j);
            if (val != 0.0)
                K.Add(i, j, val);
        }
    }
    K.Finalize();
    M.Finalize();


    // ===============================================================
    //   Define block structure of the system.
    // ===============================================================
    // Notes: - This defines an array of offsets for each variable. The last component of the Array is the sum of the dimensions
    //          of each block.
    Array<int> block_offsets(2); // The argument should be the number of variables + 1
    block_offsets[0] = 0;
    block_offsets[1] = N_AR;
    block_offsets.PartialSum();

    
    // ===============================================================
    //   Define the solution and right-hand-side (block) vectors and initialize them.
    // ===============================================================
    BlockVector sol_BLK(block_offsets), b_BLK(block_offsets);
    sol_BLK = 0.0;
    b_BLK = 0.0;
    
    // Make references grid functions to the various dependent variables in the block solution. NOTE: Editing these gridfunctions will edit the block solutions
    GridFunction sol_BLK0_ref, b_BLK0_ref;
    sol_BLK0_ref.MakeRef(sol_BLK.GetBlock(0), block_offsets[0]);
    b_BLK0_ref.MakeRef(b_BLK.GetBlock(0), block_offsets[0]);
    
    // Use the reference grid functions to project the initial conditions
    //InitialCondition *IC = new InitialCondition();
    //sol_BLK0_ref.ProjectCoefficient(*IC);


    // ===============================================================
    //   Solve the system in time.
    // ===============================================================
    double t = 0.0, fi, fip1;
    
    // Define the time stepping operator for ODE systems of the form   M * du/dt + K * u = b, and set the initial time
    LinearTimeDependentOperator oper(M, K, b_BLK, dt);
    oper.SetTime(t);
    oper.PrepareImplicitEuler();
    

    // Initialize the vectors for saving the average solutions
    Vector sol_avg(N_AR); // For interfacing with MFEM entities
    vector<vector<double>> avg_sol; for (int i = 0; i < sol_avg.Size(); i++) { avg_sol.push_back({}); } // For saving the avg solution in time. Should be avg_sol[i_AR, i_timestep]
    
    // Initialize block vector for hosting the solution
    BlockVector du_dt(block_offsets), delta_u(block_offsets);
    du_dt = 0.0;
    delta_u = 0.0;


    // Save the initial condition
    for (int i = 0; i < N_AR; i++) { avg_sol[i].push_back(sol_BLK.GetBlock(0).Elem(i)); }


    // Start the time loop
    Vector fip1_coef = Mf;
    fip1_coef *= -1/dt;
    Vector fi_coef = Mf;
    fi_coef *= 1/dt;
    //fi_coef -= Kf; // Explicit
    fip1_coef -= Kf; // Implicit
    for (int i_step = 1; i_step < N_steps + 1; i_step++)
    {
        // Increase the time and print the step
        t += dt;
        cout << "Time step: " << i_step << ", Time: " << t << endl;
        

        // Reset the RHS/b-vector
        b_BLK = 0.0;
        
        // Set the time on the inlet/source boundary condition and project the boundary condition to the
        // appropriate entries of sol_BLK through sol_BLK0_ref, and to those of b_BLK through b_BLK0_ref.
        // This sets the solution to the BC at the correct places.
        inlet_BC_func_T(t, fi);
        inlet_BC_func_T(t + dt, fip1);
        
        // Update the b vector
        Vector fip1_coef2 = fip1_coef;
        fip1_coef2 *= fip1;
        Vector fi_coef2 = fi_coef;
        fi_coef2 *= fi;
        b_BLK.GetBlock(0) += fip1_coef2;
        b_BLK.GetBlock(0) += fi_coef2;
        
        // Apply operator to solve for the time derivative
        //oper.Mult(sol_BLK, du_dt);
        oper.ImplicitSolve(dt, sol_BLK, du_dt);
        
        // Multiply derivative by dt and add to solution to get the solution at the next time step
        delta_u = du_dt;
        delta_u *= dt;
        sol_BLK += delta_u;


        // Save the solution if required
        if (i_step % output_interval == 0)
        {
            // Save the pore-scale solution as a gridfunction, and the average pore-scale solution in the vector structure
            for (int i = 0; i < N_AR; i++) { avg_sol[i].push_back(sol_BLK.GetBlock(0).Elem(i)); }
        }
        
        // Print the maximum value to the consol (for debugging purposes...)
        cout << "Max avg c: " << sol_BLK.GetBlock(0).Max() << ", avg c(0): " << sol_BLK.GetBlock(0).Elem(0) << ", fi: " << fi << endl;
    }

    for (int i = 0; i < sol_BLK.GetBlock(0).Size(); i++)
    {
        cout << sol_BLK.GetBlock(0).Elem(i) << endl;
    }
    

    // ===============================================================
    //   Save the solutions.
    // ===============================================================
    // Save the average solutions to the structure, and then the structure to a text file
    JSONDict avg_sol_struct, avg_c_struct, box_data, other_dict;
    
    for (int i = 0; i < N_AR; i++)
    {
        box_data[to_string(i)] = avg_sol[i];
    }
    avg_c_struct[average_solution_keys[0]] = &box_data;
    avg_sol_struct["average solutions"] = &avg_c_struct;

    other_dict["N_AR"] = N_AR;
    other_dict["N_steps"] = (int)avg_sol[0].size();
    other_dict["N_sol"] = (int)average_solution_keys.size();
    other_dict["simulation keys"] = average_solution_keys;
    avg_sol_struct["other"] = &other_dict;

    avg_sol_struct.saveToFile(output_file_path);
    

    // ===============================================================
    //   Free the used memory by deleting the pointers.
    // ===============================================================
    //delete IC;
    
    return 0;
}










// Define the time-dependent part of the inlet boundary condition c(x,t) = X(x)T(t)
void inlet_BC_func_T(const double &t, double &F)
{
    double a1 = 0.5;
    double b1 = a1 * M_PI * globalVars.epsilon * globalVars.epsilon * 0.1;
    F = a1 - a1 * cos(t * M_PI / b1);
}
