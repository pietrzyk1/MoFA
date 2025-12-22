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
// Description: This code solves the scalar transport MoFA model. This model is of the form
//
//                          M * d<c>/dt + K * <c> + Mf * df/dt + Kf * f = 0
//
//              Here, <c> is a vector of all average concentrations, M and K are matrices, Mf
//              and Kf are vectors, and f is a scalar. In general, the closure residuals obtained
//              from the closure problems fill M, K, Mf, and Kf. 
//

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include "JSON_IO.h"
#include "mfem_util.h"
#if __has_include( "mfem_util_STAGED.h" )
#   include "mfem_util_STAGED.h"
constexpr bool genResFormAvailable = true;
#else
constexpr bool genResFormAvailable = false;
class GeneralizedResidualManager
{
public:
    GeneralizedResidualManager();
    GeneralizedResidualManager(double alpha_, double beta_, double gamma_);
    void SetParams(double alpha_, double beta_, double gamma_);
    void ConstructParamVecs(int N_AR_global, const vector<int> &AR_inds_loc2glob);
    void ImplementAlpha(BilinearForm *&varf);
    void ImplementBeta(SparseMatrix &cavg);
    void ImplementGamma(SparseMatrix &BM_avgavg, const vector<int> &AR_inds_loc2glob);
    void AlterGammaVec(const vector<int> &AR_tags, const vector<vector<int>> &AR_neighbors, int active_AR_global, int procedureID = -1);
    vector<double>* GetParamVec(string param_name);
    bool LoadParamDicts(JSONDict &eq_dict);
    void ImplementGeneralizedResidual(vector<double> &temp, vector<string> keys, bool isKMat = false);
    void CreateTransform(JSONDict &closure_residuals_dict, vector<double> &porosities, int N_AR, double epsilon, double omega);
    void ModifyTransform(vector<double> &temp, vector<string> keys, bool addOne = false);
    void ApplyTransform(vector<vector<double>> &avg_sol, const Vector &temp_sols, int N_AR, double fip1, double fi, double dt);
    bool GetUseGenResForm();
    DenseMatrix* GetTransform();
    Vector GetTransf();
    Vector GetTransdfdt();
};
#endif



using namespace std;
using namespace mfem;



// Declare BC functions. They are defined after "main".
void inlet_BC_func_T(const double &t, double &F);

// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    string FILENAME = "Upscaled_Solver_Scalar.cpp";
    double epsilon = 0.1;
    double BC_frequency_scale = 1.0;
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
    unique_ptr<string> mesh_info_dir = make_unique<string>(mesh_dir);
    string mesh_info_file_name = "mesh_info.txt";

    int N_steps = 100;
    double dt = 1.0e-6;
    int output_interval = N_steps; // Number of time steps before saving the pore-scale and averaged pore-scale solutions
    
    string output_dir = "./";
    string output_file_name_prefix = "upscaled_sol_";
    string output_file_name_suffix = ".gf";
    
    unique_ptr<string> average_output_dir = make_unique<string>(output_dir);
    string average_output_file_name = "upscaled_sol.txt";
    vector<string> average_solution_keys = {"avg_c"};
    
    unique_ptr<string> mesh_output_dir = make_unique<string>(output_dir);
    string mesh_output_file_name = "mesh.mesh";
    

    unique_ptr<string> residual_dir = make_unique<string>(output_dir);
    string residual_file_name_prefix = "a_sol";
    string residual_file_name_suffix = ".txt";

    int N_iter = 2;
    double omega = 1.0;

    int N_AR_x;
    int N_AR_y;


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
        sub_dict.getValue("epsilon", globalVars.epsilon);
        JSONDict subsub_dict = *sub_dict["N_AR"];
        subsub_dict.getValue("x", N_AR_x);
        subsub_dict.getValue("y", N_AR_y);
        

        JSONDict upscaled_dict = *configData["upscaled"];
        
        sub_dict = *upscaled_dict["output path"];
        sub_dict.getValue("directory", output_dir);
        sub_dict.getValue("file name prefix", output_file_name_prefix);
        sub_dict.getValue("file name suffix", output_file_name_suffix);
        
        sub_dict = *upscaled_dict["average output path"];
        sub_dict.getValue("directory", average_output_dir);
        sub_dict.getValue("file name", average_output_file_name);
        
        
        sub_dict = *upscaled_dict["mesh output path"];
        sub_dict.getValue("directory", mesh_output_dir);
        sub_dict.getValue("file name", mesh_output_file_name);
        
        sub_dict = *upscaled_dict["simulation parameters"];
        sub_dict.getValue("N time steps", N_steps);
        sub_dict.getValue("dt", dt);
        sub_dict.getValue("output interval", output_interval);
        sub_dict.getValue("BC frequency scale", globalVars.BC_frequency_scale);
        

        JSONDict closure_dict = *configData["scalar closure"];

        sub_dict = *closure_dict["closure parameters"];
        sub_dict.getValue("max recursion iterations", N_iter);
        
        sub_dict = *closure_dict["physics parameters"];
        sub_dict.getValue("omega", omega);
        
        sub_dict = *closure_dict["residual path"];
        sub_dict.getValue("directory", residual_dir);
        sub_dict.getValue("file name prefix", residual_file_name_prefix);
        sub_dict.getValue("file name suffix", residual_file_name_suffix);
    }
    else
    {
        cerr << globalVars.FILENAME << ": main: Error in loading config file." << endl;
        return 1;
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_file_path = mesh_dir + mesh_file_name;
    string mesh_info_file_path = *mesh_info_dir + mesh_info_file_name;
    string residual_file_path = *residual_dir + residual_file_name_prefix + residual_file_name_suffix;
    
    string average_output_file_path = *average_output_dir + average_output_file_name;
    string mesh_output_file_path = *mesh_output_dir + mesh_output_file_name;

    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    //args.AddOption(&output_file_path, "-O", "--output_file_path", "Output file path for the upscaled solution (use config file for more detailed file names, directories, etc.).");
    args.ParseCheck();
    

    // ===============================================================
    //   Get the mesh and define a finite element space for plotting the avg solution
    // ===============================================================
    // Load the mesh info file, and the relevant quantities from it
    JSONDict mesh_info;
    mesh_info.loadFromFile(mesh_info_file_path);
    int N_AR = (int)(*mesh_info["AR"])["total_number"];
    vector<int> AR_tags = (*mesh_info["AR"])["tags"];
    string avg_area_type = (*mesh_info["AR"])["area_type"];

    // Define pointers for the mesh, finite element collection and space, and avg concentration gridfunction
    std::unique_ptr<Mesh> mesh;
    
    // Create a superficial mesh (i.e., mesh of ARs), a function space for that mesh, and a gridfunction on
    // that mesh in case the averaging area type calls for it. If the averaging area type is "AR", it means
    // the superficial average (i.e., average over the AR, not the pore-space) will be used. If it is "pore",
    // then the pore-space average will be used
    if (avg_area_type == "AR") { mesh = std::make_unique<Mesh>(Mesh::MakeCartesian2D(N_AR_x, N_AR_y, Element::Type::QUADRILATERAL, false, N_AR_x, N_AR_y)); }
    // Use the provided mesh file to define the mesh
    else { mesh = std::make_unique<Mesh>(mesh_file_path); }
    
    // Create the finite element space for concentration (this is for plotting the average in the AR or pore space depending on the porosities avaiable)
    FiniteElementCollection *fec_c = new DG_FECollection(1, mesh->Dimension());
    FiniteElementSpace *fespace_c = new FiniteElementSpace(mesh.get(), fec_c, 1);
    
    // Define a grid function for plotting the average concentration within the AR or pore space depending on the porosities avaiable)
    GridFunction *avg_sol_GF = new GridFunction(fespace_c);
    
    
    // ===============================================================
    // 
    // ===============================================================
    // Calculate the porosities from the averaging region areas and pore areas in the averaging regions
    vector<double> AR_areas = (vector<double>)(*mesh_info["AR"])["total_areas"];
    vector<double> AR_pore_areas = (vector<double>)(*mesh_info["AR"])["pore_areas"];;
    vector<double> porosities;
    
    // See if the AR pore areas are available in the mesh info file
    //JSONDict mesh_info_dict_AR = *mesh_info["AR"];
    //mesh_info_dict_AR.getValue("pore_areas", AR_pore_areas);
    
    // If not, use the. AveragingOperator class to get the AR pore areas
    //if (AR_pore_areas.size() == 0)
    //{
    //    cout << "UTS:   Creating averaging operator to compute pore areas... ";
    //    
    //    AveragingOperator *avgOp = new AveragingOperator(fespace_c, AR_tags);
    //    vector<double> dummy_vec(avgOp->GetAR_areas_vector());
    //    AR_pore_areas = dummy_vec;
    //    delete avgOp;
    //    
    //    cout << "Complete." << endl;
    //}
    //else { cout << "UTS:   Obtained AR pore areas from mesh info file. " << endl; }
    
    // Compute the porosities
    cout << "UTS:   Computing porosities... ";

    assert(AR_areas.size() == AR_pore_areas.size());
    for (int i = 0; i < AR_areas.size(); i++) { porosities.push_back( AR_pore_areas[i]/AR_areas[i] ); }
    
    cout << "Complete." << endl;


    // ===============================================================
    // 
    // ===============================================================
    // The assumed MoFA model is: M * d<c>/dt + K * <c> + Mf * df/dt + Kf * f = 0. Here, <c> is a vector of all
    // average concentrations, M and K are matrices, Mf and Kf are vectors, and f is a scalar
    //
    // Initialize variables for loading closure residuals into K and M tensors, and Kf and Mf vectors
    vector<double> K_col_major, M_col_major, temp;
    Vector Kf(N_AR), Mf(N_AR);    

    
    // ===============================================================
    // Load the residual map
    JSONDict closure_residuals_dict; closure_residuals_dict.loadFromFile(residual_file_path);
    
    
    // ===============================================================
    //   Load residuals into the K matrix
    
    string sim_type = "avg_c";
    int i_iter = 0; string i_iter_str = to_string(i_iter);
    string eq_key = "transport eq";
    string comp_key = "component 0";
    
    // Load the dictionary holding the residuals (and generalized closure residual parameters)
    JSONDict eq_dict = *(*(*closure_residuals_dict[sim_type])["iter_" + i_iter_str])[eq_key];
    
    // Obtain the residual dictionary for the specified simulation type, time derivative iteration, and equation
    JSONDict res_dict_dummy; eq_dict.getValue("residuals", res_dict_dummy);

    // If the code for implementing the general residual form is available, and any of the corresponding parameters
    // are different than what would collapse to the canonical residual form, implement the general residual form
    GeneralizedResidualManager gen_res_manager;
    bool useGenResForm = false;
    if (genResFormAvailable) { useGenResForm = gen_res_manager.LoadParamDicts(eq_dict); }
    
    // Go through the cases of active AR
    for (int i_AR = 0; i_AR < N_AR; i_AR++) {
        // Get the raw residual vector
        string i_AR_str = to_string(i_AR);
        temp = (*(res_dict_dummy[comp_key]))[i_AR_str];
        // If generalized-residual parameters were found in eq_dict by gen_res_manager, appy them to temp
        if (useGenResForm) { gen_res_manager.ImplementGeneralizedResidual(temp, {comp_key, i_AR_str}, true); }
        // Divide by the porosity (to potentially make the averages "superficial averages")
        for (int j_AR = 0; j_AR < N_AR; j_AR++) { temp[j_AR] *= 1/porosities[i_AR]; }
        // Add the temp vector to the K matrix
        K_col_major.insert(K_col_major.end(), temp.begin(), temp.end());
    }


    // ===============================================================
    //   Load residuals into the K vector

    sim_type = "inlet";
    i_iter = 0; i_iter_str = to_string(i_iter);
    eq_key = "transport eq";
    comp_key = "component 0";
    

    // Load the dictionary holding the residuals (and generalized closure residual parameters)
    eq_dict = *(*(*closure_residuals_dict[sim_type])["iter_" + i_iter_str])[eq_key];
    
    // Obtain the residual dictionary for the specified simulation type, time derivative iteration, and equation
    res_dict_dummy.clear(); eq_dict.getValue("residuals", res_dict_dummy);

    // If the code for implementing the general residual form is available, and any of the corresponding parameters
    // are different than what would collapse to the canonical residual form, implement the general residual form
    useGenResForm = false;
    if (genResFormAvailable) { useGenResForm = gen_res_manager.LoadParamDicts(eq_dict); }
    
    // Get the raw residual vector
    temp = (*res_dict_dummy[comp_key])["0"];
    
    // If generalized-residual parameters were found in eq_dict by gen_res_manager, appy them to temp
    if (useGenResForm) { gen_res_manager.ImplementGeneralizedResidual(temp, {comp_key, "0"}, false); }
    
    // Assign the temp vector to the Kf vector
    for (int i = 0; i < N_AR; i++) { Kf[i] = temp[i]; }
    

    // ===============================================================
    //   Load residuals into the M matrix

    sim_type = "avg_c";
    i_iter = 1; i_iter_str = to_string(i_iter);
    eq_key = "transport eq";
    comp_key = "component 0";
    
    
    // Load the dictionary holding the residuals (and generalized closure residual parameters)
    eq_dict = *(*(*closure_residuals_dict[sim_type])["iter_" + i_iter_str])[eq_key];
    
    // Obtain the residual dictionary for the specified simulation type, time derivative iteration, and equation
    res_dict_dummy; eq_dict.getValue("residuals", res_dict_dummy);

    // If the code for implementing the general residual form is available, and any of the corresponding parameters
    // are different than what would collapse to the canonical residual form, implement the general residual form
    useGenResForm = false;
    if (genResFormAvailable) { useGenResForm = gen_res_manager.LoadParamDicts(eq_dict); }
    
    // Go through the cases of active AR
    for (int i_AR = 0; i_AR < N_AR; i_AR++) {
        // Get the raw residual vector
        string i_AR_str = to_string(i_AR);
        temp = (*(res_dict_dummy[comp_key]))[i_AR_str];
        // If generalized-residual parameters were found in eq_dict by gen_res_manager, appy them to temp
        if (useGenResForm) { gen_res_manager.ImplementGeneralizedResidual(temp, {comp_key, i_AR_str}, false); }
        // Divide by the porosity (to potentially make the averages "superficial averages")
        for (int j_AR = 0; j_AR < N_AR; j_AR++) { temp[j_AR] *= 1/porosities[i_AR]; }
        // Add the temp vector to the M matrix
        M_col_major.insert(M_col_major.end(), temp.begin(), temp.end());
    }


    // ===============================================================
    //   Load residuals into the M vector

    sim_type = "inlet";
    i_iter = 1; i_iter_str = to_string(i_iter);
    eq_key = "transport eq";
    comp_key = "component 0";
    

    // Load the dictionary holding the residuals (and generalized closure residual parameters)
    eq_dict = *(*(*closure_residuals_dict[sim_type])["iter_" + i_iter_str])[eq_key];
    
    // Obtain the residual dictionary for the specified simulation type, time derivative iteration, and equation
    res_dict_dummy.clear(); eq_dict.getValue("residuals", res_dict_dummy);

    // If the code for implementing the general residual form is available, and any of the corresponding parameters
    // are different than what would collapse to the canonical residual form, implement the general residual form
    useGenResForm = false;
    if (genResFormAvailable) { useGenResForm = gen_res_manager.LoadParamDicts(eq_dict); }
    
    // Get the raw residual vector
    temp = (*res_dict_dummy[comp_key])["0"];
    
    // If generalized-residual parameters were found in eq_dict by gen_res_manager, appy them to temp
    if (useGenResForm) { gen_res_manager.ImplementGeneralizedResidual(temp, {comp_key, "0"}, false); }
    
    // Assign the temp vector to the Mf vector
    for (int i = 0; i < N_AR; i++) { Mf[i] = temp[i]; }


    // ===============================================================
    //   Create the actual K and M matrices

    // Create the dense mass and stiffness matrices
    DenseMatrix M_dense(M_col_major.data(), N_AR, N_AR), K_dense(K_col_major.data(), N_AR, N_AR);    
    K_dense *= 1/globalVars.epsilon/globalVars.epsilon /omega/omega;
    Kf *= 1/globalVars.epsilon/globalVars.epsilon /omega/omega;
    SparseMatrix M(N_AR), K(N_AR);

    for (int i = 0; i < N_AR; i++) {
        for (int j = 0; j < N_AR; j++) {
            double val = M_dense(i, j);
            if (val != 0.0) { M.Add(i, j, val); }
            val = K_dense(i, j);
            if (val != 0.0) { K.Add(i, j, val); }
        }
    }
    K.Finalize();
    M.Finalize();

    
    // ===============================================================
    //   Create the transform if the code for the general residual
    //   form is available
    // ===============================================================
    if (genResFormAvailable) { gen_res_manager.CreateTransform(closure_residuals_dict, porosities, N_AR, globalVars.epsilon, omega); }


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
    CreateAvgSolGridFunction(sol_BLK.GetBlock(0), *avg_sol_GF, avg_area_type);
    avg_sol_GF->Save((output_dir + output_file_name_prefix + to_string(avg_sol[0].size() - 1) + output_file_name_suffix).c_str());
    

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
            // Save the average solution in the vector structure. If the code for the generalized residual form is available, use it
            Vector avg_sol_saved(N_AR);
            if (genResFormAvailable) {
                Vector temp_sols(2*N_AR);
                for (int i_temp = 0; i_temp < N_AR; i_temp++) { temp_sols.Elem(i_temp) = sol_BLK.GetBlock(0).Elem(i_temp); }
                for (int i_temp = 0; i_temp < N_AR; i_temp++) { temp_sols.Elem(i_temp + N_AR) = du_dt.GetBlock(0).Elem(i_temp); }
                gen_res_manager.ApplyTransform(avg_sol, temp_sols, N_AR, fip1, fi, dt);
            }
            else { for (int i = 0; i < N_AR; i++) { avg_sol[i].push_back(sol_BLK.GetBlock(0).Elem(i)); } }
            
            // Save the average solution as a gridfunction
            for (int i = 0; i < N_AR; i++) { avg_sol_saved[i] = avg_sol[i][avg_sol[i].size() - 1]; }
            //CreateAvgSolGridFunction(sol_BLK.GetBlock(0), *avg_sol_GF, avg_area_type);
            CreateAvgSolGridFunction(avg_sol_saved, *avg_sol_GF, avg_area_type);
            avg_sol_GF->Save((output_dir + output_file_name_prefix + to_string(avg_sol[0].size() - 1) + output_file_name_suffix).c_str());
        }
        
        // Print the maximum value to the consol (for debugging purposes...)
        cout << "Max avg c: " << sol_BLK.GetBlock(0).Max() << ", avg c(0): " << sol_BLK.GetBlock(0).Elem(0) << ", fi: " << fi << endl;
    }

    for (int i = 0; i < sol_BLK.GetBlock(0).Size(); i++)
    {
        cout << sol_BLK.GetBlock(0).Elem(i) << endl;
    }
    

    // ===============================================================
    //   Save the mesh and solutions.
    // ===============================================================
    // Save the mesh
    mesh->Save(mesh_output_file_path);

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
    if (avg_area_type == "AR") { other_dict["average_type"] = "superficial"; }
    else { other_dict["average_type"] = "intrinsic"; }
    avg_sol_struct["other"] = &other_dict;

    avg_sol_struct.saveToFile(average_output_file_path);
    

    // ===============================================================
    //   Free the used memory by deleting the pointers.
    // ===============================================================
    //delete IC;
    delete avg_sol_GF;
    delete fespace_c;
    delete fec_c;
    //delete mesh;
    
    return 0;
}










// Define the time-dependent part of the inlet boundary condition c(x,t) = X(x)T(t)
void inlet_BC_func_T(const double &t, double &F)
{
    double a1 = 0.5;
    double b1 = a1 * M_PI * globalVars.epsilon * globalVars.epsilon * globalVars.BC_frequency_scale;// * 0.75;
    F = a1 - a1 * cos(t * M_PI / b1);
}
