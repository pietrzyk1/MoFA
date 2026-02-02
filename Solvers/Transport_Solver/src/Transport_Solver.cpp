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

//                                      Transport Solver
//
// Description: This code solves the dimensionless scalar transport equation. This
//              equation is written as the following:
//
//              Advective-Diffusive Transport:
//                           dc/dt + u . grad(c) - div(grad(c)) = 0
//
//              where u is assumed to be an incompressible, steady flow field. We
//              assign a Dirichlet condition on the inlet boundary, which oscillates
//              in time. All other boundries are no-flux. 2D/3D domains are handled.
//

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
void inlet_BC_func_X(const Vector &p, double &F);

// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    string FILENAME = "Transport_Solver.cpp";
    vector<double> L = {1., 1., 0.};
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
    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";
    
    string fluid_velocity_dir = "./";
    string fluid_velocity_file_name = "u_sol.gf";
    string fluid_velocity_mesh_file_name = "mesh.mesh";

    int order = 2;
    vector<int> isPeriodic = {0, 0, 0};
    int active_advection = 1;
    int N_steps = 5000;
    double dt = 0.00001;
    int output_interval = 1000; // Number of time steps before saving the pore-scale and averaged pore-scale solutions
    int importFluidMesh = 0;

    double Pe = 10.0;
    double omega = 1.0;

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
        

        JSONDict stokes_dict = *configData["stokes"];
        
        sub_dict = *stokes_dict["output path"];
        sub_dict.getValue("directory", fluid_velocity_dir);
        sub_dict.getValue("velocity file name", fluid_velocity_file_name);
        sub_dict.getValue("mesh file name", fluid_velocity_mesh_file_name);


        JSONDict porescale_dict = *configData["scalar porescale"];
        
        sub_dict = *porescale_dict["simulation parameters"];
        sub_dict.getValue("order", order);
        sub_dict.getValue("active advection", active_advection);
        sub_dict.getValue("isPeriodic", isPeriodic);
        sub_dict.getValue("N time steps", N_steps);
        sub_dict.getValue("dt", dt);
        sub_dict.getValue("output interval", output_interval);
        sub_dict.getValue("BC frequency scale", globalVars.BC_frequency_scale);
        sub_dict.getValue("import fluid mesh", importFluidMesh);
        
        sub_dict = *porescale_dict["physical parameters"];
        sub_dict.getValue("Pe", Pe);
        sub_dict.getValue("omega", omega);
        
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
        cerr << globalVars.FILENAME << ": main: Error in loading config file." << endl;
        return 1;
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_file_path = mesh_dir + mesh_file_name;
    string mesh_info_file_path = mesh_info_dir + mesh_info_file_name;
    string fluid_velocity_file_path = fluid_velocity_dir + fluid_velocity_file_name;
    string fluid_velocity_mesh_file_path = fluid_velocity_dir + fluid_velocity_mesh_file_name;
    
    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    args.AddOption(&mesh_file_path, "-M", "--mesh_file_path", "Mesh file path (best to define in the config file).");
    args.AddOption(&mesh_info_file_path, "-m", "--mesh_info_file_path", "Mesh attributes file path (best to define in the config file).");
    //args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree) of scalar closure variable.");
    args.AddOption(&output_dir, "-O", "--output_dir", "Output directory for the pore-scale solution (use config file for more detailed file names, directories, etc.).");
    args.AddOption(&fluid_velocity_file_path, "-v", "--fluid_velocity_file_path", "Fluid velocity solution file path (best to define in the config file).");
    args.AddOption(&fluid_velocity_mesh_file_path, "-V", "--fluid_velocity_mesh_file_path", "Fluid velocity mesh file path (best to define in the config file).");
    args.AddOption(&active_advection, "-A", "--active_advection", "Toggle whether to consider advection or not.");
    args.AddOption(&Pe, "-P", "--Peclet_Number", "Peclet Number to use.");
    args.ParseCheck();

    
    // ===============================================================
    //   Define variables based on the options provided in the parser
    // ===============================================================
    // When using the parallel solver, if advection is active/the fluid velocity was previous solved and used here, the mesh saved with the fluid velocity must be used. Also, use consistent BCs (i.e., both fluid and transport simulations should be non-periodic)
    if (active_advection == 1 && importFluidMesh == 1) { mesh_file_path = fluid_velocity_mesh_file_path; }

    
    // ===============================================================
    //   Get the mesh, the number of averaging regions, and the defined boundary attributes.
    // ===============================================================
    // Get the boundary attributes for the mesh (NOTE: the provided attribute tags are 1 more than the actual Array index used to assign essential boundary conditions (SEE BELOW))
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_file_path);
    vector<int> AR_tags = (*mesh_info["AR"])["tags"];
    string avg_area_type = (*mesh_info["AR"])["area_type"];
    globalVars.L = (*mesh_info["geometry"])["L"];

    // Use the provided mesh file path to define the mesh 
    MeshManager mesh_manager(mesh_file_path);
    
    // Make the mesh periodic if required by isPeriodic
    mesh_manager.MakePeriodic(isPeriodic, globalVars.L);
    
    // Get a pointer to the mesh being used
    Mesh *mesh = mesh_manager.GetMesh();
    

    // ===============================================================
    //   Define finite element spaces for the concentration.
    // ===============================================================
    // Finite element space for concentration
    FiniteElementCollection *fec_c = new H1_FECollection(order, mesh->Dimension());
    FiniteElementSpace *fespace_c = new FiniteElementSpace(mesh, fec_c, 1);
    
    // Print the number of unknowns for concentration and total
    cout << globalVars.FILENAME << ": Number of unknowns in concentration: " << fespace_c->GetTrueVSize() << endl;
    cout << globalVars.FILENAME << ": Number of unknowns in total: " << fespace_c->GetTrueVSize() << endl;

    
    // ===============================================================
    //   Load the fluid velocity for advection
    // ===============================================================
    // Declare pointers for loading/preparing the fluid velocity grid function coefficient
    std::unique_ptr<VectorGridFunctionCoefficientManager> fluid_velocity_vgfc_manager;
    VectorGridFunctionCoefficient *fluid_velocity = nullptr;
    
    // If advection is considered...
    if (active_advection == 1)
    {
        // Create the fluid velocity vector grid function coefficient manager
        cout << globalVars.FILENAME << ": Loading fluid velocity... ";
        fluid_velocity_vgfc_manager = std::make_unique<VectorGridFunctionCoefficientManager>(fluid_velocity_file_path, fluid_velocity_mesh_file_path, Pe);
        cout << "Complete." << endl;
    
        // Create the fluid velocity VectorGridFunctionCoefficient
        cout << globalVars.FILENAME << ": Creating fluid velocity VectorGridFunctionCoefficient... ";
        fluid_velocity_vgfc_manager->MakeVectorGridFunctionCoefficient();
        fluid_velocity = fluid_velocity_vgfc_manager->GetVectorGridFunctionCoefficient();
        cout << "Complete." << endl;
    }


    // ===============================================================
    //   Define/Prepare boundary conditions.
    // ===============================================================
    // Here, we are defining two marker arrays to mark 1.) the inlet, and 2.) the outlet.
    Array<int> marker_inlet_BC(mesh->bdr_attributes.Max()); // Define a marker array for the inlet BC
    marker_inlet_BC = 0; // Initialize marker array to "don't apply any essential BC"
    marker_inlet_BC[(int)(*mesh_info["scalar_closure"])["inlet1"] - 1] = 1; // Set the first bdr attr group (i.e., physical group defined as 1 in gmsh file) as being the inlet
    DirichlettBC *inlet_BC = new DirichlettBC(*fespace_c, marker_inlet_BC);
    inlet_BC->SetTemporallyDependentFunction( inlet_BC_func_T );
    inlet_BC->SetSpatiallyDependentFunction( inlet_BC_func_X );
    

    // ===============================================================
    //   Define block structure of the system.
    // ===============================================================
    // Notes: - This defines an array of offsets for each variable. The last component of the Array is the sum of the dimensions
    //          of each block.
    Array<int> block_offsets(2); // The argument should be the number of variables + 1
    block_offsets[0] = 0;
    block_offsets[1] = fespace_c->GetVSize();
    block_offsets.PartialSum();
    

    // ===============================================================
    //   Define the solution and right-hand-side (block) vectors and initialize them.
    // ===============================================================
    BlockVector sol_BLK(block_offsets), b_BLK(block_offsets);
    sol_BLK = 0.0;
    b_BLK = 0.0;
    
    // Make references grid functions to the various dependent variables in the block solution. NOTE: Editing these gridfunctions will edit the block solutions
    GridFunction sol_BLK0_ref, b_BLK0_ref;
    sol_BLK0_ref.MakeRef(fespace_c, sol_BLK.GetBlock(0), block_offsets[0]);
    b_BLK0_ref.MakeRef(fespace_c, b_BLK.GetBlock(0), block_offsets[0]);
    
    // Use the reference grid functions to project the initial conditions
    //InitialCondition *IC = new InitialCondition();
    //sol_BLK0_ref.ProjectCoefficient(*IC);


    // ===============================================================
    //   Create the bilinear form for the diffusion term.
    // ===============================================================
    // Initiate the bilinear form of the diffusion term
    BilinearForm *varf_cdiff(new BilinearForm(fespace_c));
    varf_cdiff->AddDomainIntegrator(new DiffusionIntegrator);
    if (active_advection == 1) { varf_cdiff->AddDomainIntegrator(new ConvectionIntegrator(*fluid_velocity)); }
    varf_cdiff->Assemble();
    
    // Before eliminating the rows and cols for the essential BCs, use the variational form to get a matrix for the vdofs that the BC affects
    inlet_BC->SetOperator(*varf_cdiff);

    // Apply the inlet BC to the diffusion bilinear form
    varf_cdiff->EliminateEssentialBC(marker_inlet_BC); // Sets the rows and cols in "varf_viscous" corresponding to DOFs marked by "marker_inlet_BC" to 0.0 and assigns a diagonal as 1.0. Then, for the same rows, sets the components of b_BLK to those of c_ref (if given this input).
    
    // Finish the bilinear form by "finalizing"
    varf_cdiff->Finalize();


    // ===============================================================
    //   Create the bilinear form for the mass matrix.
    // ===============================================================
    // Initiate the bilinear form of the mass matrix
    BilinearForm *varf_mass(new BilinearForm(fespace_c));
    ConstantCoefficient dimlessTimeScale(omega);
    varf_mass->AddDomainIntegrator(new MassIntegrator(dimlessTimeScale));
    varf_mass->Assemble();
    // varf_cdiff->EliminateEssentialBC(marker_inlet_BC, Operator::DIAG_ONE); // This might need to be here
    varf_mass->Finalize();
    
    
    // ===============================================================
    //   Create the bilinear form for the averaging operator/forcing unknowns.
    // ===============================================================
    // Use the AveragingOperator class to obtain a matrix that can be multiplied by the solution vector to obtain the average solution in each AR
    AveragingOperator *avgOp = new AveragingOperator(fespace_c, AR_tags);

    // Get the AR pore areas from the averaging operator
    vector<double> AR_pore_areas(avgOp->GetAR_areas_vector());
    
    // Get the AR total areas from the mesh info and assert that there is the same number of AR pore areas
    vector<double> AR_areas = (vector<double>)(*mesh_info["AR"])["total_areas"];
    assert(AR_areas.size() == AR_pore_areas.size());

    // Compute the porosities
    vector<double> porosities;
    AveragingOperator::ComputePorosities(AR_pore_areas, AR_areas, porosities);    
    
    
    // ===============================================================
    //   Solve the system in time.
    // ===============================================================
    double t = 0.0;
    
    // Initiate pointers to the sparse matrices of the diffusion bilinear form and mass matrix bilinear form
    SparseMatrix &BM_cdiff(varf_cdiff->SpMat()), &BM_mass(varf_mass->SpMat());

    // Define the time stepping operator for ODE systems of the form   M * du/dt + K * u = b, and set the initial time
    LinearTimeDependentOperator oper(BM_mass, BM_cdiff, b_BLK, dt);
    oper.PrepareImplicitEuler();
    //oper.PrepareExplicitEuler();
    oper.SetTime(t);
    
    // Initialize the vectors for saving the average solutions
    Vector sol_avg(avgOp->GetNAR()); // For interfacing with MFEM entities
    vector<vector<double>> avg_sol; for (int i = 0; i < sol_avg.Size(); i++) { avg_sol.push_back({}); } // For saving the avg solution in time. Should be avg_sol[i_AR, i_timestep]
    
    // Initialize block vector for hosting the solution
    BlockVector du_dt(block_offsets), delta_u(block_offsets);
    du_dt = 0.0;
    delta_u = 0.0;


    // Save the initial condition for the pore-scale and average pore-scale results
    int save_count = 0;
    sol_BLK0_ref.Save((output_dir + output_file_name_prefix + to_string(save_count) + output_file_name_suffix).c_str());
    save_count += 1;
    avgOp->ApplyAvgOperator(sol_BLK, sol_avg, porosities);
    for (int i = 0; i < sol_avg.Size(); i++) { avg_sol[i].push_back(sol_avg.Elem(i)); }


    StopWatch timer;
    // Start the time loop
    for (int i_step = 1; i_step < N_steps + 1; i_step++)
    {
        // Increase the time and print the step
        t += dt;
        cout << "Time step: " << i_step << ", Time: " << t << endl;
        
        
        // Reset the RHS/b-vector
        b_BLK = 0.0;
        
        // // Set the time on the inlet/source boundary condition and project the boundary condition to the
        // appropriate entries of sol_BLK through sol_BLK0_ref, and to those of b_BLK through b_BLK0_ref.
        // This sets the solution to the BC at the correct places.
        inlet_BC->SetTime(t);
        inlet_BC->ProjectToGridFunction(sol_BLK0_ref);
        inlet_BC->ProjectToGridFunction(b_BLK0_ref);
        // Update the peripheral vdofs in the b vector (i.e., the vdofs that are affected by the Dirichlett BC). Note, this ADDS the contribution to b_BLK; it does not rewrite the entries
        inlet_BC->UpdateLinearFormVector(sol_BLK, b_BLK);
        
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
            sol_BLK0_ref.Save((output_dir + output_file_name_prefix + to_string(save_count) + output_file_name_suffix).c_str());
            save_count += 1;
            avgOp->ApplyAvgOperator(sol_BLK, sol_avg, porosities);
            for (int i = 0; i < sol_avg.Size(); i++) {
                if (std::abs(sol_avg.Elem(i)) <= 1.0E-300) { avg_sol[i].push_back( 0.0 ); }
                else { avg_sol[i].push_back(sol_avg.Elem(i)); }
            }
        }
        

        // Print the maximum value to the consol (for debugging purposes...)
        cout << "Max c: " << sol_BLK.GetBlock(0).Max() << endl;
    }

    // Print the final average solutions (for debugging purposes...)
    avgOp->ApplyAvgOperator(sol_BLK, sol_avg, porosities);
    cout << "The averaged concentrations :" << endl;
    for (int i = 0; i < sol_avg.Size(); i++)
    {
        cout << sol_avg.Elem(i) << endl;
    }


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

    
    // Save the AR pore areas back to the mesh info file
    JSONDict AR_dict_temp = *mesh_info["AR"];
    AR_dict_temp["pore_areas"] = AR_pore_areas;
    mesh_info["AR"] = &AR_dict_temp;
    mesh_info.saveToFile(mesh_info_file_path);
    
    // Save the mesh
    mesh->Save(mesh_output_dir + mesh_output_file_name);


    // ===============================================================
    //   Free the used memory by deleting the pointers.
    // ===============================================================
    delete varf_cdiff;
    //delete IC;
    delete inlet_BC;
    delete fespace_c;
    delete fec_c;
    
    return 0;
}










// Define the time-dependent part of the inlet boundary condition c(x,t) = X(x)T(t)
void inlet_BC_func_T(const double &t, double &F)
{
    string inlet_function_type = "sinusoidal";

    if (inlet_function_type == "sinusoidal") {
        // For a sinusoidal input
        double a1 = 0.5;
        double b1 = a1 * M_PI * globalVars.epsilon * globalVars.epsilon * globalVars.BC_frequency_scale; // chosen to oscillate according to the maximum allowable amount (\epsilon^{-2})
        F = a1 - a1 * cos(t * M_PI / b1);
    }
    else if (inlet_function_type == "exponential") {
        // For a exponential input
        F = 1 - exp(-t / (globalVars.epsilon * globalVars.epsilon * globalVars.BC_frequency_scale));
    }
    else {
        cerr << globalVars.FILENAME << ": CRITICAL ERROR: Unrecognized inlet function." << endl;
        exit(1);
    }
}

// Define the space-dependent part of the inlet boundary condition c(x,t) = X(x)T(t)
void inlet_BC_func_X(const Vector &p, double &F)
{
   int dim = p.Size();
   //int vdim = F.Size();

   real_t x = p(0);
   real_t y = p(1);
   real_t z = (dim == 3) ? p(2) : 0.0;
   
   F = 1.0;
}
