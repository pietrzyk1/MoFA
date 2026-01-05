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



// Declare a class for saving the solutions
class ParSaveManager
{
private:
    int N_vdofs, rank, save_count = 0;

    ParGridFunction sol_gf_nonTrue;
    Vector sol_vec_nonTrue;

    string prefix, suffix, output_dir;
    bool IsSetPrefix = false, IsSetSuffix = false, IsSetOutputDir = false;


    // Define private class methods
    
    // Define a function for saving a solution vector in parallel (i.e., individual outputs for each rank)
    void SaveParSol()
    {
        ostringstream sol_name;
        sol_name << output_dir << prefix + to_string(save_count) + suffix << "." << setfill('0') << setw(6) << rank;
        ofstream sol_ofs(sol_name.str());
        sol_ofs.precision(8);
        sol_gf_nonTrue.Save(sol_ofs);
    }

    // Define a function for saving a solution vector in serial (i.e., one output that combines the solutions of all ranks)
    void SaveSerialSol(Mesh &mesh_serial)
    {
        // Have rank 0 save the entire serial/global solution as one gridfunction
        ofstream sol_ofs_asone((output_dir + prefix + to_string(save_count) + suffix).c_str()); //+ ".serial").c_str());
        GridFunction sol_serial = sol_gf_nonTrue.GetSerialGridFunction(0, mesh_serial);
        sol_serial.Save(sol_ofs_asone);
    }

public:
    // Define the class constructors
    ParSaveManager(ParFiniteElementSpace *fespace, int N_vdofs_, int rank_) : N_vdofs(N_vdofs_), rank(rank_)
    {
        sol_vec_nonTrue = Vector(N_vdofs);
        sol_gf_nonTrue = ParGridFunction(fespace);
        sol_gf_nonTrue.MakeRef(fespace, sol_vec_nonTrue, 0);
    }


    // Define functions for setting the variables of the class
    void SetPrefixAndSuffix(string &prefix_, string &suffix_) { SetPrefix(prefix_); SetSuffix(suffix_); }
    void SetPrefix(string &prefix_) { prefix = prefix_; IsSetPrefix = true; }
    void SetSuffix(string &suffix_) { suffix = suffix_; IsSetSuffix = true; }
    void SetOutputDir(string &output_dir_) { output_dir = output_dir_; IsSetOutputDir = true; }
    
    // Define a function for saving a provided solution vector as a gridfunction. The solution is saved by individual ranks, as well as combined to be saved as a serial solution
    void Save(Vector &tsol, Mesh &mesh_serial)
    {
        assert (IsSetPrefix && IsSetSuffix && IsSetOutputDir);
        sol_gf_nonTrue.Distribute(tsol);
        SaveParSol();
        SaveSerialSol(mesh_serial);
        save_count += 1;
    }
};


// Declare BC functions. They are defined after "main".
void inlet_BC_func_T(const double &t, double &F);
void inlet_BC_func_X(const Vector &p, double &F);

// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    string FILENAME = "Transport_Solver_p.cpp";
    vector<double> L = {1., 1., 0.};
    double epsilon = 0.1;
    double BC_frequency_scale = 1.0;
};
static Params globalVars;






int main(int argc, char *argv[])
{
    // ===============================================================
    //   Initialize MPI and HYPRE
    // ===============================================================
    int rank = 0;
    #ifdef MPI_BUILD
        Mpi::Init(argc, argv);
        //int N_ranks = Mpi::WorldSize(); // Don't need yet
        rank = Mpi::WorldRank();
        Hypre::Init();
    #endif

    // Print from each rank to the consol
    cout << globalVars.FILENAME << ": Hello from rank " << rank << "!" << endl;


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
    for (int i = 1; i < argc; i++) {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc) {
            config_path = argv[i + 1];
            if (rank == 0) { cout << globalVars.FILENAME << ": Configuration path obtained from parser options: " << config_path << endl; }
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
        cerr << globalVars.FILENAME << ": main(): CRITICAL ERROR: Could not load config file." << endl;
        exit(1);
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


    // Initialize variables from parser inputs
    if (active_advection == 1) { mesh_file_path = fluid_velocity_mesh_file_path; }
    
    
    // ===============================================================
    //   Get the mesh, the number of averaging regions, and the defined boundary attributes.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Obtaining mesh and mesh info... "; }

    // Get the boundary attributes for the mesh (NOTE: the provided attribute tags are 1 more than the actual Array index used to assign essential boundary conditions (SEE BELOW))
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_file_path);
    vector<int> AR_tags = (*mesh_info["AR"])["tags"];
    string avg_area_type = (*mesh_info["AR"])["area_type"];
    globalVars.L = (*mesh_info["geometry"])["L"];

    
    // Use the provided mesh file path to define the mesh
    MeshManager mesh_manager(mesh_file_path);
    
    // Make the mesh periodic if required by isPeriodic
    mesh_manager.MakePeriodic(isPeriodic, globalVars.L);
    
    // Make the parallel-partitioned mesh and use it
    mesh_manager.MakeParallelMesh(MPI_COMM_WORLD);

    // Get a pointer to the mesh being used
    ParMesh *mesh = mesh_manager.GetParMesh();

    if (rank == 0) {
        cout << "Complete." << endl;
        cout << globalVars.FILENAME << ":   Mesh is " << mesh->Dimension() << "D." << endl;
    }
    

    // ===============================================================
    //   Save the ParMeshes and a serial version of the global mesh.
    // ===============================================================
    {
        // Create the file names for each rank
        ostringstream mesh_name;
        mesh_name << mesh_output_file_name << "." << setfill('0') << setw(6) << rank;

        // Have each rank save their ParMesh
        ofstream mesh_ofs((output_dir + mesh_name.str()).c_str());
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);
    }

    // Have rank 0 also save the entire serial/global mesh for 
    Mesh mesh_serial = mesh->GetSerialMesh(0);
    ofstream mesh_out((output_dir + mesh_output_file_name).c_str()); //+ ".serial").c_str());
    mesh_serial.Print(mesh_out);


    // ===============================================================
    //   Define the finite element space for the concentration.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Defining finite element spaces... "; }

    // Finite element space for concentration
    FiniteElementCollection *fec_c = new H1_FECollection(order, mesh->Dimension());
    ParFiniteElementSpace *fespace_c = new ParFiniteElementSpace(mesh, fec_c, 1);
    
    if (rank == 0) {
        cout << "Complete." << endl;
        
        // Print the number of unknowns for concentration and total
        cout << globalVars.FILENAME << ":   Rank 0 number of concentration unknowns: " << fespace_c->GetTrueVSize() << endl;
        cout << globalVars.FILENAME << ":   Rank 0 number of unknowns in total: " << fespace_c->GetTrueVSize() << endl;
        //#ifdef MPI_BUILD
        //    cout << globalVars.FILENAME << ":   Global number of concentration unknowns: " << fespace_c->GlobalTrueVSize() << endl;
        //    cout << globalVars.FILENAME << ":   Global number of unknowns in total: " << fespace_c->GlobalTrueVSize() << endl;
        //#endif
    }


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
        if (rank == 0) { cout << globalVars.FILENAME << ":   Loading fluid velocity... "; }
        fluid_velocity_vgfc_manager = std::make_unique<VectorGridFunctionCoefficientManager>(fluid_velocity_file_path, fluid_velocity_mesh_file_path, Pe);
        if (rank == 0) { cout << "Complete." << endl; }
    
        // Create the fluid velocity VectorGridFunctionCoefficient
        if (rank == 0) { cout << globalVars.FILENAME << ":   Creating fluid velocity VectorGridFunctionCoefficient... "; }
        fluid_velocity_vgfc_manager->MakeParGridFunction(mesh, mesh_manager.GetPartitioning());
        fluid_velocity_vgfc_manager->MakeParVectorGridFunctionCoefficient();
        fluid_velocity = fluid_velocity_vgfc_manager->GetVectorGridFunctionCoefficient();
        if (rank == 0) { cout << "Complete." << endl; }
        
        /*
        // For debugging
        GridFunction u_gf_save(fluid_velocity_vgfc_manager->GetParGridFunctionFES());
        u_gf_save = 0.0;
        u_gf_save.ProjectCoefficient(*fluid_velocity);
        u_gf_save.Save(("fluid_velocity_verification_plot" + to_string(rank) + ".gf").c_str());
        
        ostringstream mesh_name2;
        mesh_name2 << mesh_output_file_name << "2." << setfill('0') << setw(6) << rank;
        ofstream mesh_ofs2(mesh_name2.str());
        mesh_ofs2.precision(8);
        mesh->Print(mesh_ofs2);
        exit(1);
        */
        /*
        // For debugging
        //GridFunction u_gf_save(fluid_velocity_vgfc_manager->GetGridFunctionFES());
        //u_gf_save = 0.0;
        //u_gf_save.ProjectCoefficient(*fluid_velocity);

        //int* partitioning2 = mesh_manager.GetPartitioning(); //fluid_velocity_vgfc_manager->GetGridFunctionMesh()->GeneratePartitioning(Mpi::WorldSize());
        //ParMesh par_mesh_test(MPI_COMM_WORLD, *fluid_velocity_vgfc_manager->GetGridFunctionMesh(), partitioning2);
        
        //ParGridFunction kp; //(&par_mesh_test, &u_gf_save, partitioning2); // mesh_manager.GetPartitioning());
        //kp = ParGridFunction(&par_mesh_test, &u_gf_save, partitioning2);
        //kp.Save(("fluid_velocity_verification_plot" + to_string(rank) + ".gf").c_str());
        
        
        //ostringstream mesh_name;
        //mesh_name << mesh_output_file_name << "." << setfill('0') << setw(6) << rank;
        //ofstream mesh_ofs(mesh_name.str());
        //mesh_ofs.precision(8);
        //fluid_velocity_vgfc_manager->GetGridFunctionMesh()->Print(mesh_ofs);
        //ParMesh par_mesh_test(MPI_COMM_WORLD, *fluid_velocity_vgfc_manager->GetGridFunctionMesh(), mesh_manager.GetPartitioning());
        //par_mesh_test.Print(mesh_ofs);
        */
    }
    

    // ===============================================================
    //   Define/Prepare boundary conditions.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Defining boundary condition functions... "; }

    // Here, we are defining two marker arrays to mark 1.) the inlet, and 2.) the outlet.
    Array<int> marker_inlet_BC(mesh->bdr_attributes.Max()); // Define a marker array for the inlet BC
    marker_inlet_BC = 0; // Initialize marker array to "don't apply any essential BC"
    marker_inlet_BC[(int)(*mesh_info["scalar_closure"])["inlet1"] - 1] = 1; // Set the first bdr attr group (i.e., physical group defined as 1 in gmsh file) as being the inlet
    ParDirichlettBC *inlet_BC = new ParDirichlettBC(marker_inlet_BC);
    inlet_BC->SetTemporallyDependentFunction( inlet_BC_func_T );
    inlet_BC->SetSpatiallyDependentFunction( inlet_BC_func_X );

    if (rank == 0) { cout << "Complete." << endl; }
    

    // ===============================================================
    //   Define block structure of the system.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Defining the block system structure... "; }

    // Notes: - This defines an array of offsets for each variable. The last component of the Array is the sum of the dimensions
    //          of each block.
    Array<int> block_offsets(2); // The argument should be the number of variables + 1
    block_offsets[0] = 0;
    block_offsets[1] = fespace_c->GetVSize();
    block_offsets.PartialSum();
    
    Array<int> block_offsets_true(2); // The argument should be the number of variables + 1
    block_offsets_true[0] = 0;
    block_offsets_true[1] = fespace_c->GetTrueVSize();
    block_offsets_true.PartialSum();

    if (rank == 0) { cout << "Complete." << endl; }
    

    // ===============================================================
    //   Define the solution and right-hand-side (block) vectors and initialize them.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Initializing the solution vector and RHS vector... "; }

    BlockVector sol_BLK(block_offsets), b_BLK(block_offsets);
    BlockVector sol_BLK_true(block_offsets_true), b_BLK_true(block_offsets_true);
    sol_BLK = 0.0; sol_BLK_true = 0.0;
    b_BLK = 0.0; b_BLK_true = 0.0;
    
    // Make references grid functions to the various dependent variables in the block solution. NOTE: Editing these gridfunctions will edit the block solutions
    ParGridFunction sol_BLK0_ref, b_BLK0_ref;
    sol_BLK0_ref.MakeRef(fespace_c, sol_BLK.GetBlock(0), 0);
    b_BLK0_ref.MakeRef(fespace_c, b_BLK.GetBlock(0), 0);
    
    if (rank == 0) { cout << "Complete." << endl; }


    // ===============================================================
    //   Create the bilinear forms.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Assembling bilinear forms, and stiffness and mass matrices... "; }

    // Create the bilinear form of the diffusion term
    ParBilinearForm *varf_cdiff(new ParBilinearForm(fespace_c));
    varf_cdiff->AddDomainIntegrator(new DiffusionIntegrator);
    if (active_advection == 1) { varf_cdiff->AddDomainIntegrator(new ConvectionIntegrator(*fluid_velocity)); }
    varf_cdiff->Assemble();
    varf_cdiff->Finalize();

    // Get the Hypre matrix of the bilinear form for diffusion, and set the inlet BC operator with the diffusion Hypre matrix
    HypreParMatrix *M_cdiff = varf_cdiff->ParallelAssemble();
    inlet_BC->SetOperator(*M_cdiff, *varf_cdiff);
    

    // Create the bilinear form of the mass matrix
    ParBilinearForm *varf_mass(new ParBilinearForm(fespace_c));
    ConstantCoefficient dimlessTimeScale(omega);
    varf_mass->AddDomainIntegrator(new MassIntegrator(dimlessTimeScale));
    varf_mass->Assemble();
    varf_mass->Finalize();

    // Get the Hypre matrix of the bilinear form for "mass", or the time derivative
    HypreParMatrix *M_mass = varf_mass->ParallelAssemble();
    
    if (rank == 0) { cout << "Complete." << endl; }
    
    
    // ===============================================================
    //   Solve the system in time.
    // ===============================================================
    double t = 0.0;
    int save_count = 0;
    
    // Define the time stepping operator for ODE systems of the form   M * du/dt + K * u = b, and set the initial time
    ParLinearTimeDependentOperator oper(*M_mass, *M_cdiff, b_BLK_true, dt, MPI_COMM_WORLD);
    oper.PrepareImplicitEuler();
    //oper.PrepareExplicitEuler();
    oper.SetTime(t);
    
    // Initialize block vector for hosting the solution
    BlockVector du_dt(block_offsets), delta_u(block_offsets);
    BlockVector du_dt_true(block_offsets_true), delta_u_true(block_offsets_true);
    du_dt = 0.0; delta_u = 0.0; du_dt_true = 0.0; delta_u_true = 0.0;


    // Initialize the manager for saving the solution in parallel
    ParSaveManager saver(fespace_c, sol_BLK.GetBlock(0).Size(), rank);
    saver.SetPrefixAndSuffix(output_file_name_prefix, output_file_name_suffix);
    saver.SetOutputDir(output_dir);
    
    // Save the initial condition in both parallel and serial formats
    saver.Save(sol_BLK_true.GetBlock(0), mesh_serial);

    
    // Start the time loop
    for (int i_step = 1; i_step < N_steps + 1; i_step++)
    {
        // Increase the time and print the step
        t += dt;
        if (rank == 0) { cout << "Time step: " << i_step << ", Time: " << t << endl; }        

        // Reset the RHS/b-vector
        b_BLK_true = 0.0;
        
        // Set the time on the inlet/source boundary condition
        inlet_BC->SetTime(t);
        
        // Project the boundary condition to the solution vector sol_BLK through sol_BLK0_ref
        inlet_BC->ProjectToGridFunction(sol_BLK0_ref);
        
        // Project the solution vector to the "true" solution vector, which will be used to solve the system
        fespace_c->GetRestrictionMatrix()->Mult(sol_BLK.GetBlock(0), sol_BLK_true.GetBlock(0));
        
        // Update the peripheral true vdofs in the b vector (i.e., the true vdofs that are affected by the Dirichlett BC). Note, this ADDS the contribution to b_BLK_true; it does not rewrite the entries
        inlet_BC->UpdateLinearFormVector(sol_BLK_true.GetBlock(0), b_BLK_true.GetBlock(0));

        // Apply operator to solve for the time derivative
        //oper.Mult(sol_BLK, du_dt);
        oper.ImplicitSolve(dt, sol_BLK_true, du_dt_true);
        
        // Multiply derivative by dt and add to solution to get the solution at the next time step
        delta_u_true = du_dt_true;
        delta_u_true *= dt;
        sol_BLK_true.GetBlock(0) += delta_u_true;

        // The solution is obtained for the true-DOF-representation. This needs to be distributed to the non-true-DOF-representation on each rank
        sol_BLK0_ref.Distribute(&(sol_BLK_true.GetBlock(0)));


        // Save the solution if required
        if (i_step % output_interval == 0) {
            // Save the pore-scale solution in parallel and serial formats
            saver.Save(sol_BLK_true.GetBlock(0), mesh_serial);
        }


        // Print the value of the BC's time-dependent portion to the consol (for debugging purposes)
        if (rank == 0) {
            double BC_val;
            inlet_BC_func_T(t, BC_val);
            cout << "BC val: " << BC_val << endl;
        }
    }


    // ===============================================================
    //   Free the used memory by deleting the pointers.
    // ===============================================================
    delete varf_cdiff;
    //delete IC;
    delete inlet_BC;
    delete fespace_c;
    delete fec_c;

    
    // ===============================================================
    //   Finalize the MPI and exit.
    // ===============================================================
    MPI_Finalize();
    
    return 0;
}










// Define the time-dependent part of the inlet boundary condition c(x,t) = X(x)T(t)
void inlet_BC_func_T(const double &t, double &F)
{
    double a1 = 0.5;
    double b1 = a1 * M_PI * globalVars.epsilon * globalVars.epsilon * globalVars.BC_frequency_scale; //* 0.75; // chosen to oscillate according to the maximum allowable amount (\epsilon^{-2})
    F = a1 - a1 * cos(t * M_PI / b1);
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
