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

//                                      Stokes Solver
//
// Description: This code solves the dimensionless, incompressible Stokes equations
//              in MFEM. These equations are described as the following:
//
//              Momentum:
//                                A * div(grad(u)) - grad(p) = F
//
//              Mass Continuity:
//                                        -div(u) = 0
//
//              (Note: the negative allows for the same block matrix to be used for
//                     -div(u) and -grad(p), which makes the block system symmetric
//                     for the MINRES iterative solver)
//
//              where F is the body force vector. A velocity Dirichlet condition can
//              be assigned at the labeled inlet, and p = 0 is assigned to the outlet.
//              A no-slip condition is used on all other boundaries. 2D/3D domains
//              can be handled. Periodic boundary conditions and can also be handled.
//              These work best when a constant body force is applied to the momentum
//              equation. 
//

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include "JSON_IO.h"
#include "mfem_util.h"



using namespace std;
using namespace mfem;



// Declare BC functions for velocity. They are defined after "main".
void u_inlet_func(const Vector &p, Vector &F);
void u_no_slip_func(const Vector &p, Vector &F);

// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    string FILENAME = "Stokes_Solver_p.cpp";
    vector<double> L = {1., 1., 0.};
    double u_inlet_max = 1.0;
    string sim_type = "";
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
    cout << globalVars.FILENAME << ":   Hello from rank " << rank << "!" << endl;


    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string config_path = "./";

    string mesh_dir = "./";
    string mesh_file_name = "mesh.msh";
    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";

    string output_dir = "./";
    string u_output_file_name = "u_sol.gf";
    string p_output_file_name = "p_sol.gf";
    string mesh_output_file_name = "mesh.mesh";
    
    int u_order = 2;
    int p_order = 1;
    vector<int> isPeriodic = {0, 0, 0};
    int useInlet = 0;
    
    double A_coef_val = 1.0;
    vector<double> const_body_force = {0.0, 0.0};
    
    int maxIter = 10000;
    double rtol = 1.0e-7;
    double atol = 0.0;

    
    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++) {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc) {
            config_path = argv[i + 1];
            if (rank == 0) { cout << globalVars.FILENAME << ":   Configuration path obtained from parser options: " << config_path << endl; }
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
        
        
        JSONDict stokes_dict = *configData["stokes"];
        
        sub_dict = *stokes_dict["output path"];
        sub_dict.getValue("directory", output_dir);
        sub_dict.getValue("velocity file name", u_output_file_name);
        sub_dict.getValue("pressure file name", p_output_file_name);
        sub_dict.getValue("mesh file name", mesh_output_file_name);
        
        sub_dict = *stokes_dict["simulation parameters"];
        sub_dict.getValue("u order", u_order);
        sub_dict.getValue("p order", p_order);
        sub_dict.getValue("isPeriodic", isPeriodic);
        sub_dict.getValue("use inlet", useInlet);

        sub_dict = *stokes_dict["physics parameters"];
        sub_dict.getValue("A", A_coef_val);
        sub_dict.getValue("body force", const_body_force);
        sub_dict.getValue("u inlet max", globalVars.u_inlet_max);

        sub_dict = *stokes_dict["solver parameters"];
        sub_dict.getValue("max iterations", maxIter);
        sub_dict.getValue("rel tol", rtol);
        sub_dict.getValue("abs tol", atol);
    }
    
    
    // ===============================================================
    //   Parse command line options.
    // ===============================================================
    // Parse options
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    args.AddOption(&mesh_dir, "-M", "--mesh_dir", "Mesh-containing directory.");
    args.AddOption(&mesh_file_name, "-m", "--mesh_file_name", "File name to use for the mesh.");
    args.AddOption(&mesh_info_file_name, "-i", "--mesh_info_file_name", "File name to use for the mesh boundary attributes (assumed to be in the mesh directory).");
    args.AddOption(&output_dir, "-O", "--output_dir", "Output directory.");
    args.AddOption(&u_order, "-u", "--order", "Finite element order (polynomial degree) of velocity (pressure is one less).");
    args.AddOption(&A_coef_val, "-A", "--A_coef", "Value of the dimensionless number A.");
    args.ParseCheck();

    
    // ===============================================================
    //   Define variables based on the options provided in the parser
    // ===============================================================
    // Define output velocity, pressure, and mesh save paths
    string u_output_path = output_dir + u_output_file_name;
    string p_output_path = output_dir + p_output_file_name;
    string mesh_output_path = output_dir + mesh_output_file_name;
    
    // Ensure the element order of pressure is one less than that of the velocity
    if (p_order != u_order - 1) { p_order = u_order - 1; }
    
    // Enable hardware devices, such as GPUs, and programming models like CUDA, OCCA, RAJA, and OpenMP
    //Device device("cpu");
    //if (Mpi::Root()) { device.Print(); }

    
    // ===============================================================
    //   Get the mesh and setup periodicity.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Obtaining mesh and mesh info... "; }

    // Get mesh info and load the necessary details
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_dir + mesh_info_file_name);
    globalVars.L = (*mesh_info["geometry"])["L"];

    // Use the provided mesh file path to define the mesh 
    MeshManager mesh_manager(mesh_dir + mesh_file_name);
    
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
    //   Define finite element spaces for the velocity and pressure.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Defining finite element spaces... "; }

    // Finite element space for velocity
    FiniteElementCollection *fec_u = new H1_FECollection(u_order, mesh->Dimension()); // People have also used RT_FECollection
    ParFiniteElementSpace *fespace_u = new ParFiniteElementSpace(mesh, fec_u, mesh->Dimension());
    
    // Finite element space for pressure
    FiniteElementCollection *fec_p = new H1_FECollection(p_order, mesh->Dimension()); // People have also used L2_FECollection
    ParFiniteElementSpace *fespace_p = new ParFiniteElementSpace(mesh, fec_p, 1);

    if (rank == 0) { cout << "Complete." << endl; }


    // Print the number of unknowns for velocity, pressure, and total
    if (rank == 0) {
        cout << globalVars.FILENAME << ":   Rank 0 number of velocity unknowns: " << fespace_u->GetTrueVSize() << endl;
        cout << globalVars.FILENAME << ":   Rank 0 number of pressure unknowns: " << fespace_p->GetTrueVSize() << endl;
        cout << globalVars.FILENAME << ":   Rank 0 number of unknowns in total: " << fespace_u->GetTrueVSize() + fespace_p->GetTrueVSize() << endl;
        //#ifdef MPI_BUILD
        //    cout << globalVars.FILENAME << ":   Global number of velocity unknowns: " << fespace_u->GlobalTrueVSize() << endl;
        //    cout << globalVars.FILENAME << ":   Global number of pressure unknowns: " << fespace_p->GlobalTrueVSize() << endl;
        //    cout << globalVars.FILENAME << ":   Global number of unknowns in total: " << fespace_u->GlobalTrueVSize() + fespace_p->GlobalTrueVSize() << endl;
        //#endif
    }
    
    
    // ===============================================================
    //   Define/Prepare boundary conditions.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Defining boundary condition functions... "; }

    // Here, we are defining two marker arrays to mark 1.) the inlet, and 2.) the surfaces where the no-slip BC is applied.
    Array<int> marker_inlet_BC(mesh->bdr_attributes.Max()); // Define a marker array for the inlet BC.
    marker_inlet_BC = 0; // Initialize marker array to "don't apply any essential BC".
    if (useInlet == 1) { marker_inlet_BC[(int)(*mesh_info["stokes"])["inlet1"] - 1] = 1; } // Set the first bdr attr group (i.e., physical group defined as 1 in gmsh file) as being the inlet.
    VectorFunctionCoefficient inlet_BC(mesh->SpaceDimension(), u_inlet_func); // Define a function coefficient that will apply "u_inlet_func", which varies in space.
    
    Array<int> marker_no_slip_BC(mesh->bdr_attributes.Max()); // Define a marker array for the no-slip BC.
    marker_no_slip_BC = 0; // Initialize marker array to "don't apply any essential BC".
    marker_no_slip_BC[(int)(*mesh_info["stokes"])["noslip1"] - 1] = 1; // Set the third bdr attr group (i.e., physical group defined as 3 in gmsh file) as being the no-slip condition.
    VectorFunctionCoefficient no_slip_BC(mesh->SpaceDimension(), u_no_slip_func); // Define a function coefficient that will apply "u_inlet_func", which varies in space.

    if (rank == 0) { cout << "Complete." << endl; }


    // ===============================================================
    //   Define block structure of the system.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Defining the block system structure... "; }

    // Notes: - This defines an array of offsets for each variable. The last component of the Array is the sum of the dimensions
    //          of each block.
    Array<int> block_offsets(3); // The argument should be the number of variables + 1
    block_offsets[0] = 0;
    block_offsets[1] = fespace_u->GetVSize();
    block_offsets[2] = fespace_p->GetVSize();
    block_offsets.PartialSum();
    
    Array<int> block_offsets_true(3); // The argument should be the number of variables + 1
    block_offsets_true[0] = 0;
    block_offsets_true[1] = fespace_u->GetTrueVSize();
    block_offsets_true[2] = fespace_p->GetTrueVSize();
    block_offsets_true.PartialSum();

    if (rank == 0) { cout << "Complete." << endl; }
    

    // ===============================================================
    //   Define the solution and right-hand-side (block) vectors and initialize them.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Initializing the solution vector and RHS vector... "; }

    BlockVector sol_up_BLK(block_offsets), b_BLK(block_offsets);
    BlockVector sol_up_BLK_true(block_offsets_true), b_BLK_true(block_offsets_true);
    sol_up_BLK = 0.0, sol_up_BLK_true = 0.0;
    b_BLK = 0.0, b_BLK_true = 0.0;

    if (rank == 0) { cout << "Complete." << endl; }

    
    // ===============================================================
    //   Create the linear form.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Defining the RHS vector (body force function)... "; }

    ParLinearForm varf_body_forc(fespace_u);
    Vector body_force_vec(const_body_force.data(), const_body_force.size());
    VectorConstantCoefficient body_force(body_force_vec);
    varf_body_forc.AddDomainIntegrator(new VectorDomainLFIntegrator(body_force));
    varf_body_forc.Assemble();
    
    // Set the RHS blocks using the linear form
    b_BLK.GetBlock(0) = varf_body_forc;

    if (rank == 0) { cout << "Complete." << endl; }
        

    // ===============================================================
    //   Create the bilinear forms.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Assembling bilinear forms and stiffness matrix... "; }

    // Notes: - We only create the bilinear form for the viscous term and the incompressible continuity equation. This is because the
    //          transpose of the incompressible continuity equation matrix operator will be identical to the gradient of pressure
    //          operator in the momentum equation.
    //
    // Initiate the bilinear form of the viscous term
    ParBilinearForm varf_viscous(fespace_u);
    ConstantCoefficient A_coef(A_coef_val);
    varf_viscous.AddDomainIntegrator(new VectorDiffusionIntegrator(A_coef));
    varf_viscous.Assemble();
    varf_viscous.Finalize();
    
    // Initiate the bilinear form of the incompressible continuity equation
    ParMixedBilinearForm varf_incomp(fespace_u, fespace_p);
    ConstantCoefficient neg_one(-1);
    varf_incomp.AddDomainIntegrator(new VectorDivergenceIntegrator(neg_one));
    varf_incomp.Assemble();
    varf_incomp.Finalize();
    
    // Initiate pointers to the HyprePar matrices of the bilinear forms
    HypreParMatrix *M_viscous = varf_viscous.ParallelAssemble();
    HypreParMatrix *M_incomp = varf_incomp.ParallelAssemble();

    // Apply the essential (Dirichlet) boundary conditions (i.e., the inlet and no-slip condition). The order matters: no-slip takes priority at points defined at the BC intersections
    ParGridFunction u_BC_apply;
    u_BC_apply.MakeRef(fespace_u, sol_up_BLK.GetBlock(0), 0);
    if (useInlet == 1) { u_BC_apply.ProjectBdrCoefficient(inlet_BC, marker_inlet_BC); }
    u_BC_apply.ProjectBdrCoefficient(no_slip_BC, marker_no_slip_BC);
    

    // The essential (Dirichlet) boundary conditions and body forces were applied to the non-true vdofs; however, hypre entities use true vdofs.
    // As such, move the solution and b vectors to the corresponding true vdofs vectors, and continue using true vdofs (hypre uses true vdofs)
    fespace_u->GetRestrictionMatrix()->Mult(sol_up_BLK.GetBlock(0), sol_up_BLK_true.GetBlock(0));
    fespace_p->GetRestrictionMatrix()->Mult(sol_up_BLK.GetBlock(1), sol_up_BLK_true.GetBlock(1));
    fespace_u->GetRestrictionMatrix()->Mult(b_BLK.GetBlock(0), b_BLK_true.GetBlock(0));
    fespace_p->GetRestrictionMatrix()->Mult(b_BLK.GetBlock(1), b_BLK_true.GetBlock(1));

    
    // Apply the inlet BC to both bilinear forms
    if (useInlet == 1) {
        // Apply inlet BC to varf_viscous
        HypreParMatrix *M_viscous_BC_elim = varf_viscous.ParallelEliminateEssentialBC(marker_inlet_BC, *M_viscous); // separate-out the rows and columns (into M_viscous_BC_elim) of M_viscous that correspond to the true vdofs at the inlet
        Vector b_BLK0_dummy(b_BLK_true.GetBlock(0).Size()); b_BLK0_dummy = 0.0;
        Array<int> tdof_list; fespace_u->GetEssentialTrueDofs(marker_inlet_BC, tdof_list); // get the true vdofs at the inlet
        EliminateBC(*M_viscous, *M_viscous_BC_elim, tdof_list, sol_up_BLK_true.GetBlock(0), b_BLK0_dummy); // move the contributions of the true vdofs at the inlet into the b vector
        b_BLK_true.GetBlock(0) += b_BLK0_dummy;
        delete M_viscous_BC_elim;
        
        HypreParMatrix *M_incomp_col_elim = M_incomp->EliminateCols(tdof_list); // separate-out the columns (into M_incomp_col_elim) of M_incomp that correspond to the true vdofs at the inlet
        Vector b_BLK1_dummy(b_BLK_true.GetBlock(1).Size()); b_BLK1_dummy = 0.0;
        M_incomp_col_elim->Mult(sol_up_BLK_true.GetBlock(0), b_BLK1_dummy); // multiply M_incomp_col_elim by the velocity solution (which has the inlet projected onto it) to get the contributions to the b vector
        b_BLK_true.GetBlock(1) -= b_BLK1_dummy; // move the contributions to the rhs (subtract) and combine with b_BLK_true
        delete M_incomp_col_elim;
    }

    // Apply no slip to varf_viscous
    HypreParMatrix *M_viscous_NS_elim = varf_viscous.ParallelEliminateEssentialBC(marker_no_slip_BC, *M_viscous); // separate-out the rows and columns (into M_viscous_BC_elim) of M_viscous that correspond to the true vdofs at the inlet
    Vector b_BLK0_dummy(b_BLK_true.GetBlock(0).Size()); b_BLK0_dummy = 0.0;
    Array<int> tdof_list; fespace_u->GetEssentialTrueDofs(marker_no_slip_BC, tdof_list); // get the true vdofs at the inlet
    EliminateBC(*M_viscous, *M_viscous_NS_elim, tdof_list, sol_up_BLK_true.GetBlock(0), b_BLK0_dummy); // move the contributions of the true vdofs at the inlet into the b vector
    b_BLK_true.GetBlock(0) += b_BLK0_dummy;
    delete M_viscous_NS_elim;

    HypreParMatrix *M_incomp_col_elim = M_incomp->EliminateCols(tdof_list); // separate-out the columns (into M_incomp_col_elim) of M_incomp that correspond to the true vdofs at the inlet
    Vector b_BLK1_dummy(b_BLK_true.GetBlock(1).Size()); b_BLK1_dummy = 0.0;
    M_incomp_col_elim->Mult(sol_up_BLK_true.GetBlock(0), b_BLK1_dummy); // multiply M_incomp_col_elim by the velocity solution (which has the inlet projected onto it) to get the contributions to the b vector
    b_BLK_true.GetBlock(1) -= b_BLK1_dummy; // move the contributions to the rhs (subtract) and combine with b_BLK_true
    delete M_incomp_col_elim;
    
    if (rank == 0) { cout << "Complete." << endl; }
    
    
    // ===============================================================
    //   Define the block operator for the system.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Assembling the block system... "; }
    
    // Initiate block operator
    BlockOperator stokesOp_BLK(block_offsets_true);

    // Create the transpose of the continuity equation's bilinear form (this will be the pressure gradient's block matrix for the momentum equations)
    TransposeOperator *M_incomp_T = NULL;
    M_incomp_T = new TransposeOperator(M_incomp);

    // Set the blocks of the Stokes/incompressible continuity equation block operator
    stokesOp_BLK.SetBlock(0, 0, M_viscous);
    stokesOp_BLK.SetBlock(0, 1, M_incomp_T);
    stokesOp_BLK.SetBlock(1, 0, M_incomp);

    if (rank == 0) { cout << "Complete." << endl; }
    

    // ===============================================================
    //   Construct the preconditioner operator.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Assembling the preconditioners... "; }
    
    // Notes: - Here, we use Symmetric Gauss-Seidel to approximate the inverse of the pressure Schur complement.
    //
    //                      P = [ BM_viscous      0    ] = [ BM_viscous                    0                   ]
    //                          [      0      S_actual ] = [     0       BM_incomp BM_viscous^(-1) BM_incomp_T ]
    //
    //    We use this ->      ~ [ BM_viscous  0 ] = [ BM_viscous                        0                     ]
    //                          [      0      S ]   [     0       BM_incomp diag(BM_viscous)^(-1) BM_incomp_T ]
    //
    //SaddlePointBlockPreconditioner stokesPC(block_offsets);
    //stokesPC.BuildPreconditioner(BM_viscous, BM_incomp);
    
    HypreParMatrix *MinvBt = NULL;
    HypreParVector *Md = NULL;
    HypreParMatrix *S = NULL;
    Vector Md_PA;
    Solver *invM, *invS;
    
    Md = new HypreParVector(MPI_COMM_WORLD, M_viscous->GetGlobalNumRows(), M_viscous->GetRowStarts());
    M_viscous->GetDiag(*Md);

    MinvBt = M_incomp->Transpose();
    MinvBt->InvScaleRows(*Md);
    S = ParMult(M_incomp, MinvBt);

    invM = new HypreDiagScale(*M_viscous);
    invS = new HypreSmoother(*S, 6); // 6 is for the GSSmoother

    invM->iterative_mode = false;
    invS->iterative_mode = false;

    BlockDiagonalPreconditioner *stokesPC = new BlockDiagonalPreconditioner(block_offsets_true);
    stokesPC->SetDiagonalBlock(0, invM);
    stokesPC->SetDiagonalBlock(1, invS);

    if (rank == 0) { cout << "Complete." << endl; }

    
    // ===============================================================
    //   Solve the system with MINRES.
    // ===============================================================
    if (rank == 0) { cout << globalVars.FILENAME << ":   Solving" << endl; }

    MINRESSolver solver(MPI_COMM_WORLD);
    solver.SetAbsTol(atol);
    solver.SetRelTol(rtol);
    solver.SetMaxIter(maxIter);
    solver.SetOperator(stokesOp_BLK);
    solver.SetPreconditioner(*stokesPC);
    solver.SetPrintLevel(1);
    solver.Mult(b_BLK_true, sol_up_BLK_true);
    
    if (rank == 0) {
        if (solver.GetConverged()) {
            cout << "MINRES converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << ".\n";
        }
        else {
            cout << "MINRES did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
                << ".\n";
        }
    }


    // ===============================================================
    //   Extract and save the solutions and mesh.
    // ===============================================================
    ParGridFunction u_sol, p_sol;
    
    // Distribute the velocity true DOFs to all processors
    u_sol.MakeRef(fespace_u, sol_up_BLK.GetBlock(0), 0); // In GLVis, you can use "u" to cycle through x- and y-components
    u_sol.Distribute(&(sol_up_BLK_true.GetBlock(0)));
    
    // Distribute the pressure true DOFs to all processors
    p_sol.MakeRef(fespace_p, sol_up_BLK.GetBlock(1), 0);
    p_sol.Distribute(&(sol_up_BLK_true.GetBlock(1)));

    {
        // Create the file names for each rank
        ostringstream mesh_name, u_name, p_name;
        mesh_name << mesh_output_file_name << "." << setfill('0') << setw(6) << rank;
        u_name << u_output_file_name << "." << setfill('0') << setw(6) << rank;
        p_name << p_output_file_name << "." << setfill('0') << setw(6) << rank;

        // Have each rank save their ParMesh
        ofstream mesh_ofs((output_dir + mesh_name.str()).c_str());
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);
        // Have rank 0 also save the entire serial/global mesh for 
        Mesh mesh_serial = mesh->GetSerialMesh(0);
        ofstream mesh_out((output_dir + mesh_output_file_name).c_str()); //+ ".serial").c_str());
        mesh_serial.Print(mesh_out);
        
        // Save the velocity as a vector field
        ofstream u_ofs((output_dir + u_name.str()).c_str());
        u_ofs.precision(8);
        u_sol.Save(u_ofs);
        // Have rank 0 also save the entire serial/global velocity solution as one gridfunction 
        ofstream u_ofs_asone((output_dir + u_output_file_name).c_str()); //+ ".serial").c_str());
        GridFunction u_sol_serial = u_sol.GetSerialGridFunction(0, mesh_serial);
        u_sol_serial.Save(u_ofs_asone);
        
        // Save the pressure field
        ofstream p_ofs((output_dir + p_name.str()).c_str());
        p_ofs.precision(8);
        p_sol.Save(p_ofs);
    }


    // ===============================================================
    //   Free the used memory.
    // ===============================================================
    delete M_incomp_T;
    delete fespace_p;
    delete fec_p;
    delete fespace_u;
    delete fec_u;


    // ===============================================================
    //   Finalize the MPI and exit.
    // ===============================================================
    MPI_Finalize();

    return 0;
}





// Define the inlet velocity condition
void u_inlet_func(const Vector &p, Vector &F)
{
   int dim = p.Size();

   real_t x = p(0);
   real_t y = p(1);
   real_t z = (dim == 3) ? p(2) : 0.0;
   
   double a = globalVars.L[1]/2;
   double b = globalVars.L[2]/2;

   F(0) = (dim == 3) ? -globalVars.u_inlet_max/(a*a*b*b) * (y - a) * (y + a) * (z - b) * (z + b)
                     : -globalVars.u_inlet_max/(a*a) * (y - a) * (y + a);
   F(1) = 0.0;
   if (dim == 3) { F(2) = 0.0; }
}

// Define the no-slip velocity condition
void u_no_slip_func(const Vector &p, Vector &F)
{
   int dim = p.Size();

   real_t x = p(0);
   real_t y = p(1);
   real_t z = (dim == 3) ? p(2) : 0.0;
   
   F(0) = 0.0;
   F(1) = 0.0;
   if (dim == 3) { F(2) = 0.0; }
}

