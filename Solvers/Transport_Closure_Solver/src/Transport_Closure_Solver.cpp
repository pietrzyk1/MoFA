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



using namespace std;
using namespace mfem;



// Define functions that are used in "main", but are defined after "main".
double inlet_BC_func(const Vector &p, int inlet_value_toggle);
//void outlet_BC_func(const Vector &p, Vector &F);

// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    vector<double> L = {1., 1., 0.};
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

    int is3D = 0;
    
    int order = 2;
    int active_advection = 1;
    vector<int> isPeriodic = {0, 0, 0};
    int useInlet = 0;
    
    int active_AR = -1;
    int active_inlet = 0;
    
    double Pe_s = 1.0;
    double omega = 1.0;

    int maxIter = 100000;
    double rtol = 1.0e-10;
    double atol = 0.0;

    int recursive_iter = 0;

    string output_dir = "./"; // For general output; this is adjusted by the config file for the different variations
    
    unique_ptr<string> fluid_velocity_dir = make_unique<string>(output_dir);
    string fluid_velocity_file_name = "u_sol.gf";

    unique_ptr<string> residual_output_dir = make_unique<string>(output_dir);
    string residual_output_file_name = "a_sol.txt";
    
    unique_ptr<string> closure_output_dir = make_unique<string>(output_dir);
    string closure_output_file_name = "c_sol.gf";
    
    unique_ptr<string> mesh_output_dir = make_unique<string>(output_dir);
    string mesh_output_file_name = "mesh.mesh";

    unique_ptr<string> forcing_function_dir = make_unique<string>(output_dir);
    string forcing_function_file_name = "None";
    
    unique_ptr<string> forcing_function2_dir = make_unique<string>(output_dir);
    string forcing_function2_file_name = "None";


    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++)
    {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc)
        {
            config_path = argv[i + 1];
            cout << "Transport_Closure_Solver.cpp: Configuration path obtained from parser options: " << config_path << endl;
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

        sub_dict = *mesh_dict["details"];
        sub_dict.getValue("is3D", is3D);

        
        JSONDict stokes_dict = *configData["stokes"];
        
        sub_dict = *stokes_dict["output path"];
        sub_dict.getValue("directory", fluid_velocity_dir);
        sub_dict.getValue("velocity file name", fluid_velocity_file_name);


        JSONDict closure_dict = *configData["scalar closure"];

        sub_dict = *closure_dict["simulation parameters"];
        sub_dict.getValue("order", order);
        sub_dict.getValue("active advection", active_advection);
        sub_dict.getValue("isPeriodic", isPeriodic);
        sub_dict.getValue("use inlet", useInlet);
        
        sub_dict = *closure_dict["closure parameters"];
        sub_dict.getValue("active averaging region", active_AR);
        sub_dict.getValue("active inlet", active_inlet);
        
        sub_dict = *closure_dict["physics parameters"];
        sub_dict.getValue("Pe_s", Pe_s);
        sub_dict.getValue("omega", omega);

        sub_dict = *closure_dict["solver parameters"];
        sub_dict.getValue("max iterations", maxIter);
        sub_dict.getValue("rel tol", rtol);
        sub_dict.getValue("abs tol", atol);

        sub_dict = *closure_dict["residual path"];
        sub_dict.getValue("directory", residual_output_dir);
        sub_dict.getValue("file name", residual_output_file_name);

        sub_dict = *closure_dict["closure path"];
        sub_dict.getValue("directory", closure_output_dir);
        sub_dict.getValue("file name", closure_output_file_name);
        
        sub_dict = *closure_dict["mesh output path"];
        sub_dict.getValue("directory", mesh_output_dir);
        sub_dict.getValue("file name", mesh_output_file_name);

        sub_dict = *closure_dict["forcing function path"];
        sub_dict.getValue("directory", forcing_function_dir);
        sub_dict.getValue("file name", forcing_function_file_name);
        
        //sub_dict = *closure_dict["forcing function 2 path"];
        //sub_dict.getValue("directory", forcing_function2_dir);
        //sub_dict.getValue("file name", forcing_function2_file_name);
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_file_path = mesh_dir + mesh_file_name;
    string mesh_info_file_path = *mesh_info_dir + mesh_info_file_name;
    string fluid_velocity_file_path = *fluid_velocity_dir + fluid_velocity_file_name;
    string residual_output_file_path = *residual_output_dir + residual_output_file_name;
    
    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    args.AddOption(&mesh_file_path, "-M", "--mesh_file_path", "Mesh file path (best to define in the config file).");
    args.AddOption(&mesh_info_file_path, "-m", "--mesh_info_file_path", "Mesh info file path (best to define in the config file).");
    //args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree) of scalar closure variable.");
    args.AddOption(&output_dir, "-O", "--output_dir", "Output directory for closure, residual, and mesh solutions. The specified forcing function will also be looked for in this directory (use config file for more detailed directories).");
    args.AddOption(&active_AR, "-a", "--active_AR", "The averaging region index where the closure variable averages to 1.");
    args.AddOption(&active_inlet, "-i", "--active_inlet", "Toggle whether the inlet dirichlett condition should be 1 or 0.");
    args.AddOption(&fluid_velocity_file_path, "-v", "--fluid_velocity_file_path", "Fluid velocity solution file path (best to define in the config file).");
    args.AddOption(&active_advection, "-A", "--active_advection", "Toggle whether to consider advection or not.");
    args.AddOption(&Pe_s, "-P", "--Peclet_Number", "The O(1) Peclet Number to use (i.e., this is Pe*\epsilon).");
    args.AddOption(&closure_output_file_name, "-N", "--closure_output_file_name", "Declare the file name of the output file for the closure variable solution (Note: the file name should end in '.gf').");
    args.AddOption(&forcing_function_file_name, "-F", "--forcing_function_file_name", "Declare the file name of the forcing function (Note: the file name should end in '.gf').");
    //args.AddOption(&forcing_function2_file_name, "-f", "--forcing_function2_file_name", "Declare the file name of the second forcing function (Note: the file name should end in '.gf').");
    args.AddOption(&recursive_iter, "-I", "--recursive_iter", "Provide the int that describes the iteration of this run.");
    args.ParseCheck();


    // ===============================================================
    //   Get the mesh, the number of averaging regions, and the defined boundary attributes.
    // ===============================================================
    cout << "TCS:   Obtaining mesh and mesh info... ";

    // Use the provided mesh file to define the mesh 
    Mesh *mesh = new Mesh(mesh_file_path);
    //mesh->UniformRefinement(); // Refine the mesh once uniformly

    // Get the number of averaging regions defined in the mesh file
    int N_AR = mesh->attributes.Max();

    // Get the boundary attributes for the mesh (NOTE: the provided attribute tags are 1 more than the actual Array index used to assign essential boundary conditions (SEE BELOW))
    JSONDict mesh_info;
    mesh_info.loadFromFile(mesh_info_file_path);
    globalVars.L = (*mesh_info["geometry"])["L"];

    // Create a periodic mesh if required
    Mesh periodic_mesh;
    vector<Vector> translations = {};
    Vector translation;
    for (int i = 0; i < mesh->Dimension(); i++)
    {
        if (isPeriodic[i])
        {
            translation = Vector({0.0, 0.0, 0.0});
            translation[i] = globalVars.L[i];
            translations.push_back( translation );
        }
    }
    if (translations.size() > 0)
    {
        periodic_mesh = Mesh::MakePeriodic(*mesh, mesh->CreatePeriodicVertexMapping(translations));
        delete mesh;
        mesh = &periodic_mesh;
    }

    cout << "Complete." << endl;
    cout << "Number of averaging regions: " << N_AR << endl;


    // ===============================================================
    //   Define finite element spaces for the concentration.
    // ===============================================================
    cout << "TCS:   Defining finite element spaces... ";

    // Finite element space for concentration
    FiniteElementCollection *fec_c = new H1_FECollection(order, mesh->Dimension());
    FiniteElementSpace *fespace_c = new FiniteElementSpace(mesh, fec_c, 1);
    
    cout << "Complete." << endl;
    
    // Print the number of unknowns for concentration and total
    cout << "Number of unknowns in concentration: " << fespace_c->GetTrueVSize() << endl;
    cout << "Number of unknowns in total: " << fespace_c->GetTrueVSize() << endl;

    
    // ===============================================================
    //   Load Advection
    // ===============================================================
    // Declare a grid function for the fluid velocity
    GridFunction *u_gf;
        
    // Use a previously saved gridfunction for the velocity field, or define the field as zero everywhere
    if (active_advection == 1)
    {
        cout << "TCS:   Loading fluid velocity... ";
    
        // Create an ifgzstream pointer to read in the file containing the fluid velocity
        ifgzstream *u_ifgz_stream = NULL;
        u_ifgz_stream = new ifgzstream(fluid_velocity_file_path);
        if (!(*u_ifgz_stream))
        {
            cerr << "Transport_Closure_Solver.cpp: CRITICAL ERROR: Cannot open fluid velocity file.\n";
            exit(1);
        }

        // Define the grid function with the fluid velocity data in the stream
        u_gf = new GridFunction(mesh, *u_ifgz_stream);

        // Scale the grid function by the Peclet number for the advection term
        *u_gf *= Pe_s;

        cout << "Complete." << endl;
    }
    else
    {
        // Create a finite element space on which to define the grid function for the fluid velocity
        FiniteElementCollection *fec_u = new H1_FECollection(order, mesh->Dimension());
        FiniteElementSpace *fespace_u = new FiniteElementSpace(mesh, fec_u, mesh->Dimension());
        u_gf = new GridFunction(fespace_u);
        
        // Initialize the fluid velocity as zero everywhere
        *u_gf = 0.0;
    }
    
    // Create a VectorGridFunctionCoefficient from the gridfunction for the advection term
    VectorGridFunctionCoefficient fluid_velocity(u_gf);


    // ===============================================================
    //   Load the forcing function
    // ===============================================================
    // Declare a grid function for the forcing function
    GridFunction *force_gf;
    
    // Create a finite element space on which to define the grid function for the forcing function (and later, the linear form)
    FiniteElementCollection *fec_f = new H1_FECollection(order, mesh->Dimension());
    FiniteElementSpace *fespace_f = new FiniteElementSpace(mesh, fec_f, 1);
    
    // Use a previously saved gridfunction for the forcing function, or define the field as zero everywhere
    if (forcing_function_file_name != "None")
    {
        cout << "TCS:   Loading forcing function... ";

        // Create an ifgzstream pointer to read in the file containing the forcing function
        ifgzstream *force_ifgz_stream = NULL;
        force_ifgz_stream = new ifgzstream(*forcing_function_dir + forcing_function_file_name);
        if (!(*force_ifgz_stream))
        {
            cerr << "Transport_Closure_Solver.cpp: CRITICAL ERROR: Cannot open forcing function file.\n";
            exit(1);
        }

        // Define the grid function with the forcing function data in the stream
        force_gf = new GridFunction(mesh, *force_ifgz_stream);
        
        // Scale the grid function by negative 1, since we assume it to be on the left-hand-side but move it to the righ-hand-side for the linear form
        *force_gf *= -1.0;

        cout << "Complete." << endl;
    }
    else
    {
        // Create the grid function for the forcing function
        force_gf = new GridFunction(fespace_f);
        
        // Initialize the forcing function as zero everywhere
        *force_gf = 0.0;
    }
    
    /*
    // ===============================================================
    //   Load the second forcing function
    // ===============================================================
    // Declare a grid function for the forcing function
    GridFunction *force2_gf;
    GridFunction *tester = new GridFunction(fespace_f);
        
    // Use a previously saved gridfunction for the forcing function, or define the field as zero everywhere
    if (forcing_function_file_name != "None" && forcing_function2_file_name != "None")
    {
        cout << "TCS:   Loading forcing function... ";

        // Create an ifgzstream pointer to read in the file containing the forcing function
        ifgzstream *force2_ifgz_stream = NULL;
        force2_ifgz_stream = new ifgzstream(*forcing_function2_dir + forcing_function2_file_name);
        if (!(*force2_ifgz_stream))
        {
            cerr << "Transport_Closure_Solver.cpp: CRITICAL ERROR: Cannot open the second forcing function file.\n";
            exit(1);
        }
        
        // Define the grid function with the forcing function data in the stream
        force2_gf = new GridFunction(mesh, *force2_ifgz_stream);
        
        // Multiply the grid functions together
        GridFunctionCoefficient forcing_function_dummy(force_gf);
        GridFunctionCoefficient forcing_function2(force2_gf);
        force_gf->Save((*closure_output_dir + "test1.gf").c_str());
        force2_gf->Save((*closure_output_dir + "test2.gf").c_str());
        ProductCoefficient gf_product(forcing_function_dummy, forcing_function2);
        
        tester->ProjectCoefficient(gf_product);
        force_gf->ProjectCoefficient(gf_product);
        tester->Save((*closure_output_dir + "polyO2_forcing_function.gf").c_str());
        
        cout << "Complete." << endl;
    }
    */
    
    // Create a GridFunctionCoefficient or VectorGridFunctionCoefficient from the gridfunction for the advection term
    GridFunctionCoefficient forcing_function(force_gf);


    // ===============================================================
    //   Define/Prepare boundary conditions.
    // ===============================================================
    cout << "TCS:   Defining boundary condition functions... ";
    
    // Here, we are defining two marker arrays to mark 1.) the inlet, and 2.) the outlet.
    Array<int> marker_inlet_BC(mesh->bdr_attributes.Max()); // Define a marker array for the inlet BC
    marker_inlet_BC = 0; // Initialize marker array to "don't apply any essential BC"
    if (useInlet == 1) { marker_inlet_BC[(int)(*mesh_info["scalar_closure"])["inlet1"] - 1] = 1; } // Set the first bdr attr group (i.e., physical group defined as 1 in gmsh file) as being the inlet
    int BC_toggle;
    if (active_inlet == 1 && forcing_function_file_name == "None") { BC_toggle = 1; } else { BC_toggle = 0; }
    FunctionCoefficient inlet_BC( [&](const Vector &x) -> double { return inlet_BC_func(x, BC_toggle); } ); // Define a function coefficient that will apply "inlet_BC_func", which varies in space
    
    // This can be modified to specify an Dirichlett condition at the outlet
    //Array<int> marker_outlet_BC(mesh->bdr_attributes.Max()); // Define a marker array for the no-slip BC.
    //marker_outlet_BC = 0; // Initialize marker array to "don't apply any essential BC".
    //marker_outlet_BC[mesh_info["scalar_closure"]["outlet1"] - 1] = 1; // Set the third bdr attr group (i.e., physical group defined as 3 in gmsh file) as being the no-slip condition.
    //FunctionCoefficient outlet_BC( [&](const Vector &x) -> double { return outlet_BC_func(x); } ); // Define a function coefficient that will apply "inlet_BC_func", which varies in space
    
    cout << "Complete." << endl;


    // ===============================================================
    //   Define block structure of the system.
    // ===============================================================
    cout << "TCS:   Defining the block system structure... ";

    // Notes: - This defines an array of offsets for each variable. The last component of the Array is the sum of the dimensions
    //          of each block.
    Array<int> block_offsets(3); // The argument should be the number of variables + 1
    block_offsets[0] = 0;
    block_offsets[1] = fespace_c->GetVSize();
    block_offsets[2] = N_AR;
    block_offsets.PartialSum();
    
    cout << "Complete." << endl;


    // ===============================================================
    //   Define the solution and right-hand-side (block) vectors and initialize them.
    // ===============================================================
    cout << "TCS:   Initializing the solution vector and RHS vector... ";

    BlockVector sol_BLK(block_offsets), b_BLK(block_offsets);
    sol_BLK = 0.0;
    b_BLK = 0.0;

    cout << "Complete." << endl;
    

    // ===============================================================
    //   Create the linear form for the forcing function and initialize b_BLK with it.
    // ===============================================================
    if (forcing_function_file_name != "None")
    {
        cout << "TCS:   Defining the RHS vector... ";

        LinearForm *varf_force = new LinearForm(fespace_f);
        varf_force->AddDomainIntegrator(new DomainLFIntegrator(forcing_function));
        varf_force->Assemble();
        b_BLK.GetBlock(0).SetVector(*varf_force, block_offsets[0]);

        cout << "Complete." << endl;
    }
    

    // ===============================================================
    //   Create the bilinear form for the diffusion term.
    // ===============================================================
    cout << "TCS:   Assembling bilinear forms and stiffness matrix... ";

    // Initiate the bilinear form of the diffusion term
    BilinearForm *varf_cdiff(new BilinearForm(fespace_c));
    ConstantCoefficient Diff_coef(omega);
    varf_cdiff->AddDomainIntegrator(new DiffusionIntegrator(Diff_coef));
    if (active_advection == 1) { varf_cdiff->AddDomainIntegrator(new ConvectionIntegrator(fluid_velocity)); }
    varf_cdiff->Assemble();
    
    // Define a reference grid function to the concentration solution. This is used to edit the matrix of the bilinear form so that the Dirichlett boundary conditions are applied 
    GridFunction BC_projection;
    BC_projection.MakeRef(fespace_c, sol_BLK.GetBlock(0), block_offsets[0]);
    
    // Apply the inlet BC to the diffusion bilinear form
    if (useInlet == 1)
    {
        //BC_projection.ProjectCoefficient(inlet_BC);
        BC_projection.ProjectBdrCoefficient(inlet_BC, marker_inlet_BC);
        varf_cdiff->EliminateEssentialBC(marker_inlet_BC, BC_projection, b_BLK.GetBlock(0)); // Sets the rows in "varf_viscous" corresponding to DOFs marked by "marker_inlet_BC" to 0.0 and assigns a diagonal as 1.0. Then, for the same rows, sets the components of b_BLK to those of u_BC_apply.
    }
    
    // Apply the outlet BC to the diffusion bilinear form
    //BC_projection.ProjectCoefficient(outlet_BC);
    //varf_cdiff->EliminateEssentialBC(marker_outlet_BC, BC_projection, b_BLK.GetBlock(0));
    
    // Finish the bilinear form by "finalizing"
    varf_cdiff->Finalize();

    cout << "Complete." << endl;


    // ===============================================================
    //   Create the bilinear form for the averaging operator/forcing unknowns.
    // ===============================================================
    cout << "TCS:   Assembling averaging operator... ";

    // Notes: - We only create the bilinear form matrix for the averaging operator, because its transpose is identically
    //          the bilinear form of the forcing unknowns in the mass conservation equation.
    //
    // Use the AveragingOperator class to obtain a matrix that can be multiplied by the solution vector to obtain the average solution in each AR
    AveragingOperator *avgOp = new AveragingOperator(fespace_c);
    SparseMatrix BM_cavg(avgOp->Getavg_mat());
    Array<double> AR_areas(avgOp->GetAR_areas());
    
    cout << "Complete." << endl;
    

    // ===============================================================
    //   Define the block operator for the system.
    // ===============================================================
    cout << "TCS:   Assembling the block system... ";

    // Initiate block operator
    BlockOperator closureOp_BLK(block_offsets);
    
    // Initiate pointer to the sparse matrix of the diffusion bilinear form
    SparseMatrix &BM_cdiff(varf_cdiff->SpMat());

    // Create the transpose of the averaging operator matrix (this will be the bilinear form/operator of the unknown forcing terms in the mass conservation equation)
    TransposeOperator *BM_cavg_T = NULL;
    BM_cavg_T = new TransposeOperator(&BM_cavg);
    
    // Set the blocks of the closure system block operator
    closureOp_BLK.SetBlock(0, 0, &BM_cdiff);
    closureOp_BLK.SetBlock(0, 1, BM_cavg_T);
    closureOp_BLK.SetBlock(1, 0, &BM_cavg);

    // Set the averaging operator to be 1 in the specified averaging region
    if (active_AR != -1 && forcing_function_file_name == "None")
    {
        b_BLK.GetBlock(1).Elem(active_AR) += AR_areas[active_AR];
    }

    cout << "Complete." << endl;
    

    // ===============================================================
    //   Construct the preconditioner operator.
    // ===============================================================
    cout << "TCS:   Assembling the preconditioners... ";

    // Notes: - Here, we use Symmetric Gauss-Seidel to approximate the inverse of the "a" Schur complement.
    //
    //                      P = [ BM_cdiff      0    ] = [ BM_cdiff                 0                ]
    //                          [    0      S_actual ] = [    0      BM_cavg BM_cdiff^(-1) BM_cavg_T ]
    //
    //    We use this ->      ~ [ BM_cdiff  0 ] = [ BM_cdiff                    0                   ]
    //                          [     0     S ]   [    0      BM_cavg diag(BM_cdiff)^(-1) BM_cavg_T ]
    //
    SaddlePointBlockPreconditioner closurePC(block_offsets);
    if (is3D) { closurePC.BuildPreconditioner(BM_cdiff, BM_cavg, "DSmoother", "GSSmoother"); }
    else { closurePC.BuildPreconditioner(BM_cdiff, BM_cavg, "BlockILU", "BlockILU"); }
    
    cout << "Complete." << endl;
    


    /*
    if (forcing_function_file_name != "None" && forcing_function2_file_name != "None")
    {
        tester->Save((*closure_output_dir + "Meep.gf").c_str());
        Vector tester2(10);
        BM_cavg.Mult(*tester, tester2);
        for (int i = 0; i < tester2.Size(); i++)
        {
            cout << tester2[i] << endl;
        }
        //tester2.Save((*closure_output_dir + "Applied_L_to_force.gf").c_str());
    }
    */
    /*
    if (forcing_function_file_name != "None" && forcing_function2_file_name != "None")
    {
        Array<int> block_offsets2(2); // The argument should be the number of variables + 1
        block_offsets2[0] = 0;
        block_offsets2[1] = fespace_c->GetVSize();
        block_offsets2.PartialSum();
        BlockOperator test_op(block_offsets2);
        test_op.SetBlock(0, 0, &BM_cdiff);
        BlockDiagonalPreconditioner test_PC(block_offsets2);
        invBM_cdiff = new GSSmoother(BM_cdiff);
        invBM_cdiff->iterative_mode = false;
        test_PC.SetDiagonalBlock(0, invBM_cdiff);

        GMRESSolver solver;
        solver.SetAbsTol(atol);
        solver.SetRelTol(rtol);
        solver.SetMaxIter(maxIter);
        solver.SetOperator(test_op);
        solver.SetPreconditioner(test_PC);
        solver.SetPrintLevel(1);
        solver.Mult(b_BLK.GetBlock(0), sol_BLK.GetBlock(0));
    }
    */
    
    // ===============================================================
    //   Solve the system with MINRES.
    // ===============================================================
    cout << "TCS:   Solving" << endl;

    //MINRESSolver solver;
    GMRESSolver solver;
    solver.SetAbsTol(atol);
    solver.SetRelTol(rtol);
    solver.SetMaxIter(maxIter);
    solver.SetOperator(closureOp_BLK);
    solver.SetPreconditioner(closurePC);
    solver.SetPrintLevel(1);
    solver.Mult(b_BLK, sol_BLK);

    if (solver.GetConverged())
    {
        cout << "GMRES converged in " << solver.GetNumIterations()
            << " iterations with a residual norm of "
            << solver.GetFinalNorm() << ".\n";
    }
    else
    {
        cout << "GMRES did not converge in " << solver.GetNumIterations()
            << " iterations. Residual norm is " << solver.GetFinalNorm()
            << ".\n";
    }


    if (forcing_function_file_name != "None" && forcing_function2_file_name != "None")
    {
        Vector tester2(10);
        BM_cavg.Mult(sol_BLK.GetBlock(0), tester2);
        for (int i = 0; i < tester2.Size(); i++)
        {
            cout << tester2[i] << endl;
        }
        //tester2.Save((*closure_output_dir + "Applied_L_to_force.gf").c_str());
    }

    

    // ===============================================================
    //   Extract and save the solutions and mesh.
    // ===============================================================
    // Extract and save the closure variable solution
    GridFunction c_sol;
    c_sol.MakeRef(fespace_c, sol_BLK.GetBlock(0), 0);
    c_sol.Save((*closure_output_dir + closure_output_file_name).c_str());
    
    // Extract the residual data from the block solution vector (and print the data to cout for viewing)
    vector<double> a_sol_cpp_vec;
    cout << "Closure residuals:" << endl;
    for (int i_a = 0; i_a < N_AR; i_a++)
    {
        a_sol_cpp_vec.push_back( -sol_BLK.GetBlock(1).Elem(i_a) ); // The negative is because we solve with "a" on the left, but expect it on the right for the model
        cout << a_sol_cpp_vec[i_a] << endl;
    }
    

    // Load the previous residual map (if there is one) to append/edit previous residual data
    JSONDict residual_dict, time_func_dict, iter_dict, res_vec_dict, closure_file_name_dict;
    
    // Define the key strings based on the simulation
    string sim_key;
    string AR_number;
    if (active_inlet == 1) { sim_key = "inlet"; AR_number = "0"; }
    else { sim_key = "avg_c"; AR_number = to_string(active_AR); }

     // If the residual file exists, load it into the dictionaries
    ifstream file(residual_output_file_path);
    if (file)
    {
        residual_dict.loadFromFile(residual_output_file_path);
        residual_dict.getValue(sim_key, time_func_dict);
        time_func_dict.getValue("iter_" + to_string(recursive_iter), iter_dict);
        iter_dict.getValue("residuals", res_vec_dict);
        iter_dict.getValue("closure file name", closure_file_name_dict);
    }
    
    // Create the new dictionary
    res_vec_dict[AR_number] = a_sol_cpp_vec;
    closure_file_name_dict[AR_number] = closure_output_file_name;
    
    // Save the new dictionary into the loaded dictionaries
    iter_dict["residuals"] = &res_vec_dict;
    iter_dict["closure file name"] = &closure_file_name_dict;
    time_func_dict["iter_" + to_string(recursive_iter)] = &iter_dict;
    residual_dict[sim_key] = &time_func_dict;
    residual_dict.saveToFile(residual_output_file_path);

    // Save the pore areas to the mesh_info.txt file (i.e., mesh_info)
    JSONDict AR_dict_temp = *mesh_info["AR"];
    vector<double> AR_areas_temp(AR_areas.GetData(), AR_areas.GetData() + AR_areas.Size());
    AR_dict_temp["pore_areas"] = AR_areas_temp;
    mesh_info["AR"] = &AR_dict_temp;
    mesh_info.saveToFile(mesh_info_file_path);
    
    // Save the mesh
    mesh->Save(*mesh_output_dir + mesh_output_file_name);


    // ===============================================================
    // 17. Free the used memory.
    // ===============================================================
    delete BM_cavg_T;
    delete avgOp;
    delete varf_cdiff;
    // delete fespace_c; // Deleted with avgOp.
    delete fec_c;
    if (translations.size() == 0) { delete mesh; }
    
    return 0;
}





// Define the inlet velocity condition
double inlet_BC_func(const Vector &p, int inlet_value_toggle)
{
   int dim = p.Size();

   real_t x = p(0);
   real_t y = p(1);
   real_t z = (dim == 3) ? p(2) : 0.0;
   
   if (inlet_value_toggle == 1) { return 1.0; }
   else { return 0.0; }
}


// Define the no-slip velocity condition
void outlet_BC_func(const Vector &p, Vector &F)
{
   int dim = p.Size();

   real_t x = p(0);
   real_t y = p(1);
   real_t z = (dim == 3) ? p(2) : 0.0;
   
   F(0) = 0.0;
}

