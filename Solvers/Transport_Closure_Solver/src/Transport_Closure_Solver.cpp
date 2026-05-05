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
//

#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <vector>
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
    constexpr GeneralizedResidualManager() {};
    GeneralizedResidualManager(double alpha_, double beta_, double gamma_);
    void SetParams(double alpha_, double beta_, double gamma_);
    constexpr void ConstructParamVecs(int N_AR_global, const vector<int> &AR_inds_loc2glob) {};
    constexpr void ConstructParamVecs(int N_AR_global, const vector<int> &AR_inds_loc2glob, bool test) {};
    constexpr void ImplementAlpha(BilinearForm *&varf) {};
    constexpr void ImplementBeta(SparseMatrix &cavg) {};
    constexpr void ImplementGamma(SparseMatrix &BM_avgavg, const vector<int> &AR_inds_loc2glob) {};
    void AlterGammaVec(const vector<int> &AR_tags, const vector<vector<int>> &AR_neighbors, int active_AR_global, int procedureID = -1);
    static constexpr vector<double>* GetParamVec(const char* param_name) { return nullptr; };
    bool LoadParamDicts(JSONDict &eq_dict);
    void ImplementGeneralizedResidual(vector<double> &temp, vector<string> keys, int N_AR, bool isKMat = false);
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




class NormalDotVelocityCoefficient : public Coefficient
{
private:
    VectorCoefficient &u_;

public:
    NormalDotVelocityCoefficient(VectorCoefficient &u) : u_(u) { }

    double Eval(ElementTransformation &T,
                const IntegrationPoint &ip) override
    {
        // Evaluate the vector
        Vector uval_;
        u_.Eval(uval_, T, ip);
        
        // Compute physical normal from face Jacobian
        Vector n_(u_.GetVDim());
        mfem::CalcOrtho(T.Jacobian(), n_);
        
        // Normalize for unit normal
        const double nn = n_.Norml2();
        if (nn > 0.0) { n_ /= nn; }

        // Return vector dot normal
        return uval_ * n_;
    }
};

class NegativeVectorCoefficient : public VectorCoefficient
{
private:
    VectorGridFunctionCoefficient &u_;

public:
    // Inherit the constructors
    using VectorCoefficient::VectorCoefficient;
    
    // Define the class constructors
    NegativeVectorCoefficient(VectorGridFunctionCoefficient &u) : VectorCoefficient(u.GetVDim()), u_(u) { }

    void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
    {
        // Evaluate the vector
        u_.Eval(V, T, ip);
        // Make it negative
        V *= -1.0;
    }
};




// Define functions that are used in "main", but are defined after "main".
double inlet_BC_func(const Vector &p, int inlet_value_toggle);
//void outlet_BC_func(const Vector &p, Vector &F);

// Define a global "parameters struct" so that the BC functions can take in parameters
struct Params
{
    string FILENAME = "Transport_Closure_Solver.cpp";
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
    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";

    vector<string> reaction_pg_names; // Reaction physical group names (defined with gmsh)

    int is3D = 0;

    // Multi-threaded/Parallel controls
    string solve_mode = "serial";
    int save_mesh = 1;
    
    int order = 2;
    int active_advection = 1;
    vector<int> isPeriodic = {0, 0, 0};
    int useInlet = 0;
    int useReactions = 0;
    int importFluidMesh = 0;
    
    int active_AR = -1;
    int active_inlet = 0;
    vector<int> active_reactions;
    int active_reaction = -2;
    int active_adjoint = 0;

    // Physical equation parameters
    double Pe_s = 1.0;
    double omega = 1.0;
    vector<double> Da_s;

    // Linear solver controls 
    int maxIter = 100000;
    double rtol = 1.0e-10;
    double atol = 0.0;
    
    int recursive_iter = 0;

    // Variables for changing the residual form to localize the solution
    double resAvg_alpha = 0.0;
    double resAvg_beta = 1.0;
    double resAvg_gamma = 0.0;

    // Variables for reducing the mesh size to around the localized closure solution
    int N_neighbor_layers;
    int useLocalMesh = 0;
    int saveLocalMesh = 0;
    
    

    string output_dir = "./"; // For general output; this is adjusted by the config file for the different variations

    string fluid_velocity_dir = "./";
    string fluid_velocity_file_name = "u_sol.gf";
    string fluid_velocity_mesh_file_name = "mesh.mesh";

    string residual_output_dir = "./";
    string residual_output_file_name_prefix = "a_sol";
    string residual_output_file_name_suffix = ".txt";

    string closure_output_dir = "./";
    string closure_output_file_name_prefix = "chi";
    string closure_output_file_name_suffix = ".gf";
    string closure_output_file_name_ID = "";
    
    //string log_file_file_name_prefix = "log";
    //string log_file_file_name_suffix = ".txt";
    
    string mesh_output_dir = "./";
    string submesh_output_dir_extn = "submeshes/";
    //string mesh_output_file_name = "mesh.mesh";
    string mesh_output_file_name_prefix = "mesh";
    string mesh_output_file_name_suffix = ".mesh";
    
    string forcing_function_dir = "./";
    string forcing_function_file_name = "None";
    
    string forcing_function2_dir = "./";
    string forcing_function2_file_name = "None";


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

        sub_dict = *mesh_dict["details"];
        sub_dict.getValue("is3D", is3D);

        //mesh_dict.getValue("cut physical groups", sub_dict);
        sub_dict = *mesh_dict["cut physical groups"];
        sub_dict.getValue("physical group cut names", reaction_pg_names);

        
        JSONDict stokes_dict = *configData["stokes"];
        
        sub_dict = *stokes_dict["output path"];
        sub_dict.getValue("directory", fluid_velocity_dir);
        sub_dict.getValue("velocity file name", fluid_velocity_file_name);
        sub_dict.getValue("mesh file name", fluid_velocity_mesh_file_name);


        JSONDict closure_dict = *configData["scalar closure"];

        sub_dict = *closure_dict["simulation parameters"];
        sub_dict.getValue("order", order);
        sub_dict.getValue("active advection", active_advection);
        sub_dict.getValue("isPeriodic", isPeriodic);
        sub_dict.getValue("use inlet", useInlet);
        sub_dict.getValue("use reactions", useReactions);
        sub_dict.getValue("import fluid mesh", importFluidMesh);
        sub_dict.getValue("active adjoint", active_adjoint);
        
        sub_dict = *closure_dict["residual parameters"];
        sub_dict.getValue("resAvg_alpha", resAvg_alpha);
        sub_dict.getValue("resAvg_beta", resAvg_beta);
        sub_dict.getValue("resAvg_gamma", resAvg_gamma);

        sub_dict = *closure_dict["mesh localization parameters"];
        sub_dict.getValue("use local mesh", useLocalMesh);
        sub_dict.getValue("N_neighbor_layers", N_neighbor_layers);
        sub_dict.getValue("save local mesh", saveLocalMesh);
        
        sub_dict = *closure_dict["closure parameters"];
        sub_dict.getValue("active averaging region", active_AR);
        sub_dict.getValue("active inlet", active_inlet);
        sub_dict.getValue("active reactions", active_reactions);
        
        sub_dict = *closure_dict["physics parameters"];
        sub_dict.getValue("Pe_s", Pe_s);
        sub_dict.getValue("omega", omega);
        sub_dict.getValue("Da_s", Da_s);

        sub_dict = *closure_dict["solver parameters"];
        sub_dict.getValue("max iterations", maxIter);
        sub_dict.getValue("rel tol", rtol);
        sub_dict.getValue("abs tol", atol);

        sub_dict = *closure_dict["residual path"];
        sub_dict.getValue("directory", residual_output_dir);
        sub_dict.getValue("file name prefix", residual_output_file_name_prefix);
        sub_dict.getValue("file name suffix", residual_output_file_name_suffix);

        sub_dict = *closure_dict["closure path"];
        sub_dict.getValue("directory", closure_output_dir);
        sub_dict.getValue("file name prefix", closure_output_file_name_prefix);
        sub_dict.getValue("file name suffix", closure_output_file_name_suffix);
        
        sub_dict = *closure_dict["mesh output path"];
        sub_dict.getValue("directory", mesh_output_dir);
        //sub_dict.getValue("file name", mesh_output_file_name);
        sub_dict.getValue("file name prefix", mesh_output_file_name_prefix);
        sub_dict.getValue("file name suffix", mesh_output_file_name_suffix);

        sub_dict = *closure_dict["forcing function path"];
        sub_dict.getValue("directory", forcing_function_dir);
        sub_dict.getValue("file name", forcing_function_file_name);
        
        //sub_dict = *closure_dict["forcing function 2 path"];
        //sub_dict.getValue("directory", forcing_function2_dir);
        //sub_dict.getValue("file name", forcing_function2_file_name);
    }

    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_file_path = mesh_dir + mesh_file_name;
    string mesh_info_file_path = mesh_info_dir + mesh_info_file_name;
    string fluid_velocity_file_path = fluid_velocity_dir + fluid_velocity_file_name;
    string fluid_velocity_mesh_file_path = fluid_velocity_dir + fluid_velocity_mesh_file_name;
    string residual_output_file_name = residual_output_file_name_prefix + residual_output_file_name_suffix;
    //string closure_output_file_name = closure_output_file_name_prefix + closure_output_file_name_suffix;

    //string log_file_file_name = log_file_file_name_prefix + log_file_file_name_suffix;
    //string log_file_file_path = residual_output_dir + log_file_file_name;


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
    args.AddOption(&active_reaction, "-R", "--active_reaction", "Toggle which heterogeneous reaction is active. Provide 0, 1, ... to say which is active.");
    args.AddOption(&fluid_velocity_file_path, "-v", "--fluid_velocity_file_path", "Fluid velocity solution file path (best to define in the config file).");
    args.AddOption(&fluid_velocity_mesh_file_path, "-V", "--fluid_velocity_mesh_file_path", "Fluid velocity mesh file path (best to define in the config file).");
    args.AddOption(&active_advection, "-A", "--active_advection", "Toggle whether to consider advection or not.");
    args.AddOption(&Pe_s, "-P", "--Peclet_Number", "The O(1) Peclet Number to use (i.e., this is Pe*\epsilon).");
    //args.AddOption(&closure_output_file_name, "-N", "--closure_output_file_name", "Declare the file name of the output file for the closure variable solution (Note: the file name should end in '.gf').");
    args.AddOption(&closure_output_file_name_ID, "-N", "--closure_output_file_name_ID", "Declare the file name ID of the output file for the closure variable solution.");
    args.AddOption(&residual_output_file_name, "-r", "--residual_output_file_name", "Declare the file name of the output file containing the closure residuals (Note: the file name should end in '.txt').");
    //args.AddOption(&log_file_file_name, "-L", "--log_file_file_name", "Declare the file name of the log file (Note: the file name should end in '.txt').");
    args.AddOption(&forcing_function_file_name, "-F", "--forcing_function_file_name", "Declare the file name of the forcing function (Note: the file name should end in '.gf').");
    //args.AddOption(&forcing_function2_file_name, "-f", "--forcing_function2_file_name", "Declare the file name of the second forcing function (Note: the file name should end in '.gf').");
    args.AddOption(&recursive_iter, "-I", "--recursive_iter", "Provide the int that describes the iteration of this run.");
    args.AddOption(&save_mesh, "-s", "--save_mesh", "Toggle whether to save the mesh or not. 0 = No, 1 = Yes.");
    args.AddOption(&solve_mode, "-S", "--solve_mode", "Define how the code is being run. Options: 'serial' or 'ensemble' (i.e., closure problems are being solved as an ensamble in parallel)");
    args.AddOption(&active_adjoint, "-j", "--active_adjoint", "Toggle whether the adjoint closure problem will be solved. 0 = No, 1 = Yes.");
    args.ParseCheck();

    
    // ===============================================================
    //   Define variables based on the options provided in the parser
    // ===============================================================
    // When using the parallel solver, if advection is active/the fluid velocity was previous solved and used here, the mesh saved with the fluid velocity must be used. Also, use consistent BCs (i.e., both fluid and transport simulations should be non-periodic)
    if (active_advection == 1 && importFluidMesh == 1) { mesh_file_path = fluid_velocity_mesh_file_path; }

    // Define the closure output file name with the prefix, ID, and suffix
    string closure_output_file_name = closure_output_file_name_prefix + closure_output_file_name_ID + closure_output_file_name_suffix;
    // Define the mesh output file name depending on whether the local mesh is being used/saved
    string mesh_output_file_name;
    if (useLocalMesh == 1 && saveLocalMesh == 1) {
        if (closure_output_file_name_ID == "") { cout << globalVars.FILENAME << ": WARNING: No closure file ID has been provided, but solver is saving the local mesh. If multiple closure problems are being solved, the local mesh file will be overwritten, and results will not be able to be viewed." << endl; }
        mesh_output_file_name = mesh_output_file_name_prefix + closure_output_file_name_ID + mesh_output_file_name_suffix; }
    else { mesh_output_file_name = mesh_output_file_name_prefix + mesh_output_file_name_suffix; }
    // Define the residual output file path with the corresponding directory and file name
    string residual_output_file_path = residual_output_dir + residual_output_file_name;
    
    // Determine which of the reactions is active in the curret closure problem
    if (active_reaction != -2) {
        active_reactions.clear();
        for (int i = 0; i < active_reaction + 1; i++) {
            if (i == active_reaction) { active_reactions.push_back( 1 ); }
            else { active_reactions.push_back( 0 ); } } }
    
    // If the code for implementing the general residual form is available, and any of the corresponding parameters
    // are different than what would collapse to the canonical residual form, implement the general residual form
    bool useGenResForm = false;
    GeneralizedResidualManager gen_res_manager;
    if (genResFormAvailable && (resAvg_alpha != 0.0 || resAvg_beta != 1.0 || resAvg_gamma != 0.0)) {
        useGenResForm = true;
        gen_res_manager.SetParams(resAvg_alpha, resAvg_beta, resAvg_gamma);
    }
    
    // Define whether the simulation is periodic in at least one direction
    bool meshIsPeriodic = false;
    for (int i_prdc = 0; i_prdc < isPeriodic.size(); i_prdc++) { if (isPeriodic[i_prdc] == 1) { meshIsPeriodic = true; break; } }
    
    
    // ===============================================================
    //   Create a log file (if ensemble solving, output to cout because each thread is obtaining the cout of this code)
    // ===============================================================
    //ofstream LOG(log_file_file_path);
    //streambuf* cout_backup = nullptr;
    //if (solve_mode != "serial")
    //{
    //    // Back up the original streambuf of cout
    //    cout_backup = cout.rdbuf();
    //
    //    // Redirect cout's streambuf to the file's streambuf
    //    cout.rdbuf(LOG.rdbuf());
    //}


    // ===============================================================
    //   Get the mesh, mesh information, and corresponding variables.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Obtaining mesh and mesh info... ";

    // Get "global" mesh info from the mesh info file (this is information about the full/parent/global mesh)
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_file_path);
    globalVars.L = (*mesh_info["geometry"])["L"];

    JSONDict AR_dict = *mesh_info["AR"];
    int N_AR_global = (int)(AR_dict["total_number"]); // Get the number of averaging regions defined in the mesh file
    vector<int> AR_tags = AR_dict["tags"];
    vector<vector<int>> AR_neighbors = AR_dict["neighbors"];
    vector<vector<int>> AR_neighbor_macroIDs; AR_dict.getValue("neighbor macroIDs", AR_neighbor_macroIDs);
    vector<vector<int>> macroIDCoordsKey; AR_dict.getValue("macroID-coords key", macroIDCoordsKey);
    

    // Define additional variables using the global mesh information
    int N_AR = N_AR_global;
    int active_AR_global = active_AR;
    vector<int> AR_inds_loc2glob(N_AR_global); for (int i = 0; i < N_AR_global; i++) { AR_inds_loc2glob[i] = i; }
    vector<int> AR_macroIDs; // Would need an additional method for definition while using a non-local mesh

    // Use the provided mesh file to define the mesh 
    MeshManager mesh_manager(mesh_file_path);
    
    // Make the mesh periodic if required by isPeriodic
    mesh_manager.MakePeriodicMesh(isPeriodic, globalVars.L, true);
    
    // Make a local mesh around active_AR/the inlet, if called for by useLocalMesh
    if (useLocalMesh == 1) {
        // Make the local mesh
        if (active_inlet == 1) {
            vector<int> inlet_AR_neighbors = (*mesh_info["scalar_closure"])["inlet1 AR neighbors"]; for (int i = 0; i < inlet_AR_neighbors.size(); i++) { inlet_AR_neighbors[i] -= 1; }
            if (meshIsPeriodic) { mesh_manager.MakeLocalMesh(inlet_AR_neighbors, N_neighbor_layers, AR_tags, AR_neighbors, AR_neighbor_macroIDs, macroIDCoordsKey, isPeriodic); }
            else { mesh_manager.MakeLocalMesh(inlet_AR_neighbors, N_neighbor_layers, AR_tags, AR_neighbors); }
        }
        else {
            if (meshIsPeriodic) { mesh_manager.MakeLocalMesh({active_AR}, N_neighbor_layers, AR_tags, AR_neighbors, AR_neighbor_macroIDs, macroIDCoordsKey, isPeriodic); }
            else { mesh_manager.MakeLocalMesh({active_AR}, N_neighbor_layers, AR_tags, AR_neighbors); }
            active_AR = mesh_manager.get_central_AR_inds_local()[0];
        }
        // Reassign variable "mesh" to the local mesh
        mesh_manager.UseLocalMesh();
        // If the solve mode is not serial, and the save mesh is not 1, delete the parent mesh; it is no longer needed
        if (solve_mode != "serial" && save_mesh != 1) { mesh_manager.DeleteParentMesh(); }
        // Update variables based on the local mesh used
        N_AR = mesh_manager.get_N_AR_local();
        AR_tags = mesh_manager.get_AR_kept();
        AR_inds_loc2glob = mesh_manager.get_AR_inds_loc2glob();
        AR_macroIDs = mesh_manager.get_AR_macroIDs();
    }

    // Get a pointer to the mesh being used
    Mesh *mesh = mesh_manager.GetMesh();

    cout << "Complete." << endl;
    cout << globalVars.FILENAME << ":   Total number of averaging regions: " << N_AR_global << endl;
    cout << globalVars.FILENAME << ":   Number of averaging regions considered in closure problem: " << N_AR << endl;



    // There is currently no adjoint problem for an inlet force (i.e., Dirichlet BC = 1). Therefore, just save the
    // AR_inds and macro IDs in the residuals file, as they are used in the conversion of adjoint residuals to
    // generalized residuals
    if (active_adjoint == 1 && active_inlet == 1) {
        // Define the sim_key and AR_number
        string sim_key = "inlet"; string AR_number = "0";
        
        // Create the saver for the residuals and save them to the text file
        ResultsSaver Saver;
        Saver.SetSolveMode(solve_mode);
        Saver.SetVariables(sim_key, recursive_iter, {"transport eq"}, {closure_output_file_name}, AR_number);
        Saver.InitializeJSONDicts_EqLevel({1});
        
        Saver.StoreResidualParameters("AR_inds", {{AR_inds_loc2glob}});
        if (AR_macroIDs.size() > 0) { Saver.StoreResidualParameters("AR_macroIDs", {{AR_macroIDs}}); }
        
        Saver.SaveResidualDictionary(residual_output_file_path);

        return 0;
    }



    // ===============================================================
    //   Define finite element spaces for the concentration.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Defining finite element spaces... ";

    // Finite element space for concentration
    FiniteElementCollection *fec_c = new H1_FECollection(order, mesh->Dimension());
    FiniteElementSpace *fespace_c = new FiniteElementSpace(mesh, fec_c, 1);
    
    cout << "Complete." << endl;
    
    // Print the number of unknowns for concentration and total
    cout << globalVars.FILENAME << ":   Number of unknowns in concentration: " << fespace_c->GetTrueVSize() << endl;
    cout << globalVars.FILENAME << ":   Number of unknowns in total: " << fespace_c->GetTrueVSize() << endl;

    
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
        cout << globalVars.FILENAME << ":   Loading fluid velocity... ";
        fluid_velocity_vgfc_manager = std::make_unique<VectorGridFunctionCoefficientManager>(fluid_velocity_file_path, fluid_velocity_mesh_file_path, Pe_s);
        cout << "Complete." << endl;

        // Create the local fluid velocity grid function if necessary
        if (useLocalMesh == 1) {
            cout << globalVars.FILENAME << ":   Obtaining the local fluid velocity... ";
            fluid_velocity_vgfc_manager->MakeLocalGridFunction(mesh_manager.GetLocalMesh()); // Local mesh (i.e., the submesh) is needed to make the TransferMap
            fluid_velocity_vgfc_manager->UseLocalGridFunction();
            cout << "Complete." << endl;
        }
        
        // Create the fluid velocity VectorGridFunctionCoefficient
        cout << globalVars.FILENAME << ":   Creating fluid velocity VectorGridFunctionCoefficient... ";
        fluid_velocity_vgfc_manager->MakeVectorGridFunctionCoefficient();
        fluid_velocity = fluid_velocity_vgfc_manager->GetVectorGridFunctionCoefficient();
        cout << "Complete." << endl;


        
        // ======== For Debugging ========
        /*
        //Mesh *mesh2 = fluid_velocity_vgfc_manager->GetGridFunctionMesh();
        FiniteElementCollection *fec_u_save = new H1_FECollection(2, mesh->Dimension());
        FiniteElementSpace *fespace_u_save = new FiniteElementSpace(mesh, fec_u_save, mesh->Dimension());
        GridFunction u_gf_save(fespace_u_save);

        //u_gf_save.ProjectGridFunction(*fluid_velocity_vgfc_manager->GetGridFunction());
        u_gf_save.ProjectCoefficient(*fluid_velocity);

        u_gf_save.Save((mesh_output_dir + "fluid_velocity_verification_plot.gf").c_str());
        
        ostringstream mesh_name; mesh_name << "fluid_velocity_verification_mesh.mesh";
        ofstream mesh_ofs((mesh_output_dir + mesh_name.str()).c_str());
        mesh_ofs.precision(8);
        mesh->Print(mesh_ofs);
        exit(1);
        */
        // ===============================
    }

    
    // ===============================================================
    //   Load the forcing function (i.e., the closure solution from the previous closure problem iteration)
    // ===============================================================
    // Declare pointers for loading/preparing the forcing function grid function coefficient
    std::unique_ptr<GridFunctionCoefficientManager> forcing_function_gfc_manager;
    GridFunctionCoefficient *forcing_function = nullptr;
    
    // If a forcing function is provided...
    if (forcing_function_file_name != "None")
    {
        // Create the forcing function grid function coefficient manager (here, we use the current mesh---parent or sub---, as it is likely that the saved forcing function of the same mesh. If not, a procedure similar to the velocity will have to be done)
        cout << globalVars.FILENAME << ":   Loading forcing function... ";
        double neg_one = -1.0;
        string forcing_function_path_name = forcing_function_dir + forcing_function_file_name;
        forcing_function_gfc_manager = std::make_unique<GridFunctionCoefficientManager>(forcing_function_path_name, mesh_manager.GetMesh(), neg_one);
        cout << "Complete." << endl;
        
        // Create the forcing function GridFunctionCoefficient
        cout << globalVars.FILENAME << ":   Creating forcing function GridFunctionCoefficient... ";
        forcing_function_gfc_manager->MakeGridFunctionCoefficient();
        forcing_function = forcing_function_gfc_manager->GetGridFunctionCoefficient();
        cout << "Complete." << endl;
        //forcing_function_gfc_manager->GetGridFunction()->Save("TESTING_GF.gf");



        // ======== For Debugging ========
        //FiniteElementCollection *fec_u_save = new H1_FECollection(2, mesh->Dimension());
        //FiniteElementSpace *fespace_u_save = new FiniteElementSpace(mesh, fec_u_save, 1);
        //GridFunction u_gf_save(fespace_u_save);

        //u_gf_save.ProjectGridFunction(*forcing_function_gfc_manager->GetGridFunction());
        //u_gf_save.ProjectCoefficient(*forcing_function);

        //u_gf_save.Save((mesh_output_dir + "fluid_velocity_verification_plot.gf").c_str());
        //mesh->Save(mesh_output_dir + "closure_mesh.mesh");
        //exit(1);

        //ostringstream mesh_name; mesh_name << "closure_mesh.mesh";
        //ofstream mesh_ofs((mesh_output_dir + mesh_name.str()).c_str());
        //mesh_ofs.precision(8);
        //mesh->Print(mesh_ofs);
        //exit(1);
        // ===============================
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
    
    // Create a GridFunctionCoefficient
    //GridFunctionCoefficient forcing_function(force_gf);


    // ===============================================================
    //   Define/Prepare boundary conditions.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Defining boundary condition functions... ";
    
    // Here, we are defining two marker arrays to mark 1.) the inlet, and 2.) the outlet.
    Array<int> marker_inlet_BC; //(mesh->bdr_attributes.Max()); // Define a marker array for the inlet BC
    marker_inlet_BC.SetSize(mesh->bdr_attributes.Max()); marker_inlet_BC = 0; // Initialize marker array to "don't apply any essential BC"
    vector<int> inlet_PGs = (*mesh_info["scalar_closure"])["inlet1"];
    if (useInlet == 1) {
        for (int i_BCin = 0; i_BCin < inlet_PGs.size(); i_BCin++) { marker_inlet_BC[inlet_PGs[i_BCin] - 1] = 1; } // Set the bdr attr group for the inlet to 1
        //marker_inlet_BC[(int)(*mesh_info["scalar_closure"])["inlet1"] - 1] = 1;
    }
    int BC_toggle; if (active_inlet == 1 && forcing_function_file_name == "None") { BC_toggle = 1; } else { BC_toggle = 0; }
    FunctionCoefficient inlet_BC( [&](const Vector &x) -> double { return inlet_BC_func(x, BC_toggle); } ); // Define a function coefficient that will apply "inlet_BC_func", which varies in space
    
    
    // Create the marker that labels the outlet boundary
    Array<int> marker_outlet_BC(mesh->bdr_attributes.Max()); marker_outlet_BC = 0; // Initialize a marker array for the outlet BC
    if (active_adjoint == 1) {
        vector<int> outlet_PGs = (*mesh_info["scalar_closure"])["outlet1"];
        for (int i_BCout = 0; i_BCout < outlet_PGs.size(); i_BCout++) { marker_outlet_BC[outlet_PGs[i_BCout] - 1] = 1; } // Set the bdr attr group for the outlet to 1
    }
    
    cout << "Complete." << endl;


    // ===============================================================
    //   Define block structure of the system.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Defining the block system structure... ";

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
    cout << globalVars.FILENAME << ":   Initializing the solution vector and RHS vector... ";

    BlockVector sol_BLK(block_offsets), b_BLK(block_offsets);
    sol_BLK = 0.0;
    b_BLK = 0.0;

    cout << "Complete." << endl;
    

    // ===============================================================
    //   Create the linear form for the forcing function and initialize b_BLK with it.
    // ===============================================================
    LinearForm *varf_force = nullptr, *varf_reaction = nullptr;
    if (forcing_function_file_name != "None")
    {
        cout << globalVars.FILENAME << ":   Defining the RHS vector (body force function)... ";

        varf_force = new LinearForm(forcing_function_gfc_manager->GetGridFunctionFES());
        varf_force->AddDomainIntegrator(new DomainLFIntegrator(*forcing_function));
        varf_force->Assemble();
        b_BLK.GetBlock(0).SetVector(*varf_force, block_offsets[0]);
        
        cout << "Complete." << endl;
    }

    // TODO: check compatibility with local mesh formulation
    if (useReactions == 1)
    {
        cout << globalVars.FILENAME << ":   Defining the RHS vector (Neumann BC), if there are any... ";
        
        for (int i_r = 0; i_r < active_reactions.size(); i_r++)
        {
            if (active_reactions[i_r] == 1)
            {
                // Check that there is a physical groups name and Da_s defined for the active reaction
                assert (reaction_pg_names.size() - 1 >= i_r);
                assert (Da_s.size() - 1 >= i_r);
                
                // Define a marker array for the heterogeneous reaction (i.e., Neumann BC)
                cout << globalVars.FILENAME << ": CRITICAL ERROR: Need to reconsider the tags for reactive boundaries." << endl;
                exit(1);
                Array<int> marker_reaction_BC(mesh->bdr_attributes.Max()); marker_reaction_BC = 0; // Initialize marker array
                marker_reaction_BC[(int)(*mesh_info["scalar_closure"])[reaction_pg_names[i_r]] - 1] = 1; // Assign marker array 1 for active reaction
                
                // Create, assemble, and add the linear form to the "b" vector
                LinearForm *varf_reaction = new LinearForm(fespace_c);
                ConstantCoefficient Da_s_coeff(Da_s[i_r]);
                varf_reaction->AddBoundaryIntegrator(new BoundaryLFIntegrator(Da_s_coeff), marker_reaction_BC);
                varf_reaction->Assemble();
                b_BLK.GetBlock(0) += *varf_reaction;
                delete varf_reaction;
            }
        }
        
        cout << "Complete." << endl;
    }


    // ===============================================================
    //   Create the bilinear form for the diffusion term.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Assembling bilinear forms and stiffness matrix... ";

    // Initiate the bilinear form of the diffusion term
    BilinearForm *varf_cdiff(new BilinearForm(fespace_c));
    ConstantCoefficient Diff_coef(omega);
    varf_cdiff->AddDomainIntegrator(new DiffusionIntegrator(Diff_coef));
    
    // If advection is active, add the advection term
    std::unique_ptr<NegativeVectorCoefficient> neg_fluid_velocity;
    std::unique_ptr<NormalDotVelocityCoefficient> n_dot_fluid_velocity;
    if (active_advection == 1) {
        if (active_adjoint == 1) {
            // If adjoint is active, make the velocity negative and add the advection term
            neg_fluid_velocity = std::make_unique<NegativeVectorCoefficient>(*fluid_velocity);
            varf_cdiff->AddDomainIntegrator(new ConvectionIntegrator(*neg_fluid_velocity));

            // Add the boundary flux term
            n_dot_fluid_velocity = std::make_unique<NormalDotVelocityCoefficient>(*fluid_velocity);
            varf_cdiff->AddBoundaryIntegrator(new BoundaryMassIntegrator(*n_dot_fluid_velocity), marker_outlet_BC);
        }
        else {
            varf_cdiff->AddDomainIntegrator(new ConvectionIntegrator(*fluid_velocity));
        }
    }

    // Implement alpha for generalized closure residual (if available)
    if (useGenResForm) { gen_res_manager.ImplementAlpha(varf_cdiff); }
    
    varf_cdiff->Assemble();
    

    // Apply the inlet BC to the diffusion bilinear form
    if (useInlet == 1) {
        // Define a reference grid function to the concentration solution. This is used to edit the matrix of the bilinear form so that the Dirichlett boundary conditions are applied
        GridFunction BC_projection;
        BC_projection.MakeRef(fespace_c, sol_BLK.GetBlock(0), block_offsets[0]);
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
    cout << globalVars.FILENAME << ":   Assembling averaging operator... ";

    // Notes: - We only create the bilinear form matrix for the averaging operator, because its transpose is identically
    //          the bilinear form of the forcing unknowns in the mass conservation equation.
    //
    // Use the AveragingOperator class to obtain a matrix that can be multiplied by the solution vector to obtain the average solution in each AR
    AveragingOperator *avgOp = new AveragingOperator(fespace_c, AR_tags);
    //SparseMatrix BM_cavg(avgOp->Getavg_mat());
    SparseMatrix *BM_cavg(avgOp->Getavg_mat());
    Array<double> AR_areas(avgOp->GetAR_areas());
    
    // If adjoint is active, we normalize BM_cavg by the area, and multiply by the area of the active_AR
    if (active_adjoint == 1) {
        SparseMatrix AR_areas_mat(AR_areas.Size(), AR_areas.Size());
        for (int i = 0; i < AR_areas.Size(); i++) { AR_areas_mat.Set(i, i, AR_areas[active_AR] / AR_areas[i]); }
        AR_areas_mat.Finalize();
        mfem::Mult(AR_areas_mat, *BM_cavg, BM_cavg);
        BM_cavg->Finalize();
    }

    cout << "Complete." << endl;
    

    // ===============================================================
    //   Define the block operator for the system.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Assembling the block system... ";
    
    // Initiate block operator
    BlockOperator closureOp_BLK(block_offsets);
    
    // Initiate pointer to the sparse matrix of the diffusion bilinear form
    SparseMatrix &BM_cdiff(varf_cdiff->SpMat());
    


    // Apply beta to the avg operator and its adjoint
    if (useGenResForm) { gen_res_manager.ImplementBeta(*BM_cavg); }
    
    // Obtain the block matrix for adding "a" into the averaging operator, and multiply it by beta and gamma
    SparseMatrix BM_avgavg = avgOp->GetAR_areas_SpMat(); BM_avgavg.Finalize();
    if (useGenResForm) {
        gen_res_manager.ConstructParamVecs(N_AR_global, AR_inds_loc2glob, true);
        if (active_adjoint == 1) { for (int i = 0; i < N_AR; i++) { BM_avgavg.Set(i, i, AR_areas[active_AR] * AR_areas[active_AR] / AR_areas[i]); } } // If adjoint, make all areas the active AR area
        gen_res_manager.ImplementGamma(BM_avgavg, AR_inds_loc2glob);
        closureOp_BLK.SetBlock(1, 1, &BM_avgavg);
    }
    


    // Create the transpose of the averaging operator matrix (this will be the bilinear form/operator of the unknown forcing terms in the mass conservation equation)
    TransposeOperator *BM_cavg_T = nullptr;
    BM_cavg_T = new TransposeOperator(BM_cavg);
    
    // Set the blocks of the closure system block operator
    closureOp_BLK.SetBlock(0, 0, &BM_cdiff);
    closureOp_BLK.SetBlock(0, 1, BM_cavg_T);
    closureOp_BLK.SetBlock(1, 0, BM_cavg);
    
    // Set the averaging operator to be 1 in the specified averaging region
    if (active_AR != -1 && forcing_function_file_name == "None") {
        b_BLK.GetBlock(1).Elem(active_AR) += AR_areas[active_AR] * resAvg_beta; }

    cout << "Complete." << endl;
    

    // ===============================================================
    //   Construct the preconditioner operator.
    // ===============================================================
    cout << globalVars.FILENAME << ":   Assembling the preconditioners... ";

    // Notes: - Here, we use Symmetric Gauss-Seidel to approximate the inverse of the "a" Schur complement.
    //
    //                      P = [ BM_cdiff      0    ] = [ BM_cdiff                 0                ]
    //                          [    0      S_actual ] = [    0      BM_cavg BM_cdiff^(-1) BM_cavg_T ]
    //
    //    We use this ->      ~ [ BM_cdiff  0 ] = [ BM_cdiff                    0                   ]
    //                          [     0     S ]   [    0      BM_cavg diag(BM_cdiff)^(-1) BM_cavg_T ]
    //
    SaddlePointBlockPreconditioner closurePC(block_offsets);
    if (is3D) { closurePC.BuildPreconditioner(BM_cdiff, *BM_cavg, "DSmoother", "GSSmoother"); }
    else { closurePC.BuildPreconditioner(BM_cdiff, *BM_cavg, "BlockILU", "BlockILU"); }
    
    cout << "Complete." << endl;
    


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
    cout << globalVars.FILENAME << ":   Solving" << endl;

    //MINRESSolver solver;
    GMRESSolver solver;
    solver.SetAbsTol(atol);
    solver.SetRelTol(rtol);
    solver.SetMaxIter(maxIter);
    solver.SetOperator(closureOp_BLK);
    solver.SetPreconditioner(closurePC);
    if (solve_mode == "serial") { solver.SetPrintLevel(1); }
    else { solver.SetPrintLevel(2); }
    solver.Mult(b_BLK, sol_BLK);

    if (solver.GetConverged()) {
        cout << "GMRES converged in " << solver.GetNumIterations()
            << " iterations with a residual norm of "
            << solver.GetFinalNorm() << ".\n";
    } else {
        cout << "GMRES did not converge in " << solver.GetNumIterations()
            << " iterations. Residual norm is " << solver.GetFinalNorm()
            << ".\n";
    }

    /*
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
    */


    // ===============================================================
    //   Extract and save the solutions (as gridfunctions) and the mesh
    // ===============================================================
    // Extract and save the closure variable solution
    GridFunction c_sol;
    c_sol.MakeRef(fespace_c, sol_BLK.GetBlock(0), 0);
    c_sol.Save((closure_output_dir + closure_output_file_name).c_str());
    
    // Save the mesh (only if solve mode is serial, or if rank is explicitly told to do so)
    if (useLocalMesh && saveLocalMesh) { mesh->Save(mesh_output_dir + submesh_output_dir_extn + mesh_output_file_name); }
    if (solve_mode == "serial" || save_mesh == 1) { mesh_manager.GetParentMesh()->Save(mesh_output_dir + mesh_output_file_name_prefix + mesh_output_file_name_suffix); }


    // ===============================================================
    //   Save the residuals in the residuals text file
    // ===============================================================
    // Define the sim_key and AR_number strings based on the simulation
    string sim_key;
    string AR_number;
    if (active_inlet == 1) { sim_key = "inlet"; AR_number = "0"; }
    else { sim_key = "avg_c"; AR_number = to_string(active_AR_global); }

    // Create the saver for the residuals and save them to the text file
    sol_BLK.GetBlock(1) *= -1.0; // The negative is because we solve with "a" on the left, but expect it on the right for the upscaled model
    ResultsSaver Saver;
    Saver.SetSolveMode(solve_mode);
    Saver.SetVariables(sim_key, recursive_iter, {"transport eq"}, {closure_output_file_name}, AR_number);
    /*
    if (useLocalMesh == 1) {
        bool obtainedResiduals = false;
        if (useGenResForm) { if (gen_res_manager.GetIsSparseCoefVecs()) { Saver.ObtainClosureResiduals(sol_BLK, {1}, {fespace_c->GetVDim()}, N_AR); obtainedResiduals = true; } }
        if (!obtainedResiduals) { Saver.ObtainClosureResiduals(sol_BLK, {1}, {fespace_c->GetVDim()}, N_AR_global, AR_inds_loc2glob); }
    }
    else { Saver.ObtainClosureResiduals(sol_BLK, {1}, {fespace_c->GetVDim()}, N_AR_global); }
    */
    
    // Save the residuals (i.e., the residuals of the mesh that was used: if a full mesh was used, a residual will be saved for each AR;
    // if a local mesh was used, a residual will be saved for each AR in the local mesh). Also save the corresponding AR indices (i.e.,
    // for a local mesh, these are the indices in AR_tags that identify which ARs were in the local mesh)
    Saver.ObtainClosureResiduals(sol_BLK, {1}, {fespace_c->GetVDim()}, N_AR);
    Saver.StoreResidualParameters("AR_inds", {{AR_inds_loc2glob}});
    if (AR_macroIDs.size() > 0) { Saver.StoreResidualParameters("AR_macroIDs", {{AR_macroIDs}}); }
    
    // If the generalized residual form is available, save the corresponding parameters for each AR in the used mesh
    // NOTE: gen_res_manager.GetIsSparseCoefVecs() should always be true now
    if (useGenResForm) {
        if (gen_res_manager.GetIsSparseCoefVecs()) {
            Saver.StoreResidualParameters("alpha", {{*gen_res_manager.GetParamVec("alpha")}});
            Saver.StoreResidualParameters("beta", {{*gen_res_manager.GetParamVec("beta")}});
            Saver.StoreResidualParameters("gamma", {{*gen_res_manager.GetParamVec("gamma")}});
            //Saver.StoreResidualParameters("AR_inds", {{*gen_res_manager.GetParamARInds()}});
        }
        else {
            Saver.StoreResidualParameters("alpha", {*gen_res_manager.GetParamVec("alpha")}, N_AR_global);
            Saver.StoreResidualParameters("beta", {*gen_res_manager.GetParamVec("beta")}, N_AR_global);
            Saver.StoreResidualParameters("gamma", {{*gen_res_manager.GetParamVec("gamma")}});
        }
    }
    Saver.SaveResidualDictionary(residual_output_file_path);


    // ===============================================================
    // 17. Free the used memory.
    // ===============================================================
    delete BM_cavg_T;
    delete avgOp;
    delete varf_cdiff;
    delete varf_force;
    // delete fespace_c; // Deleted with avgOp.
    delete fec_c;

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
