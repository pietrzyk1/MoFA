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

#pragma once
#include "JSON_IO.h"



using namespace std;



// ========================================
// For creating the Mesh JSON dictionary
// ========================================
istringstream createMeshDict(const string &project_dir, const string mesh_file_name, const string sim_type,
    const int geo_toggle, const int is3D, const vector<int> periodicity, const vector<double> min_max_elem_size,
    const double epsilon, const vector<int> N_AR, const vector<double> L, const vector<double> ell,
    const int merge_ARs, const double min_area_thresh, const double max_AR_length, const double max_AR_length_ratio,
    const string geo_dir, const string geo_file_name, const double geo_scale,
    const vector<vector<int>> pg_cut_inds, const vector<string> pg_cut_names)
{
    JSONDict meshDict;

    JSONDict meshDict_output;
    meshDict_output["directory"] = project_dir + "mesh/"; // The mesh output directory 
    meshDict_output["file name"] = mesh_file_name + ".msh"; // The mesh file name
    meshDict["output path"] = &meshDict_output;

    JSONDict meshDict_info;
    meshDict_info["directory"] = meshDict_output["directory"]; // The mesh info output directory
    meshDict_info["file name"] = mesh_file_name + "_info.txt"; // The mesh info file name
    meshDict["info path"] = &meshDict_info;
    
    JSONDict meshDict_details;
    meshDict_details["gmsh model name"] = "2D_porous_media"; // The Gmsh model name (i.e., in Gmsh, you have to define model with a name)
    meshDict_details["geometry toggle"] = geo_toggle; // Declare whether to use the geometry in "geoDict". 0 = No, 1 = Yes
    meshDict_details["simulation type"] = sim_type; // The simulation type (either "upscaled" or "porescale")
    meshDict_details["is3D"] = is3D; // Declare whether the mesh is 2D or 3D. 0 = 2D, 1 = 3D
    meshDict_details["isPeriodic"] = periodicity; // Declare whether the mesh should be generated to handle periodicity B.C.s in the 3 directions. 0 = No, 1 = Yes
    meshDict_details["inlet side"] = "left"; // Declare which side of the domain is the inlet. Options: "left", "right", "top", "bottom", "none" (3D: "front", "back")
    meshDict_details["outlet side"] = "right"; // Declare which side of the domain is the outlet. Options: "left", "right", "top", "bottom", "none" (3D: "front", "back")
    meshDict_details["max element length"] = min_max_elem_size[1]; // The maximum element size
    meshDict_details["min element length"] = min_max_elem_size[0]; // The minimum element size
    meshDict["details"] = &meshDict_details;
    
    JSONDict meshDict_AR;
    meshDict_AR["epsilon"] = epsilon;
    
    JSONDict N_AR_dict; N_AR_dict["x"] = N_AR[0]; N_AR_dict["y"] = N_AR[1]; // The number of averaging regions in each direction
    if (N_AR.size() == 3) { N_AR_dict["z"] = N_AR[2]; }
    meshDict_AR["N_AR"] = &N_AR_dict;
    
    JSONDict LDict; LDict["x"] = L[0]; LDict["y"] = L[1]; // The domain length in each direction
    if (L.size() == 3) { LDict["z"] = L[2]; }
    meshDict_AR["L"] = &LDict;

    JSONDict ellDict; ellDict["x"] = ell[0]; ellDict["y"] = ell[1]; // The averaging region length in each direction
    if (ell.size() == 3) { ellDict["z"] = ell[2]; }    
    meshDict_AR["ell"] = &ellDict;

    // Complex geometries can cause the pore space within an averaging region to be disconnected. In this case, each disconnected pore space is declared an independent AR.
    // The following options toggle whether to merge the newly created ARs with an adjacent "host" AR if 1.) the newly create ARs have a small area, 2.) the resulting
    // merged AR does not become too big, and 3.) the resulting merged AR does not have a large length-width ratio.
    meshDict_AR["merge ARs"] = merge_ARs; // Toggle whether merging averaging regions is allowed. 0 = No, 1 = Yes
    meshDict_AR["min area threshold"] = min_area_thresh; // Declare the minimum area an AR can have without merging
    meshDict_AR["max AR length"] = max_AR_length; // Declare the maximum length a resulting merged AR can have
    meshDict_AR["max AR length ratio"] = max_AR_length_ratio; // Declare the maximum length-width ratio a resulting merged AR can have
    meshDict["AR"] = &meshDict_AR;
    
    JSONDict meshDict_geo;
    meshDict_geo["directory"] = geo_dir; // "Geometry/Directory/"; // The mesh geometry directory
    meshDict_geo["file name"] = geo_file_name; // "Geometry_File_Name.txt"; // The mesh geometry file name
    meshDict_geo["scale"] = geo_scale; // The value with which to scale the geometry coordinates
    meshDict["geometry path"] = &meshDict_geo;

    JSONDict meshDict_cutPhysicalGroups;
    if (pg_cut_names.size() != 0) { meshDict_cutPhysicalGroups["physical group cut names"] = pg_cut_names; }
    if (pg_cut_inds.size() != 0) { meshDict_cutPhysicalGroups["physical group cut indices"] = pg_cut_inds; }
    meshDict["cut physical groups"] = &meshDict_cutPhysicalGroups;

    return meshDict.saveToStream();
}





// ========================================
// For creating the Stokes Solver JSON dictionary
// ========================================
istringstream createStokesSolverDict(const string &project_dir,
    const string velocity_file_name, const string pressure_file_name, const string mesh_file_name,
    const int u_order, const int p_order, const vector<int> periodicity, const double A_coef,
    const int use_inlet, const vector<double> body_force, const double u_inlet_max,
    const int max_iter, const double rel_tol, const double abs_tol)
{
    JSONDict stokesDict;
    
    JSONDict stokesDict_output;
    stokesDict_output["directory"] = project_dir + "output/stokes_solution/"; // The Stokes solver output directory
    stokesDict_output["velocity file name"] = velocity_file_name; // The Stokes solver output file name for velocity
    stokesDict_output["pressure file name"] = pressure_file_name; // The Stokes solver output file name for pressure
    stokesDict_output["mesh file name"] = mesh_file_name; // The Stokes solver output file name for the mesh
    stokesDict["output path"] = &stokesDict_output;
    
    JSONDict stokesDict_sim;
    stokesDict_sim["u order"] = u_order; // The element order to use for velocity
    stokesDict_sim["p order"] = p_order; // The element order to use for pressure
    stokesDict_sim["isPeriodic"] = periodicity; // Define whether to use periodic boundary conditions in the 3 directions. 0 = No, 1 = Yes
    stokesDict_sim["use inlet"] = use_inlet; // Define whether to use apply an inlet condition to the inlet marked in the mesh. 0 = No, 1 = Yes
    stokesDict["simulation parameters"] = &stokesDict_sim;
    
    JSONDict stokesDict_phys;
    stokesDict_phys["A"] = A_coef; // The dimensionless number A (i.e., in the scaled Stokes equation)
    stokesDict_phys["body force"] = body_force; //337.84 // A constant body force to apply to the Stokes equation
    stokesDict_phys["u inlet max"] = u_inlet_max; //0.1182 // The maximum velocity value for a parabolic inlet condition
    stokesDict["physics parameters"] = &stokesDict_phys;
    
    JSONDict stokesDict_solver;
    stokesDict_solver["max iterations"] = max_iter; // The maximum number of iterations for the solver
    stokesDict_solver["rel tol"] = rel_tol; // The relative tolerance for the solver
    stokesDict_solver["abs tol"] = abs_tol; // The absolute tolerance for the solver
    stokesDict["solver parameters"] = &stokesDict_solver;

    return stokesDict.saveToStream();
}





// ========================================
// For creating the Scalar Transport Closure Solver JSON dictionary
// ========================================
istringstream createTransportClosureDict(const string &project_dir,
    const int order, const vector<int> isPeriodic, const int active_advection, const int useInlet, const int useReactions,
    const double Pe_s, const double omega, const vector<double> Da_s,
    const int useLocalMesh, const int N_neighbor_layers, const int saveLocalMesh,
    const double resAvg_alpha, const double resAvg_beta, const double resAvg_gamma,
    const int max_iter, const double rel_tol, const double abs_tol,
    const int importFluidMesh = 0, const int maxRecursionIter = 2)
{
    JSONDict closureDict;

    JSONDict closureDict_sim;
    closureDict_sim["order"] = order; // The element order to use
    closureDict_sim["active advection"] = active_advection; // Toggle for using a zero velocity field for advection. 0 = zero velocity, 1 = find and use previously solved Stokes solution
    closureDict_sim["isPeriodic"] = isPeriodic; // Define whether to use periodic boundary conditions in the 3 directions. 0 = No, 1 = Yes
    closureDict_sim["use inlet"] = useInlet; // Define whether to solve the closure problems associated with the inlet condition at the inlet marked in the mesh. 0 = No, 1 = Yes (the answer is "no" for periodic BC; "yes" for Dirichlet BC)
    closureDict_sim["use reactions"] = useReactions; // Define whether to solve the closure problems associated with the reaction surfaces defined and marked in the mesh. 0 = No, 1 = Yes
    closureDict_sim["import fluid mesh"] = importFluidMesh; // Define whether the solver should import the mesh used for solving the Stokes problem (for fluid velocity) or use the mesh create from gmsh. 0 = use gmsh mesh, 1 = use MFEM mesh saved with the fluid velocity solution
    closureDict["simulation parameters"] = &closureDict_sim;

    JSONDict closureDict_res;
    closureDict_res["resAvg_alpha"] = resAvg_alpha;
    closureDict_res["resAvg_beta"] = resAvg_beta;
    closureDict_res["resAvg_gamma"] = resAvg_gamma;
    closureDict["residual parameters"] = &closureDict_res;
    
    JSONDict closureDict_loc;
    closureDict_loc["use local mesh"] = useLocalMesh; // Define whether the closure problems should be solved on the full mesh provided (= 0), or a mesh that is local to the AR/boundary forcing the closure problem
    closureDict_loc["N_neighbor_layers"] = N_neighbor_layers; // Define the number of "neighbor layers" with which to define the local mesh. (e.g., 0 = consider just the forcing AR as the mesh, 1 = consider the domain in case 0 and all ARs that touch it, 2 = consider the domain in case 1 and all ARs that touch it, etc.)
    closureDict_loc["save local mesh"] = saveLocalMesh; // Toggle whether the local mesh for each closure problem should be saved or not. 0 = No, 1 = Yes
    closureDict["mesh localization parameters"] = &closureDict_loc;

    JSONDict closureDict_clos; // "Do not touch" Zone
    closureDict_clos["active averaging region"] = -1; // Toggle for which averaging region has the <\chi> = 1 condition. -1 = none, other identified the marked AR (NOTE: this is automatically changed by shell script)
    closureDict_clos["active inlet"] = 0; // Toggle for turning the inlet condition off or on. 0 = Off, 1 = On (NOTE: this is automatically changed by shell script)
    closureDict_clos["active reactions"] = {0}; // Toggle for turning the inlet condition off or on. 0 = Off, 1 = On (NOTE: this is automatically changed by shell script)
    closureDict_clos["max recursion iterations"] = maxRecursionIter; // The maximum number of recursive solves for the closure variable (usually due to the time derivative). NOTE: only the shell script used to run the ensemble of closure problems uses this
    closureDict["closure parameters"] = &closureDict_clos;

    JSONDict closureDict_phys;
    closureDict_phys["Pe_s"] = Pe_s; // The scaled Peclet number Pe_s = Pe*\epsilon*\omega
    closureDict_phys["omega"] = omega; // The omega factor (i.e., Pe = \epsilon^{-1} * \omega^{-1})
    closureDict_phys["Da_s"] = Da_s; // The scaled Peclet number Pe_s = Pe*\epsilon*\omega
    closureDict["physics parameters"] = &closureDict_phys;

    JSONDict closureDict_solver;
    closureDict_solver["max iterations"] = max_iter; // The maximum number of iterations for the solver
    closureDict_solver["rel tol"] = rel_tol; // The relative tolerance for the solver
    closureDict_solver["abs tol"] = abs_tol; // The absolute tolerance for the solver
    closureDict["solver parameters"] = &closureDict_solver;

    string closureDict_sol_dir = project_dir + "output/scalar_closure_solution/";
    
    JSONDict closureDict_residual;
    closureDict_residual["directory"] = closureDict_sol_dir + "closure_residual/"; // The closure residual output directory
    closureDict_residual["file name prefix"] = "a_sol"; // The prefix of the closure residual file name
    closureDict_residual["file name suffix"] = ".txt"; // The suffix fo the closure residual file name
    closureDict["residual path"] = &closureDict_residual;

    JSONDict closureDict_path;
    closureDict_path["directory"] = closureDict_sol_dir + "closure_solution/"; // The closure solution output directory
    closureDict_path["file name prefix"] = "chi"; // The prefix of the closure solution output file name
    closureDict_path["file name suffix"] = ".gf"; // The suffix of the closure solution output file name
    closureDict["closure path"] = &closureDict_path;

    JSONDict closureDict_mesh;
    closureDict_mesh["directory"] = closureDict_path["directory"]; // The closure solution mesh output directory
    //closureDict_mesh["file name"] = "closure_mesh.mesh"; // The closure solution mesh output file name
    closureDict_mesh["file name prefix"] = "closure_mesh"; // The prefix of the closure solution mesh output file name
    closureDict_mesh["file name suffix"] = ".mesh"; // The suffix of the closure solution mesh output file name
    closureDict["mesh output path"] = &closureDict_mesh;

    JSONDict closureDict_forcing;
    closureDict_forcing["directory"] = closureDict_sol_dir + "closure_solution/"; // The applied forcing function directory
    closureDict_forcing["file name"] = "None"; // The applied forcing function file name
    closureDict["forcing function path"] = &closureDict_forcing;

    return closureDict.saveToStream();
}





// ========================================
// For creating the Scalar Transport Upscaled Model Solver JSON dictionary
// ========================================
istringstream createTransportUpscaledDict(const string &project_dir,
    const double inletFreq,
    const int N_time_steps, const double dt, const int output_interval)
{
    JSONDict upscaledDict;

    JSONDict upscaledDict_output_path;
    upscaledDict_output_path["directory"] = project_dir + "output/upscaled_solution/"; // The output directory for the gridfunctions of the upscaled solution
    upscaledDict_output_path["file name prefix"] = "upscaled_sol_"; // The prefix of the output gridfunctions of the upscaled solution
    upscaledDict_output_path["file name suffix"] = ".gf"; // The suffix of the output gridfunctions of the upscaled solution
    upscaledDict["output path"] = &upscaledDict_output_path;

    JSONDict upscaledDict_avg_output_path;
    upscaledDict_avg_output_path["directory"] = upscaledDict_output_path["directory"]; // The output directory for the upscaled solution (i.e., the averaged solution)
    upscaledDict_avg_output_path["file name"] = "upscaled_sol.txt"; // The file name for the upscaled solution (i.e., the averaged solution)
    upscaledDict["average output path"] = &upscaledDict_avg_output_path;

    JSONDict upscaledDict_mesh_output;
    upscaledDict_mesh_output["directory"] = upscaledDict_output_path["directory"]; // The upscaled solution mesh output directory
    upscaledDict_mesh_output["file name"] = "upscaled_mesh.mesh"; // The upscaled solution mesh output file name
    upscaledDict["mesh output path"] = &upscaledDict_mesh_output;

    JSONDict upscaledDict_sim;
    upscaledDict_sim["N time steps"] = N_time_steps; // The number of time steps
    upscaledDict_sim["dt"] = dt; // The time step size
    upscaledDict_sim["output interval"] = output_interval; // The number of time steps before another solution is saved to the output
    upscaledDict_sim["BC frequency scale"] = inletFreq; // This value multiplies the boundary condition frequency for faster/slower oscillation. Maximum magnitude of dc/dt at the BC is this number times \epsilon^(-2). (\epsilon^(-2) is typically the maximum BC oscillation frequency allowed for by the model)
    upscaledDict["simulation parameters"] = &upscaledDict_sim;

    return upscaledDict.saveToStream();
}





// ========================================
// For creating the Scalar Transport Porescale Solver JSON dictionary
// ========================================
istringstream createTransportPorescaleDict(const string &project_dir,
    const int order, const vector<int> isPeriodic, const int active_advection, const int useInlet, const int useReactions,
    const double Pe, const double omega, const vector<double> Da_s, const double inletFreq,
    const int N_time_steps, const double dt, const int output_interval)
{
    JSONDict porescaleDict;
    
    JSONDict porescaleDict_output_path;
    porescaleDict_output_path["directory"] = project_dir + "output/porescale_solution/"; // The porescale solution output directory
    porescaleDict_output_path["file name prefix"] = "c_"; // The porescale solution file name prefix
    porescaleDict_output_path["file name suffix"] = ".gf"; // The porescale solution file name suffix
    porescaleDict["output path"] = &porescaleDict_output_path;

    JSONDict porescaleDict_avg_output_path;
    porescaleDict_avg_output_path["directory"] = project_dir + "output/avg_porescale_solution/"; // The average porescale solution output directory
    porescaleDict_avg_output_path["file name"] = "avg_sol.txt"; // The average porescale solution output file name
    porescaleDict["average output path"] = &porescaleDict_avg_output_path;

    JSONDict porescaleDict_mesh_output;
    porescaleDict_mesh_output["directory"] = porescaleDict_output_path["directory"]; // The porescale solution mesh output directory
    porescaleDict_mesh_output["file name"] = "porescale_mesh.mesh"; // The porescale solution mesh output file name
    porescaleDict["mesh output path"] = &porescaleDict_mesh_output;
    
    JSONDict porescaleDict_phys;
    porescaleDict_phys["Pe"] = Pe; // The Peclet number
    porescaleDict_phys["omega"] = omega; // The omega factor
    porescaleDict["physical parameters"] = &porescaleDict_phys;
    
    JSONDict porescaleDict_sim;
    porescaleDict_sim["order"] = order; // The element order to use for the concentration
    porescaleDict_sim["active advection"] = active_advection; // Toggle for using a zero velocity field for advection. 0 = zero velocity, 1 = find and use previously solved Stokes solution
    porescaleDict_sim["isPeriodic"] = isPeriodic; // Define whether to use periodic boundary conditions in the 3 directions. 0 = No, 1 = Yes

    porescaleDict_sim["N time steps"] = N_time_steps; // The number of time steps
    porescaleDict_sim["dt"] = dt; // The time step size
    porescaleDict_sim["output interval"] = output_interval; // The number of time steps before another solution is saved to the output
    porescaleDict_sim["BC frequency scale"] = inletFreq; // This value multiplies the boundary condition frequency for faster/slower oscillation. Maximum magnitude of dc/dt at the BC is this number times \epsilon^(-2). (\epsilon^(-2) is typically the maximum BC oscillation frequency allowed for by the model)
    porescaleDict["simulation parameters"] = &porescaleDict_sim;
    
    return porescaleDict.saveToStream();
}





// ========================================
// For creating the Error JSON dictionary
// ========================================
istringstream createErrorDict(const string &project_dir, string sim_type)
{
    JSONDict errorDict;

    JSONDict errorDict_output;
    errorDict_output["directory"] = project_dir + "output/absolute_error/"; // The absolute error output directory
    errorDict_output["file name"] = "abs_error.txt"; // The absolute error output file name
    errorDict["output path"] = &errorDict_output;

    JSONDict errorDict_avgpore;
    if (sim_type == "transport") { errorDict_avgpore["directory"] = project_dir + "output/avg_porescale_solution/"; } // The directory where the average porescale solution is stored
    else { errorDict_avgpore["directory"] = project_dir + "output/avg_unsteady_stokes_solution/"; } // The directory where the average porescale solution is stored
    errorDict_avgpore["file name"] = "avg_sol.txt"; // The average porescale solution file name
    errorDict["avg porescale path"] = &errorDict_avgpore;

    JSONDict errorDict_upscaled;
    if (sim_type == "transport") { errorDict_upscaled["directory"] = project_dir + "output/upscaled_solution/"; } // The directory where the upscaled solution is stored
    else { errorDict_upscaled["directory"] = project_dir + "output/upscaled_unsteady_stokes_solution/"; } // The directory where the upscaled solution is stored
    errorDict_upscaled["file name"] = "upscaled_sol.txt"; // The upscaled solution file name
    errorDict["upscaled path"] = &errorDict_upscaled;

    return errorDict.saveToStream();
}

