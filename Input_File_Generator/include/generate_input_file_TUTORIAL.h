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
#include "generate_input_file_TOOLS.h"



using namespace std;



// ===============================================================
//   Define a function for generating a default/tutorial config
//   file for the MoFA model formulation process.
// ===============================================================

void createTutorialUpscaledConfig(string &project_dir, string &config_file_name)
{
    
    // ===============================================================
    //   Define the config file dictionary.
    // ===============================================================
    JSONDict configFile;
    
    

    // ===============================================================
    //   Define the mesh dictionary.
    // ===============================================================
    JSONDict meshDict;

    JSONDict meshDict_output;
    meshDict_output["directory"] = project_dir + "mesh/"; // The mesh output directory 
    meshDict_output["file name"] = "upscaled_mesh.msh"; // The mesh file name
    meshDict["output path"] = &meshDict_output;

    JSONDict meshDict_info;
    meshDict_info["directory"] = meshDict_output["directory"]; // The mesh info output directory
    meshDict_info["file name"] = "upscaled_mesh_info.txt"; // The mesh info file name
    meshDict["info path"] = &meshDict_info;
    
    JSONDict meshDict_details;
    meshDict_details["gmsh model name"] = "2D_porous_media"; // The Gmsh model name (i.e., in Gmsh, you have to define model with a name)
    meshDict_details["geometry toggle"] = 1; // Declare whether to use the geometry in "geoDict". 0 = No, 1 = Yes
    meshDict_details["simulation type"] = "upscaled"; // The simulation type (either "upscaled" or "porescale")
    meshDict_details["is3D"] = 0; // Declare whether the mesh is 2D or 3D. 0 = 2D, 1 = 3D
    meshDict_details["isPeriodic"] = {1, 0, 0}; // Declare whether the mesh should be generated to handle periodicity B.C.s in the 3 directions. 0 = No, 1 = Yes
    meshDict_details["inlet side"] = "left"; // Declare which side of the domain is the inlet. Options: "left", "right", "top", "bottom", "none" (3D: "front", "back")
    meshDict_details["outlet side"] = "right"; // Declare which side of the domain is the outlet. Options: "left", "right", "top", "bottom", "none" (3D: "front", "back")
    meshDict_details["max element length"] = 0.02; // The maximum element size
    meshDict_details["min element length"] = 0.02; // The minimum element size
    meshDict["details"] = &meshDict_details;

    JSONDict meshDict_AR;
    meshDict_AR["epsilon"] = 0.1;
    JSONDict N_AR; N_AR["x"] = 10; N_AR["y"] = 1; // The number of averaging regions in each direction
    meshDict_AR["N_AR"] = &N_AR;
    JSONDict LDict; LDict["x"] = 10.0; LDict["y"] = 1.0; // The domain length in each direction
    meshDict_AR["L"] = &LDict;
    JSONDict ellDict; ellDict["x"] = 1.0; ellDict["y"] = 1.0; // The averaging region length in each direction
    meshDict_AR["ell"] = &ellDict;

    // Complex geometries can cause the pore space within an averaging region to be disconnected. In this case, each disconnected pore space is declared an independent AR.
    // The following options toggle whether to merge the newly created ARs with an adjacent "host" AR if 1.) the newly create ARs have a small area, 2.) the resulting
    // merged AR does not become too big, and 3.) the resulting merged AR does not have a large length-width ratio.
    meshDict_AR["merge ARs"] = 0; // Toggle whether merging averaging regions is allowed. 0 = No, 1 = Yes
    meshDict_AR["min area threshold"] = 0.2; // Declare the minimum area an AR can have without merging
    meshDict_AR["max AR length"] = 1.75; // Declare the maximum length a resulting merged AR can have
    meshDict_AR["max AR length ratio"] = 2.0; // Declare the maximum length-width ratio a resulting merged AR can have
    meshDict["AR"] = &meshDict_AR;

    JSONDict meshDict_geo;
    meshDict_geo["directory"] = project_dir + "data/"; // The mesh geometry directory
    meshDict_geo["file name"] = "Tutorial_Geometry.txt"; // The mesh geometry file name
    meshDict_geo["scale"] = 1.0; // The value with which to scale the geometry coordinates
    meshDict["geometry path"] = &meshDict_geo;

    JSONDict meshDict_cutPhysicalGroups;
    meshDict["cut physical groups"] = &meshDict_cutPhysicalGroups;

    configFile["mesh"] = &meshDict;



    // ===============================================================
    //   Define the stokes solver dictionary.
    // ===============================================================
    JSONDict stokesDict;

    JSONDict stokesDict_output;
    stokesDict_output["directory"] = project_dir + "output/stokes_solution/"; // The Stokes solver output directory
    stokesDict_output["velocity file name"] = "upscaled_u_sol.gf"; // The Stokes solver output file name for velocity
    stokesDict_output["pressure file name"] = "upscaled_p_sol.gf"; // The Stokes solver output file name for pressure
    stokesDict_output["mesh file name"] = "upscaled_mesh.mesh"; // The Stokes solver output file name for the mesh
    stokesDict["output path"] = &stokesDict_output;

    JSONDict stokesDict_sim;
    stokesDict_sim["u order"] = 2; // The element order to use for velocity
    stokesDict_sim["p order"] = 1; // The element order to use for pressure
    stokesDict_sim["isPeriodic"] = {1, 0, 0}; // Define whether to use periodic boundary conditions in the 3 directions. 0 = No, 1 = Yes
    stokesDict_sim["use inlet"] = 0; // Define whether to use apply an inlet condition to the inlet marked in the mesh. 0 = No, 1 = Yes
    stokesDict["simulation parameters"] = &stokesDict_sim;

    JSONDict stokesDict_phys;
    stokesDict_phys["A"] = 1.0; // The dimensionless number A (i.e., in the scaled Stokes equation)
    stokesDict_phys["body force"] = {337.84, 0.0, 0.0}; //337.84 // A constant body force to apply to the Stokes equation
    stokesDict_phys["u inlet max"] = 0.0; //0.1182 // The maximum velocity value for a parabolic inlet condition
    stokesDict["physics parameters"] = &stokesDict_phys;

    JSONDict stokesDict_solver;
    stokesDict_solver["max iterations"] = 10000; // The maximum number of iterations for the solver
    stokesDict_solver["rel tol"] = 1.0e-10; // The relative tolerance for the solver
    stokesDict_solver["abs tol"] = 0.0; // The absolute tolerance for the solver
    stokesDict["solver parameters"] = &stokesDict_solver;

    configFile["stokes"] = &stokesDict;

    

    // ===============================================================
    //   Define the closure problem solver dictionary.
    // ===============================================================
    JSONDict closureDict;

    JSONDict closureDict_sim;
    closureDict_sim["order"] = 2; // The element order to use
    closureDict_sim["active advection"] = 1; // Toggle for using a zero velocity field for advection. 0 = zero velocity, 1 = find and use previously solved Stokes solution
    closureDict_sim["isPeriodic"] = {0, 0, 0}; // Define whether to use periodic boundary conditions in the 3 directions. 0 = No, 1 = Yes
    closureDict_sim["use inlet"] = 1; // Define whether to solve the closure problems associated with the inlet condition at the inlet marked in the mesh. 0 = No, 1 = Yes (the answer is "no" for periodic BC; "yes" for Dirichlet BC)
    closureDict_sim["use reactions"] = 0; // Define whether to solve the closure problems associated with the reaction surfaces defined and marked in the mesh. 0 = No, 1 = Yes
    closureDict["simulation parameters"] = &closureDict_sim;

    JSONDict closureDict_res;
    closureDict_res["resAvg_alpha"] = 0.0;
    closureDict_res["resAvg_beta"] = 1.0;
    closureDict_res["resAvg_gamma"] = 0.0;
    closureDict["residual parameters"] = &closureDict_res;
    
    JSONDict closureDict_loc;
    closureDict_loc["use local mesh"] = 0; // Define whether the closure problems should be solved on the full mesh provided (= 0), or a mesh that is local to the AR/boundary forcing the closure problem
    closureDict_loc["N_neighbor_layers"] = 2; // Define the number of "neighbor layers" with which to define the local mesh. (e.g., 0 = consider just the forcing AR as the mesh, 1 = consider the domain in case 0 and all ARs that touch it, 2 = consider the domain in case 1 and all ARs that touch it, etc.)
    closureDict_loc["save local mesh"] = 0; // Toggle whether the local mesh for each closure problem should be saved or not. 0 = No, 1 = Yes
    closureDict["mesh localization parameters"] = &closureDict_loc;

    JSONDict closureDict_clos; // "Do not touch" Zone
    closureDict_clos["active averaging region"] = -1; // Toggle for which averaging region has the <\chi> = 1 condition. -1 = none, other identified the marked AR (NOTE: this is automatically changed by shell script)
    closureDict_clos["active inlet"] = 0; // Toggle for turning the inlet condition off or on. 0 = Off, 1 = On (NOTE: this is automatically changed by shell script)
    closureDict_clos["active reactions"] = {0}; // Toggle for turning the inlet condition off or on. 0 = Off, 1 = On (NOTE: this is automatically changed by shell script)
    closureDict_clos["max recursion iterations"] = 2; // The maximum number of recursive solves for the closure variable (usually due to the time derivative). NOTE: only the shell script used to run the ensemble of closure problems uses this
    closureDict["closure parameters"] = &closureDict_clos;

    JSONDict closureDict_phys;
    closureDict_phys["Pe_s"] = 1.0; // The scaled Peclet number Pe_s = Pe*\epsilon*\omega
    closureDict_phys["omega"] = 1.0; // The omega factor (i.e., Pe = \epsilon^{-1} * \omega^{-1})
    closureDict_phys["Da_s"] = {0.0}; // The scaled Peclet number Pe_s = Pe*\epsilon*\omega
    closureDict["physics parameters"] = &closureDict_phys;

    JSONDict closureDict_solver;
    closureDict_solver["max iterations"] = 10000; // The maximum number of iterations for the solver
    closureDict_solver["rel tol"] = 1.0e-13; // The relative tolerance for the solver
    closureDict_solver["abs tol"] = 1.0e-8; // The absolute tolerance for the solver
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

    configFile["scalar closure"] = &closureDict;



    // ===============================================================
    //   Define the upscaled solver dictionary.
    // ===============================================================
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
    upscaledDict_sim["N time steps"] = 1000; // The number of time steps
    upscaledDict_sim["dt"] = 0.0001; // The time step size
    upscaledDict_sim["output interval"] = 100; // The number of time steps before another solution is saved to the output
    upscaledDict_sim["BC frequency scale"] = 1.0; // This value multiplies the boundary condition frequency for faster/slower oscillation. Maximum magnitude of dc/dt at the BC is this number times \epsilon^(-2). (\epsilon^(-2) is typically the maximum BC oscillation frequency allowed for by the model)
    upscaledDict["simulation parameters"] = &upscaledDict_sim;

    configFile["upscaled"] = &upscaledDict;



    // ===============================================================
    //   Define the error calculator dictionary.
    // ===============================================================
    JSONDict errorDict;
    istringstream dict = createErrorDict(project_dir, "transport");
    errorDict.loadFromStream( dict );
    configFile["error calc"] = &errorDict;



    // ===============================================================
    //   Save the config file.
    // ===============================================================
    configFile.saveToFile(project_dir + "config/" + config_file_name);
}










// ===============================================================
//   Define a function for generating a default/tutorial config
//   file for the porescale model formulation process.
// ===============================================================

void createTutorialPorescaleConfig(string &project_dir, string &config_file_name)
{
    
    // ===============================================================
    //   Define the config file dictionary.
    // ===============================================================
    JSONDict configFile;
    

    
    // ===============================================================
    //   Define the mesh dictionary.
    // ===============================================================
    JSONDict meshDict;

    JSONDict meshDict_output;
    meshDict_output["directory"] = project_dir + "mesh/"; // The mesh output directory 
    meshDict_output["file name"] = "porescale_mesh.msh"; // The mesh file name
    meshDict["output path"] = &meshDict_output;

    JSONDict meshDict_info;
    meshDict_info["directory"] = meshDict_output["directory"]; //project_dir + "mesh/"; // The mesh info output directory
    meshDict_info["file name"] = "porescale_mesh_info.txt"; // The mesh info file name
    meshDict["info path"] = &meshDict_info;
    
    JSONDict meshDict_details;
    meshDict_details["gmsh model name"] = "2D_porous_media"; // The Gmsh model name (i.e., in Gmsh, you have to define model with a name)
    meshDict_details["geometry toggle"] = 1; // Declare whether to use the geometry in "geoDict". 0 = No, 1 = Yes
    meshDict_details["simulation type"] = "porescale"; // The simulation type (either "upscaled" or "porescale")
    meshDict_details["is3D"] = 0; // Declare whether the mesh is 2D or 3D. 0 = 2D, 1 = 3D
    meshDict_details["isPeriodic"] = {1, 0, 0}; // Declare whether the mesh should be generated to handle periodicity B.C.s in the 3 directions. 0 = No, 1 = Yes
    meshDict_details["inlet side"] = "left"; // Declare which side of the domain is the inlet. Options: "left", "right", "top", "bottom", "none" (3D: "front", "back")
    meshDict_details["outlet side"] = "right"; // Declare which side of the domain is the outlet. Options: "left", "right", "top", "bottom", "none" (3D: "front", "back")
    meshDict_details["max element length"] = 0.002; // The maximum element size
    meshDict_details["min element length"] = 0.002; // The minimum element size
    meshDict["details"] = &meshDict_details;

    JSONDict meshDict_AR;
    meshDict_AR["epsilon"] = 0.1;
    JSONDict N_AR; N_AR["x"] = 10; N_AR["y"] = 1; // The number of averaging regions in each direction
    meshDict_AR["N_AR"] = &N_AR;
    JSONDict LDict; LDict["x"] = 1.0; LDict["y"] = 0.1; // The domain length in each direction
    meshDict_AR["L"] = &LDict;
    JSONDict ellDict; ellDict["x"] = 0.1; ellDict["y"] = 0.1; // The averaging region length in each direction
    meshDict_AR["ell"] = &ellDict;
    meshDict["AR"] = &meshDict_AR;

    // Complex geometries can cause the pore space within an averaging region to be disconnected. In this case, each disconnected pore space is declared an independent AR.
    // The following options toggle whether to merge the newly created ARs with an adjacent "host" AR if 1.) the newly create ARs have a small area, 2.) the resulting
    // merged AR does not become too big, and 3.) the resulting merged AR does not have a large length-width ratio.
    meshDict_AR["merge ARs"] = 0; // Toggle whether merging averaging regions is allowed. 0 = No, 1 = Yes
    meshDict_AR["min area threshold"] = 0.2; // Declare the minimum area an AR can have without merging
    meshDict_AR["max AR length"] = 1.75; // Declare the maximum length a resulting merged AR can have
    meshDict_AR["max AR length ratio"] = 2.0; // Declare the maximum length-width ratio a resulting merged AR can have
    meshDict["AR"] = &meshDict_AR;

    JSONDict meshDict_geo;    
    meshDict["geometry path"] = &meshDict_geo;
    meshDict_geo["directory"] = project_dir + "data/"; // The mesh geometry directory
    meshDict_geo["file name"] = "Tutorial_Geometry.txt"; // The mesh geometry file name
    meshDict_geo["scale"] = 0.1; // The value with which to scale the geometry coordinates
    meshDict["geometry path"] = &meshDict_geo;
    
    JSONDict meshDict_cutPhysicalGroups;
    meshDict["cut physical groups"] = &meshDict_cutPhysicalGroups;

    configFile["mesh"] = &meshDict;



    // ===============================================================
    //   Define the stokes solver dictionary.
    // ===============================================================
    JSONDict stokesDict;

    JSONDict stokesDict_output;
    stokesDict_output["directory"] = project_dir + "output/stokes_solution/"; // The Stokes solver output directory
    stokesDict_output["velocity file name"] = "porescale_u_sol.gf"; // The Stokes solver output file name for velocity
    stokesDict_output["pressure file name"] = "porescale_p_sol.gf"; // The Stokes solver output file name for pressure
    stokesDict_output["mesh file name"] = "porescale_mesh.mesh"; // The Stokes solver output file name for the mesh
    stokesDict["output path"] = &stokesDict_output;

    JSONDict stokesDict_sim;
    stokesDict_sim["u order"] = 2; // The element order to use for velocity
    stokesDict_sim["p order"] = 1; // The element order to use for pressure
    stokesDict_sim["isPeriodic"] = {1, 0, 0}; // Define whether to use periodic boundary conditions in the 3 directions. 0 = No, 1 = Yes
    stokesDict_sim["use inlet"] = 0; // Define whether to use apply an inlet condition to the inlet marked in the mesh. 0 = No, 1 = Yes
    stokesDict["simulation parameters"] = &stokesDict_sim;

    JSONDict stokesDict_phys;
    stokesDict_phys["A"] = 0.01; // The dimensionless number A (i.e., in the scaled Stokes equation)
    stokesDict_phys["body force"] = {337.84, 0.0, 0.0}; // 337.84 A constant body force to apply to the Stokes equation
    stokesDict_phys["u inlet max"] = 0.0; //0.1182; // The maximum velocity value for a parabolic inlet condition
    stokesDict["physics parameters"] = &stokesDict_phys;

    JSONDict stokesDict_solver;
    stokesDict_solver["max iterations"] = 10000; // The maximum number of iterations for the solver
    stokesDict_solver["rel tol"] = 1.0e-10; // The relative tolerance for the solver
    stokesDict_solver["abs tol"] = 0.0; // The absolute tolerance for the solver
    stokesDict["solver parameters"] = &stokesDict_solver;

    configFile["stokes"] = &stokesDict;

    

    // ===============================================================
    //   Define the porescale transport solver dictionary.
    // ===============================================================
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
    porescaleDict_phys["Pe"] = 10.0; // The Peclet number
    porescaleDict_phys["omega"] = 1.0; // The omega factor
    porescaleDict["physical parameters"] = &porescaleDict_phys;
    
    JSONDict porescaleDict_sim;
    porescaleDict_sim["order"] = 1; // The element order to use for the concentration
    porescaleDict_sim["active advection"] = 1; // Toggle for using a zero velocity field for advection. 0 = zero velocity, 1 = find and use previously solved Stokes solution
    porescaleDict_sim["isPeriodic"] = {0, 0, 0}; // Define whether to use periodic boundary conditions in the 3 directions. 0 = No, 1 = Yes

    porescaleDict_sim["N time steps"] = 1000; // The number of time steps
    porescaleDict_sim["dt"] = 0.0001; // The time step size
    porescaleDict_sim["output interval"] = 100; // The number of time steps before another solution is saved to the output
    porescaleDict_sim["BC frequency scale"] = 1.0; // This value multiplies the boundary condition frequency for faster/slower oscillation. Maximum magnitude of dc/dt at the BC is this number times \epsilon^(-2). (\epsilon^(-2) is typically the maximum BC oscillation frequency allowed for by the model)
    porescaleDict["simulation parameters"] = &porescaleDict_sim;
    
    configFile["scalar porescale"] = &porescaleDict;

    

    // ===============================================================
    //   Save the config file.
    // ===============================================================
    configFile.saveToFile(project_dir + "config/" + config_file_name);
}

