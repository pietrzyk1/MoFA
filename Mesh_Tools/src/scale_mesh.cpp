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

//                                      Mesh Scaler
// Description: This code is used to create a pore-scale mesh, and a file for it's
//              corresponding information from an upscaled mesh. Generally speaking,
//              this means this code essentially scales the upscaled mesh information
//              to get the pore-scale mesh and information so that a separate
//              pore-scale mesh does not need to be created.
//

// ========================================
// Variable Nomenclature
// ========================================  
// pt(s) : point(s)
// ln(s) : lines(s)
// srfc(s) : surfaces(s)
// vol(s) : volume(s)
// AR(s) : averaging region(s)
// uc : uncut
// c : cut
// tool : the geometry used to cut a blank domain (which is created by contiguously placed ARs)
// geo : refers to the geometry; variables that end in "geo" typically have the structure [AR, lines, points, coords],
//       where the data (at the ends of the arrays) is coordinate values

#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include <gmsh.h>
#include "JSON_IO.h"
#include "gmsh_utils.h"





// ========================================
// Define global variables
// ========================================  
// Define a structure with global parameters
struct Params
{
    string FILENAME = "scale_mesh.cpp";
};
static Params globalVars; 





int main(int argc, char **argv)
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string config_path = "./";

    int is3D = 0; // If the domain is 3D, = 1. If 2D, = 0
    double scale = 1.0;

    string problem_type = "porescale";
    
    string orig_mesh_dir = "./";
    string orig_mesh_file_name = "mesh.msh";
    string orig_mesh_info_dir = "./";
    string orig_mesh_info_file_name = "mesh_info.txt";

    string scaled_mesh_dir = "./";
    string scaled_mesh_file_name = "porescale_mesh.msh";
    string scaled_mesh_info_dir = "./";
    string scaled_mesh_info_file_name = "porescale_mesh_info.txt";


    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++) {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc) {
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
        
        JSONDict sub_dict = *mesh_dict["details"];
        sub_dict.getValue("is3D", is3D);

        sub_dict = *mesh_dict["AR"];
        sub_dict.getValue("epsilon", scale);
        
        sub_dict = *mesh_dict["output path"];
        sub_dict.getValue("directory", orig_mesh_dir);
        sub_dict.getValue("file name", orig_mesh_file_name);

        sub_dict = *mesh_dict["info path"];
        sub_dict.getValue("directory", orig_mesh_info_dir);
        sub_dict.getValue("file name", orig_mesh_info_file_name);
    }
    else
    {
        cerr << globalVars.FILENAME << ": main(): Error in loading config file." << endl;
        exit(1);
    }
    
    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string orig_mesh_file_path = orig_mesh_dir + orig_mesh_file_name;
    string orig_mesh_info_file_path = orig_mesh_info_dir + orig_mesh_info_file_name;
    
    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    //OptionsParser args(argc, argv);
    //args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    //args.ParseCheck();

    scaled_mesh_dir = orig_mesh_dir;
    scaled_mesh_info_dir = orig_mesh_info_dir;
    string scaled_mesh_file_path = scaled_mesh_dir + scaled_mesh_file_name;
    string scaled_mesh_info_file_path = scaled_mesh_info_dir + scaled_mesh_info_file_name;
    

    // ========================================
    //   Initialize gmsh and open the mesh file
    // ========================================
    // Gmsh must be initialized before using any functions in the C++ API
    gmsh::initialize();
    
    // Open the mesh file to be scaled
    gmsh::open(orig_mesh_file_path);
    

    // ========================================
    //   Retrieve and scale all node corrdinates in the mesh
    // ========================================
    std::vector<std::size_t> nodeTags;
    std::vector<double> coords;
    std::vector<double> parametricCoords;

    // Retrieve all node tags and coordinates in the mesh (independent of model entities)
    gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);

    // If there are parametric nodes, throw an error. Not sure how to handle these
    if (parametricCoords.size() != 0) {
        cerr << globalVars.FILENAME << ": CRITICAL ERROR: parametricCoords.size() > 0. Have not considered parametric nodes yet." << endl;
        exit(1);
    }

    // Scale the node coordinates. Variable "coords" is [x0,y0,z0, x1,y1,z1, ...]
    for (int i = 0; i < coords.size(); i++) {
        if (i % 3 == 0) { cout << globalVars.FILENAME << ": Scaling coordinate " << i/3 << " out of " << coords.size()/3 << "." << endl; }
        coords[i] *= scale;
    }

    // Set the node coordinates back in the mesh
    for (int i = 0; i < nodeTags.size(); i++) {
        cout << globalVars.FILENAME << ": Setting node coordinate " << i << " out of " << nodeTags.size() << "." << endl;
        gmsh::model::mesh::setNode(nodeTags[i], {coords[3*i + 0], coords[3*i + 1], coords[3*i + 2]}, parametricCoords);
    }

    
    // ========================================
    // Save the mesh
    // ========================================
    cout << "scale_mesh.cpp: Saving mesh..." << endl;
    
    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2); // Make the saved version 2.2 so that it is MFEM compatible
    gmsh::write(scaled_mesh_file_path); // Save mesh with defined name
    cout << "    Complete!" << endl;
    

    // ========================================
    //   Create and save the scaled mesh information
    // ========================================
    cout << "scale_mesh.cpp: Preparing mesh info file..." << endl;

    // Load the mesh information
    JSONDict mesh_info; mesh_info.loadFromFile(orig_mesh_info_file_path);

    
    // Get the AR dictionary from the mesh information
    JSONDict AR_dict = *mesh_info["AR"];

    // Get the pore areas from the mesh information, scale them, and reassign them
    {
        std::vector<double> pore_areas;
        AR_dict.getValue("pore_areas", pore_areas);
        for (int i = 0; i < pore_areas.size(); i++) {
            if (is3D) { pore_areas[i] *= scale*scale*scale; }
            else { pore_areas[i] *= scale*scale; }
        }
        AR_dict["pore_areas"] = pore_areas;
    }

    // Get the total areas from the mesh information, scale them, and reassign them
    {
        std::vector<double> total_areas;
        AR_dict.getValue("total_areas", total_areas);
        for (int i = 0; i < total_areas.size(); i++) {
            if (is3D) { total_areas[i] *= scale*scale*scale; }
            else { total_areas[i] *= scale*scale; }
        }
        AR_dict["total_areas"] = total_areas;
    }

    // Reassign the AR dictionary into the mesh info dictionary
    mesh_info["AR"] = &AR_dict;

    
    // Get the geometry dictionary from the mesh information
    JSONDict lengths_dict = *mesh_info["geometry"];

    // Get the large length scales from the mesh information, scale them, and reassign them
    {
        std::vector<double> L;
        lengths_dict.getValue("L", L);
        for (int i = 0; i < L.size(); i++) { L[i] *= scale; }
        lengths_dict["L"] = L;
    }

    // Get the small length scales from the mesh information, scale them, and reassign them
    {
        std::vector<double> l;
        lengths_dict.getValue("l", l);
        for (int i = 0; i < l.size(); i++) { l[i] *= scale; }
        lengths_dict["l"] = l;
    }

    // Reassign the geometry dictionary into the mesh info dictionary
    mesh_info["geometry"] = &lengths_dict;


    // Get the simulation_info dictionary from the mesh information
    JSONDict sim_info_dict = *mesh_info["simulation_info"];

    // Reassign the problem type
    sim_info_dict["problem_type"] = problem_type;

    // Reassign the simulation_info dictionary into the mesh info dictionary
    mesh_info["simulation_info"] = &sim_info_dict;


    // Save the mesh information
    mesh_info.saveToFile(scaled_mesh_info_file_path);
    cout << "    Complete!" << endl;
    

    // ========================================
    // Finalize
    // ========================================
    gmsh::finalize(); // Called to finish using the Gmsh C++ API


    return 0;
}
