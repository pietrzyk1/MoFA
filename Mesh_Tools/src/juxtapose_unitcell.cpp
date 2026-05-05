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

//                                      Juxtapose Unit-cell
// Description: This code is used to juxtapose unit-cells contiguously to create a
//              periodic macroscopic domain, as well as a macroscopic mesh. This is
//              used in the periodic implementation of MoFA models. Call this
//              function with an "upscaled" mesh (i.e., a mesh scaled to 1 in length)
//              that is geometrically periodic.
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

#include <gmsh.h>
#include <cassert>
#include "JSON_IO.h"
#include "gmsh_utils.h"





// ========================================
// Define global variables
// ========================================  
// Define a structure with global parameters
struct Params
{
    string FILENAME = "juxtapose_unitcell.cpp";
};
static Params globalVars; 





int main(int argc, char **argv)
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string config_path = "./";

    int is3D = 0; // If the domain is 3D, = 1. If 2D, = 0
    double L_x, L_y, L_z;
    std::vector<int> N_macroDomains = {1, 1, 1};
    
    double min_elem_size, max_elem_size;

    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";
    string unitcell_geo_dir = "./";
    string unitcell_geo_file_name = "geo.geo_unrolled";
    string output_mesh_dir = "./";
    string output_mesh_file_name = "mesh.msh";


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
        
        JSONDict sub_dict = *mesh_dict["info path"];
        sub_dict.getValue("directory", mesh_info_dir);
        sub_dict.getValue("file name", mesh_info_file_name);
        
        sub_dict = *mesh_dict["output geo path"];
        sub_dict.getValue("directory", unitcell_geo_dir);
        sub_dict.getValue("file name", unitcell_geo_file_name);

        sub_dict = *mesh_dict["details"];
        sub_dict.getValue("is3D", is3D);
        sub_dict.getValue("max element length", max_elem_size);
        sub_dict.getValue("min element length", min_elem_size);

        sub_dict = *mesh_dict["AR"];
        JSONDict subsub_dict = *sub_dict["L"];
        subsub_dict.getValue("x", L_x);
        subsub_dict.getValue("y", L_y);
        if (is3D) { subsub_dict.getValue("z", L_z); }

        sub_dict = *mesh_dict["juxtaposed unit-cells info"];
        sub_dict.getValue("N_macroDomains", N_macroDomains);
        sub_dict.getValue("directory", output_mesh_dir);
        sub_dict.getValue("file name", output_mesh_file_name);
    }
    else
    {
        cerr << globalVars.FILENAME << ": Error in loading config file." << endl;
        exit(1);
    }
    
    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string mesh_info_file_path = mesh_info_dir + mesh_info_file_name;
    string unitcell_geo_file_path = unitcell_geo_dir + unitcell_geo_file_name;
    string output_mesh_file_path = output_mesh_dir + output_mesh_file_name;
    
    // Make vectors L, ell, and N_AR (eventually, input should just give L, ell, and N_AR as vectors, not the components)
    std::vector<double> L = {L_x, L_y}; if (is3D == 1) { L.push_back( L_z ); }
    
    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    //OptionsParser args(argc, argv);
    //args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    //args.ParseCheck();


    // ========================================
    //   Initialize gmsh and open the unit-cell geometry file
    // ========================================
    // Gmsh must be initialized before using any functions in the C++ API
    gmsh::initialize();
    
    // Open the unit-cell geometry file
    cout << globalVars.FILENAME << ": Opening unit-cell geometry file..." << endl;
    gmsh::open(unitcell_geo_file_path);
    cout << "    Complete!" << endl;

    cout << globalVars.FILENAME << ": Loading mesh/geometry information..." << endl;
    // Load the mesh info file
    JSONDict mesh_info; mesh_info.loadFromFile(mesh_info_file_path);
    // Get the AR tags from the AR dictionary in the mesh information
    JSONDict AR_dict = *mesh_info["AR"];
    std::vector<int> AR_tags; AR_dict.getValue("tags", AR_tags);
    // Get the stokes physical groups from the mesh information
    JSONDict stokes_dict = *mesh_info["stokes"];
    std::vector<int> no_slip_PGs; stokes_dict.getValue("noslip1", no_slip_PGs);
    std::vector<int> stokes_top_PGs; stokes_dict.getValue("top1", stokes_top_PGs);
    std::vector<int> stokes_bottom_PGs; stokes_dict.getValue("bottom1", stokes_bottom_PGs);
    std::vector<int> stokes_outlet_PGs; stokes_dict.getValue("outlet1", stokes_outlet_PGs);
    cout << "    Complete!" << endl;


    // ========================================
    //   Translate base unit-cell geometry
    // ========================================
    cout << globalVars.FILENAME << ": Translate the base unit-cell geometry..." << endl;

    std::vector<double> initial_trans = {-(L[0]/2.0)*(N_macroDomains[0] - 1), -(L[1]/2.0)*(N_macroDomains[1] - 1), -(L[2]/2.0)*(N_macroDomains[2] - 1)};
    std::vector<std::pair<int, int>> dimtags_to_trans;
    for (int i_AR = 0; i_AR < AR_tags.size(); i_AR++) {
        // Get the gmsh entity tags of each AR physical group (due to the AR merging, there could be multiple gmsh surface/volume tags for each physical group)
        std::vector<int> entityTags; gmsh::model::getEntitiesForPhysicalGroup(2 + is3D, AR_tags[i_AR], entityTags);
    
        // Create a vector of dimTags, dimtags_to_copy, of the entities to be copied/translated
        for (int i = 0; i < entityTags.size(); i++) { dimtags_to_trans.push_back( {2 + is3D, entityTags[i]} ); }
    }

    // Translate the AR entities
    gmsh::model::geo::translate(dimtags_to_trans, initial_trans[0], initial_trans[1], initial_trans[2]);
    gmsh::model::geo::synchronize();

    cout << "    Complete!" << endl;
    


    int forPorescaleSim = 0;


    
    std::vector<int> macro_no_slip_PGs;
    /*
    // ========================================
    //   Copy and translate the no-slip tags for Stokes
    // ========================================
    cout << globalVars.FILENAME << ": Copy and translate the no-slip tags for Stokes Solver..." << endl;
    
    // Initialize the new no slip tags vectors
    macro_no_slip_PGs.insert(macro_no_slip_PGs.end(), no_slip_PGs.begin(), no_slip_PGs.end());
    
    // Copy and translate the unit-cell geometry. Save the macrostarts and new AR tags in the process
    for (int i_MZ = 0; i_MZ < N_macroDomains[2]; i_MZ++) {
        for (int i_MY = 0; i_MY < N_macroDomains[1]; i_MY++) {
            for (int i_MX = 0; i_MX < N_macroDomains[0]; i_MX++) {
                cout << globalVars.FILENAME << ": i_MX = " << i_MX << ", i_MY = " << i_MY << ", i_MZ = " << i_MZ << endl;
                
                // Skip the very first copy, as it already exists (i.e., it is the geometry that was loaded from the geometry file)
                if (i_MX == 0 && i_MY == 0 && i_MZ == 0) { continue; }
                
                // Declare a vector for the new AR PGs created from copying and translating the unit-cell geometry
                int new_no_slip_PG;

                // Copy and translate the unit-cell geometry
                CopyAndTranslatePGGeometry(no_slip_PGs, {i_MX*L[0], i_MY*L[1], i_MZ*L[2]}, new_no_slip_PG, 1 + is3D);
                
                // Save the macrostart and new AR tags
                //macro_no_slip_PGs.insert(macro_no_slip_PGs.end(), new_no_slip_PGs.begin(), new_no_slip_PGs.end());
                macro_no_slip_PGs.push_back( new_no_slip_PG );
            }
        }
    }

    cout << "    Complete!" << endl;


    // ========================================
    //   Copy and translate the top tags
    // ========================================
    cout << globalVars.FILENAME << ": Copy and translate the top tags..." << endl;
    
    // Initialize the new top tags vectors
    std::vector<int> macro_top_PGs;
    macro_top_PGs.insert(macro_top_PGs.end(), stokes_top_PGs.begin(), stokes_top_PGs.end());
    
    // Copy and translate the unit-cell geometry. Save the macrostarts and new AR tags in the process
    for (int i_MZ = 0; i_MZ < N_macroDomains[2]; i_MZ++) {
        for (int i_MY = 0; i_MY < N_macroDomains[1]; i_MY++) {
            for (int i_MX = 0; i_MX < N_macroDomains[0]; i_MX++) {
                cout << globalVars.FILENAME << ": i_MX = " << i_MX << ", i_MY = " << i_MY << ", i_MZ = " << i_MZ << endl;
                
                // Skip the very first copy, as it already exists (i.e., it is the geometry that was loaded from the geometry file)
                if (i_MX == 0 && i_MY == 0 && i_MZ == 0) { continue; }
                
                // Declare a vector for the new AR PGs created from copying and translating the unit-cell geometry
                int new_top_PG;

                // Copy and translate the unit-cell geometry
                CopyAndTranslatePGGeometry(stokes_top_PGs, {i_MX*L[0], i_MY*L[1], i_MZ*L[2]}, new_top_PG, 1 + is3D);
                
                // Save the macrostart and new AR tags
                macro_top_PGs.push_back( new_top_PG );
            }
        }
    }

    cout << "    Complete!" << endl;


    // ========================================
    //   Copy and translate the bottom tags
    // ========================================
    cout << globalVars.FILENAME << ": Copy and translate the bottom tags..." << endl;
    
    // Initialize the new top tags vectors
    std::vector<int> macro_bottom_PGs;
    macro_bottom_PGs.insert(macro_bottom_PGs.end(), stokes_bottom_PGs.begin(), stokes_bottom_PGs.end());
    
    // Copy and translate the unit-cell geometry. Save the macrostarts and new AR tags in the process
    for (int i_MZ = 0; i_MZ < N_macroDomains[2]; i_MZ++) {
        for (int i_MY = 0; i_MY < N_macroDomains[1]; i_MY++) {
            for (int i_MX = 0; i_MX < N_macroDomains[0]; i_MX++) {
                cout << globalVars.FILENAME << ": i_MX = " << i_MX << ", i_MY = " << i_MY << ", i_MZ = " << i_MZ << endl;
                
                // Skip the very first copy, as it already exists (i.e., it is the geometry that was loaded from the geometry file)
                if (i_MX == 0 && i_MY == 0 && i_MZ == 0) { continue; }
                
                // Declare a vector for the new AR PGs created from copying and translating the unit-cell geometry
                int new_bottom_PG;

                // Copy and translate the unit-cell geometry
                CopyAndTranslatePGGeometry(stokes_bottom_PGs, {i_MX*L[0], i_MY*L[1], i_MZ*L[2]}, new_bottom_PG, 1 + is3D);
                
                // Save the macrostart and new AR tags
                macro_bottom_PGs.push_back( new_bottom_PG );
            }
        }
    }

    cout << "    Complete!" << endl;


    // ========================================
    //   Copy and translate the outlet tags
    // ========================================
    cout << globalVars.FILENAME << ": Copy and translate the outlet tags..." << endl;
    
    // Initialize the new top tags vectors
    std::vector<int> macro_outlet_PGs;
    macro_outlet_PGs.insert(macro_outlet_PGs.end(), stokes_outlet_PGs.begin(), stokes_outlet_PGs.end());
    
    // Copy and translate the unit-cell geometry. Save the macrostarts and new AR tags in the process
    for (int i_MZ = 0; i_MZ < N_macroDomains[2]; i_MZ++) {
        for (int i_MY = 0; i_MY < N_macroDomains[1]; i_MY++) {
            for (int i_MX = 0; i_MX < N_macroDomains[0]; i_MX++) {
                cout << globalVars.FILENAME << ": i_MX = " << i_MX << ", i_MY = " << i_MY << ", i_MZ = " << i_MZ << endl;
                
                // Skip the very first copy, as it already exists (i.e., it is the geometry that was loaded from the geometry file)
                if (i_MX == 0 && i_MY == 0 && i_MZ == 0) { continue; }
                
                // Declare a vector for the new AR PGs created from copying and translating the unit-cell geometry
                int new_outlet_PG;

                // Copy and translate the unit-cell geometry
                CopyAndTranslatePGGeometry(stokes_outlet_PGs, {i_MX*L[0], i_MY*L[1], i_MZ*L[2]}, new_outlet_PG, 1 + is3D);
                
                // Save the macrostart and new AR tags
                macro_outlet_PGs.push_back( new_outlet_PG );
            }
        }
    }

    cout << "    Complete!" << endl;

    */
    // ========================================
    //   Copy and translate the unit-cell geometry
    // ========================================
    cout << globalVars.FILENAME << ": Copy and translate the unit-cell geometry..." << endl;
    
    // Initialize the macrostarts and new AR tags vectors
    std::vector<int> AR_tag_macrostarts, AR_tags_macro;
    AR_tag_macrostarts.push_back( AR_tags[0] );
    AR_tags_macro.insert(AR_tags_macro.end(), AR_tags.begin(), AR_tags.end());
    
    // Copy and translate the unit-cell geometry. Save the macrostarts and new AR tags in the process
    for (int i_MZ = 0; i_MZ < N_macroDomains[2]; i_MZ++) {
        for (int i_MY = 0; i_MY < N_macroDomains[1]; i_MY++) {
            for (int i_MX = 0; i_MX < N_macroDomains[0]; i_MX++) {
                cout << globalVars.FILENAME << ": i_MX = " << i_MX << ", i_MY = " << i_MY << ", i_MZ = " << i_MZ << endl;
                
                // Skip the very first copy, as it already exists (i.e., it is the geometry that was loaded from the geometry file)
                if (i_MX == 0 && i_MY == 0 && i_MZ == 0) { continue; }
                
                // Declare a vector for the new AR PGs created from copying and translating the unit-cell geometry
                std::vector<int> new_AR_PGs;

                // Copy and translate the unit-cell geometry
                CopyAndTranslatePGGeometry(AR_tags, {i_MX*L[0], i_MY*L[1], i_MZ*L[2]}, new_AR_PGs, 2 + is3D);
                
                // Save the macrostart and new AR tags
                AR_tag_macrostarts.push_back( new_AR_PGs[0] );
                AR_tags_macro.insert(AR_tags_macro.end(), new_AR_PGs.begin(), new_AR_PGs.end());
            }
        }
    }

    cout << "    Complete!" << endl;
    

    // ========================================
    // Generate a 2D mesh
    // ========================================
    cout << globalVars.FILENAME << ": Generating mesh..." << endl;
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 1.0); //min_elem_size);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 1.0); //max_elem_size);
    if (is3D) { gmsh::model::mesh::generate(3); } // The "3" is for 3D.
    else { gmsh::model::mesh::generate(2); } // The "2" is for 2D.
    cout << "    Complete!" << endl;
    

    // ========================================
    // Save the mesh
    // ========================================
    cout << globalVars.FILENAME << ": Saving mesh..." << endl;
    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2); // Make the saved version 2.2 so that it is MFEM compatible
    gmsh::write(output_mesh_file_path); // Save mesh with defined name
    cout << "    Complete!" << endl;
    

    // ========================================
    // Finalize
    // ========================================
    gmsh::finalize(); // Called to finish using the Gmsh C++ API
    




    // ========================================
    // Save the macro information in the dictionaries
    // ========================================
    cout << globalVars.FILENAME << ": Adding new AR tag info to mesh info file..." << endl;

    // Save the AR tag macrostarts
    AR_dict["macrodomain tags"] = AR_tags_macro;
    AR_dict["macrodomain tag starts"] = AR_tag_macrostarts;
    
    // Reassign the AR dictionary into the mesh info dictionary
    mesh_info["AR"] = &AR_dict;

    // Save the new PGs for the no slip condition and reassign the stokes dictionary into the mesh info dictionary
    if (forPorescaleSim) {
        stokes_dict["macro noslip1"] = macro_no_slip_PGs;
        mesh_info["stokes"] = &stokes_dict;
    }
    
    // Add a flag that tells other codes (i.e., scale_mesh.cpp) the mesh has been juxtaposed
    mesh_info["juxtaposed"] = 1;

    // Save the mesh information
    mesh_info.saveToFile(mesh_info_file_path);
    cout << "    Complete!" << endl;


    return 0;
}
