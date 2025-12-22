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

//                                      Mesh Generator
// Description: This code and its headers are used to create a 2D or 3D mesh from a
//              rectangular domain. The domain can consider different topologies,
//              which are typically created to alternative means.
//
// Domain Schematic: (point and line numbers are not necessarily the same as those
//                   used in the code)
//
//                                      Line 4
//          Point 1 *------------------------------------------* Point 4
//                  |                  Length L                |
//           Line 1 | Width W                                  | Line 3
//                  |                                          |
//          Point 2 *------------------------------------------* Point 3
//                                      Line 2
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
//#include "mfem.hpp"
#include "JSON_IO.h"
#include "gmsh_utils.h"
#include "gmsh_import_STL.h"





// ========================================
// Define global variables
// ========================================  
// Declare whether warnings/action messages are output
constexpr bool verbose = false;





int main(int argc, char **argv)
{
    // ===============================================================
    //   Define variables (and their default values) that can be altered by the config file and command line options
    // ===============================================================
    string config_path = "./";

    string geometry_dir = "./";
    string geometry_file_name = "geo.txt";

    string gmsh_model_name = "default_mesh"; // Not really important
    int geometryToggle = 0; // Toggle whether to cut out the geometry; 0 is false (no cut), 1 is true (cut)
    string simulation_type = "porescale"; // Porescale or upscaled
    
    int is3D = 0; // If the domain is 3D, = 1. If 2D, = 0
    vector<int> isPeriodic = {0, 0, 0};
    string inletSide = "none";
    string outletSide = "none";

    int N_AR_x = 10; // Number of averaging regions in the x-direction
    double geo_scale = (simulation_type == "upscaled") ? 1. : 1./N_AR_x; // geo_scale = 1 is the upscaled model mesh, where AR are 1x1(x1)
    double L_x = N_AR_x  * geo_scale; // Domain length in x-direction
    double ell_x = L_x / N_AR_x; // Averaging region length in x-direction
    
    int N_AR_y = 1; // Number of averaging regions in the y-direction
    double L_y = N_AR_y * ell_x; // Domain length in y-direction
    double ell_y = L_y / N_AR_y; // Averaging region length in y-direction
    
    int N_AR_z = 0; // Number of averaging regions in the z-direction
    double L_z = N_AR_z * ell_x; // Domain length in z-direction
    double ell_z = ell_x; // Averaging region length in z-direction

    double max_elem_size = 0.25 * geo_scale;
    double min_elem_size = 0.025 * geo_scale;

    int merge_ARs = 0;
    double min_area_threshold = 0.2;
    double max_AR_length = 1.75;
    double max_AR_length_ratio = 2.0;

    std::vector<std::string> pg_names = {};
    std::vector<std::vector<int>> pg_cut_inds;
    
    string output_dir = "./";
    string output_file_name = "mesh.msh";
    string mesh_info_dir = "./";
    string mesh_info_file_name = "mesh_info.txt";


    // ===============================================================
    //   Search for config file path in argv (i.e., command line options) 
    // ===============================================================
    for (int i = 1; i < argc; i++)
    {
        if ((string(argv[i]) == "-C" || string(argv[i]) == "--config_path") && i + 1 < argc)
        {
            config_path = argv[i + 1];
            cout << "generate_mesh.cpp: Configuration path obtained from parser options: " << config_path << endl;
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
        
        JSONDict sub_dict = *mesh_dict["geometry path"];
        sub_dict.getValue("directory", geometry_dir);
        sub_dict.getValue("file name", geometry_file_name);
        sub_dict.getValue("scale", geo_scale);
        
        sub_dict = *mesh_dict["details"];
        sub_dict.getValue("gmsh model name", gmsh_model_name);
        sub_dict.getValue("geometry toggle", geometryToggle);
        sub_dict.getValue("simulation type", simulation_type);
        sub_dict.getValue("is3D", is3D);
        sub_dict.getValue("isPeriodic", isPeriodic);
        sub_dict.getValue("inlet side", inletSide);
        sub_dict.getValue("outlet side", outletSide);
        sub_dict.getValue("max element length", max_elem_size);
        sub_dict.getValue("min element length", min_elem_size);
        
        sub_dict = *mesh_dict["AR"];
        
        sub_dict.getValue("merge ARs", merge_ARs);
        sub_dict.getValue("min area threshold", min_area_threshold);
        sub_dict.getValue("max AR length", max_AR_length);
        sub_dict.getValue("max AR length ratio", max_AR_length_ratio);

        JSONDict subsub_dict = *sub_dict["N_AR"];
        subsub_dict.getValue("x", N_AR_x);
        subsub_dict.getValue("y", N_AR_y);
        if (is3D) { subsub_dict.getValue("z", N_AR_z); }
        
        subsub_dict = *sub_dict["L"];
        subsub_dict.getValue("x", L_x);
        subsub_dict.getValue("y", L_y);
        if (is3D) { subsub_dict.getValue("z", L_z); }

        subsub_dict = *sub_dict["ell"];
        subsub_dict.getValue("x", ell_x);
        subsub_dict.getValue("y", ell_y);
        if (is3D) { subsub_dict.getValue("z", ell_z); }

        sub_dict = *mesh_dict["output path"];
        sub_dict.getValue("directory", output_dir);
        sub_dict.getValue("file name", output_file_name);

        sub_dict = *mesh_dict["info path"];
        sub_dict.getValue("directory", mesh_info_dir);
        sub_dict.getValue("file name", mesh_info_file_name);

        sub_dict = *mesh_dict["cut physical groups"];
        sub_dict.getValue("physical group cut names", pg_names);
        sub_dict.getValue("physical group cut indices", pg_cut_inds);
    }
    else
    {
        cerr << "generate_mesh.cpp: main: Error in loading config file." << endl;
        exit(1);
    }
    
    // Initialize paths and file names from the parts taken from the defaults (either in the variable declarations or the config file)
    string geometry_file_path = geometry_dir + geometry_file_name;
    string output_file_path = output_dir + output_file_name;
    string mesh_info_file_path = mesh_info_dir + mesh_info_file_name;
    
    
    // ===============================================================
    //   Define the option parser and add options that can be changed from the command line
    // ===============================================================
    /*
    OptionsParser args(argc, argv);
    args.AddOption(&config_path, "-C", "--config_path", "The path to the configuration file.");
    args.ParseCheck();
    */


    // ========================================
    // Initialize the model
    // ========================================
    // Before using any functions in the C++ API, Gmsh must be initialized:
    gmsh::initialize();
    
    // Add a model
    gmsh::model::add(gmsh_model_name);
    

    // ========================================
    // Define averaging regions
    // ========================================
    std::vector<std::pair<int, int>> ARs_uc; // Format is {dimension of entity, integer tag of entity}.
    std::vector<std::vector<std::vector<std::vector<double>>>> ARs_uc_geo; // For 2D
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> ARs_uc_geo_3D; // For 3D
    std::vector<std::vector<int>> AR_neighbors_OLD; // Outer vector has N_AR slots; inner vector is a list of all neighboring AR numbers (Not really needed anymore? We get them later after AR merge)
    
    bool makeUniformARBackground = false;
    if (is3D) { if (makeUniformARBackground) { MakeARMesh_Uniform3DRectangular({N_AR_x, N_AR_y, N_AR_z}, {ell_x, ell_y, ell_z}, {L_x, L_y, L_z}, ARs_uc, ARs_uc_geo_3D); } }
    else { MakeARMesh_Uniform2DRectangular({N_AR_x, N_AR_y}, {ell_x, ell_y}, {L_x, L_y}, ARs_uc, ARs_uc_geo, AR_neighbors_OLD); } 
    
    
    // ========================================
    // Define the geometry to be cut into the averaging regions
    // ========================================
    std::vector<std::pair<int, int>> tool_sfs; // For 2D
    std::vector<std::vector<std::vector<std::vector<double>>>> tool_geo; // This has indices [surface, lines for those surfaces, points in those lines, coords of those points]
    std::vector<std::pair<int, int>> tool_vols; // For 3D
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> tool_geo_3D; // This has indices [volumes, surfaces for those volumes, lines for those surfaces, points in those lines, coords of those points]
    
    std::vector<std::pair<int, int>> ARs;
    std::vector<std::vector<std::pair<int, int>>> mapOut;

    if (geometryToggle == 1)
    {
        cout << "gmsh_utils.cpp: Cutting AR domain with geometry..." << endl;
        if (is3D)
        {
            // Old geometry-cutting tests. 
            /*
            // Create demo geometries
            int cut_geo1 = gmsh::model::occ::addBox(-ell_x/4, -ell_y/4, -ell_z/4, ell_x/2, ell_y/2, ell_z/2);
            gmsh::model::occ::synchronize();
            tool_vols.push_back({3, cut_geo1});
            tool_geo_3D.push_back( getVolumeCoords(cut_geo1) );
            int cut_geo2 = gmsh::model::occ::addBox(ell_x*2.1, -ell_y/3, -ell_z * 2/3, ell_x*1.8, ell_y/3, ell_z/2);
            gmsh::model::occ::synchronize();
            tool_vols.push_back({3, cut_geo2}); 
            tool_geo_3D.push_back( getVolumeCoords(cut_geo2) );
            
            // Cut the averaging regions with the geometry
            //gmsh::model::occ::cut(ARs_uc, tool_vols, ARs, mapOut)
            //gmsh::model::occ::intersect(ARs_uc, tool_vols, ARs, mapOut);
            */
            
            // TODO: allow for more than 10 boxes.
            string importAndCutMethod = "import and cut STL";
            if (importAndCutMethod == "import and cut STL") {
                std::vector<std::vector<double>> AR_planes;
                assert (N_AR_x == 10 && N_AR_y == 10 && N_AR_z == 10);
                if (simulation_type == "upscaled") {
                    AR_planes = {{-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
                                 {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0},
                                 {-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0}};
                }
                else if (simulation_type == "porescale") {
                    AR_planes = {{-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5},
                                 {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5},
                                 {-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5}};
                }
                else {
                    cerr << "generate_mesh.cpp: CRITICAL ERROR: Unknown string assigned to 'simulation_type'." << endl;
                    exit(1);
                }
                
                ARs = importAndCutSTLVolume(geometry_file_path, geo_scale, AR_planes);
                //gmsh::write("my_model.geo_unrolled"); // For debugging. Writes a (dense/heavy) file of the current gmsh model geometry
                //exit(1);
            }
            else {
                assert (makeUniformARBackground);
                // Imports geometry from an STL file (i.e., the triangles); makes the STL triangles gmsh surfaces;
                // makes a surface loop with the triangles; makes a gmsh volume with the surface loop. Then, uses
                // the gmsh::model::occ::cut() function (below) to cut the volume into AR volumes
                int cut_geo = createCutVolumeFromSTL(geometry_file_path, geo_scale);
                gmsh::model::occ::synchronize();
                tool_vols.push_back( {3, cut_geo} );
                //gmsh::write("my_model.geo_unrolled"); // For debugging. Writes a (dense/heavy) file of the current gmsh model geometry
                //exit(1);
                // Cut the averaging regions with the geometry
                gmsh::model::occ::cut(ARs_uc, tool_vols, ARs, mapOut);
                //gmsh::model::occ::intersect(ARs_uc, tool_vols, ARs, mapOut);
            }
        }
        else
        {
            // Load the cut geometry
            //load2DToolGeo(tool_sfs, tool_geo, geo_scale, -ell_x/2, 0.0, 0.0); // For a simple cut circle
            load2DToolGeo(tool_sfs, tool_geo, geo_scale, geometry_file_path);
            
            // Cut the averaging regions with the geometry
            gmsh::model::occ::cut(ARs_uc, tool_sfs, ARs, mapOut);
        }
        cout << "    Complete!" << endl;
    }
    else
    {
        if (is3D) { gmsh::model::occ::fragment(ARs_uc, tool_vols, ARs, mapOut); }
        else { gmsh::model::occ::fragment(ARs_uc, tool_sfs, ARs, mapOut); }
    }
    

    // ========================================
    // Synchronize the model
    // ========================================
    gmsh::model::occ::synchronize();
    
    
    // ========================================
    // Re-obtain the geometry curves
    // ========================================
    cout << "generate_mesh.cpp: Getting non-AR interface entity tags..." << endl;
    std::vector<std::vector<int>> entity_tags_in_AR;
    getNonARInterfaceEntityTags(ARs, entity_tags_in_AR, 2 + is3D);
    cout << "    Complete!" << endl;
    //std::vector<std::vector<std::pair<int, int>>> tool_sfs_lns(tool_geo.size()); // Re-obtained lines defining the surfaces of the cut geometry (2D)
    //std::vector<std::vector<std::pair<int, int>>> tool_vol_srfcs(tool_geo_3D.size()); // Re-obtained surfaces defining the volumes of the cut geometry (3D)
    //std::vector<std::vector<std::pair<int, int>>> ARs_lns(ARs.size()); // Re-obtained lines defining the AR boundaries (2D)
    //std::vector<std::vector<std::pair<int, int>>> ARs_srfcs(ARs.size()); // Re-obtained surfaces defining the AR boundaries (3D)
    //if (is3D) { separateARandGeometryVolumeSurfaces(ARs, tool_geo_3D, ARs_srfcs, tool_vol_srfcs); }
    //else { getNonARInterfaceEntityTags(ARs, ln_tags_in_AR, 2 + is3D); } //separateARandGeometrySurfaceLines(ARs, tool_geo, ARs_lns, tool_sfs_lns); }
    

    // ========================================
    // Users can identify certain tool entities
    // by the index in which they appear in tool_sfs,
    // tool_geo, tool_vols, and tool_geo_3D. Here, we
    // dig through the post-cut geometry and try to gain
    // the new line/surface tags that ID those tool entities,
    // so that they can be made into a physical group
    // ========================================
    std::vector<std::vector<int>> pg_ln_tags;
    
    if (pg_names.size() > 0) {
        cout << "generate_mesh.cpp: Defining boundary physical groups..." << endl;
        assert (pg_names.size() == pg_cut_inds.size());
        if (is3D) { cout << "generate_mesh.cpp: main(): NOTE: Defining boundary physical groups for 3D is not implemented yet. No boundary physical groups will be defined, except for the top, bottom, left, and right sides of the domain." << endl; }
        else {
            // Get the post-cut tags of the tool lines
            getToolGeoTags(pg_cut_inds, entity_tags_in_AR, tool_geo, {-L_x/2.0, L_x/2.0, -L_y/2.0, L_y/2.0}, pg_ln_tags);
        }
        cout << "    Complete!" << endl;
    }


    // ========================================
    // Tag the AR domains and decided if the
    // top, bottom, left, and right tags are
    // inlet, outlet, no-slip, or periodic
    // ========================================
    std::vector<int> AR_tags, no_slip_ln_tags, top_ln_tags, bottom_ln_tags, left_ln_tags, right_ln_tags, front_ln_tags, back_ln_tags, unused_ln_tags;
    int inlet_tag, outlet_tag, no_slip_tag, unused_tag;
    std::vector<int> inlet_AR_neighbors, outlet_AR_neighbors;
    bool anything_Merged = false;
    std::vector<double> AR_pore_space_areas;

    std::vector<std::vector<int>> domain_boundary_tags, boundary_neighbor_AR_tags, boundary_neighbor_PG;
    std::vector<std::vector<int>> AR_neighbors;
    
    
    if (is3D)
    {
        // Separate the tags of the surfaces at the domain boundaries (i.e., at +-L_x, +-L_y, and +-L_z) from the rest of the surface tags.
        // The remaining surface tags (in domain_boundary_tags[6]) are the cut geometry boundaries
        cout << "generate_mesh.cpp: Separating the domain boundary tags from cut surface/volume tags..." << endl;
        getDomainBoundarySurfaceTags(entity_tags_in_AR, domain_boundary_tags, {L_x, L_y, L_z});
        top_ln_tags = domain_boundary_tags[0];
        bottom_ln_tags = domain_boundary_tags[1];
        left_ln_tags = domain_boundary_tags[2];
        right_ln_tags = domain_boundary_tags[3];
        front_ln_tags = domain_boundary_tags[4];
        back_ln_tags = domain_boundary_tags[5];
        //unused_ln_tags = domain_boundary_tags[6];
        no_slip_ln_tags = domain_boundary_tags[6];

        // Add the top, bottom, front, and back domain boundaries to the no-slip condition
        no_slip_ln_tags.insert(no_slip_ln_tags.end(), top_ln_tags.begin(), top_ln_tags.end());
        no_slip_ln_tags.insert(no_slip_ln_tags.end(), bottom_ln_tags.begin(), bottom_ln_tags.end());
        no_slip_ln_tags.insert(no_slip_ln_tags.end(), front_ln_tags.begin(), front_ln_tags.end());
        no_slip_ln_tags.insert(no_slip_ln_tags.end(), back_ln_tags.begin(), back_ln_tags.end());
        cout << "    Complete!" << endl;

        
        // ========================================
        // Merge the small AR pore spaces with others (If things are merged, do not use AR average; use pore-space average in simulations) 
        // ========================================
        cout << "generate_mesh.cpp: Merging small AR with neighboring AR..." << endl;
        std::vector<std::vector<int>> AR_tags_toBeMerged; // The inner list is {merge host tag (i.e., the big AR), merging tag (i.e., the small AR)}
        if (merge_ARs == 1) { anything_Merged = createMergedARGroups(ARs, min_area_threshold, max_AR_length, max_AR_length_ratio, AR_tags_toBeMerged, AR_pore_space_areas); }
        
        if (anything_Merged) { cout << "    ARs were merged." << endl; } else { cout << "    No ARs were merged." << endl; }
        
        if (anything_Merged) {
            cout << "    Total number of AR tags: " << AR_tags_toBeMerged.size() << "." << endl;
            // Create a physical groups from the AR merge groups
            cout << "    Creating AR physical groups:" << endl;
            for (int i = 0; i < AR_tags_toBeMerged.size(); i++) {
                cout << "\r        Considering AR  i = " << i << "/" << AR_tags_toBeMerged.size() - 1 << "." << std::flush;
                AR_tags.push_back( gmsh::model::addPhysicalGroup(2 + is3D, AR_tags_toBeMerged[i]) );
            }
        }
        else {
            cout << "    Total number of AR tags: " << ARs.size() << "." << endl;
            // Clear AR_tags_toBeMerged (if merge_ARs == 1, but nothing was merged, AR_tags_toBeMerged will already be loaded, so clear it)
            AR_tags_toBeMerged.clear();
            // Create a physical group for each AR
            cout << "    Creating AR physical groups:" << endl;
            for (int i_ar = 0; i_ar < ARs.size(); i_ar++) {
                cout << "\r        Considering AR  i_ar = " << i_ar << "/" << ARs.size() - 1 << "." << std::flush;
                // Get AR information
                int dim = ARs[i_ar].first; assert (dim == 3);
                int tag = ARs[i_ar].second;
                AR_tags_toBeMerged.push_back( {tag} );
                AR_tags.push_back( gmsh::model::addPhysicalGroup(dim, AR_tags_toBeMerged[i_ar]) );
                
                // Compute AR area
                double AR_pore_space_area;
                gmsh::model::occ::getMass(dim, tag, AR_pore_space_area);
                AR_pore_space_areas.push_back( AR_pore_space_area );
            }
        }
        cout << "    Complete!" << endl;

        
        // ========================================
        // Get the neighboring ARs for each AR (this includes diagonal neighbors; the code looks for any AR that share a point/line/surface with other AR)
        // ========================================
        // Go through the AR_tags_toBeMerged, which provides the final number of AR (AR_tags_toBeMerged.size()), and the surface/volume tags to be merged to each that total of ARs
        cout << "generate_mesh.cpp: Obtaining lists of neighboring AR..." << endl;
        std::vector<std::vector<int>> AR_gpts;
        for (int i_ar = 0; i_ar < AR_tags_toBeMerged.size(); i_ar++) {
            cout << "\r    Getting data from AR  i_ar = " << i_ar << "/" << AR_tags_toBeMerged.size() - 1 << "." << std::flush;

            // Add a vector and create a reference for the AR gpts
            AR_gpts.push_back( std::vector<int>() );
            std::vector<int> &pts_vec = AR_gpts[i_ar];
            
            // Go through the tags in the AR group and collect the gmsh point tags
            for (int i_tag = 0; i_tag < AR_tags_toBeMerged[i_ar].size(); i_tag++) {
                std::vector<std::pair<int, int>> gpts;
                gmsh::model::getBoundary({{3, AR_tags_toBeMerged[i_ar][i_tag]}}, gpts, false, false, true);
                for (int i_pt = 0; i_pt < gpts.size(); i_pt++) { pts_vec.push_back( gpts[i_pt].second ); }
            }

            // Remove duplicates
            std::sort(pts_vec.begin(), pts_vec.end());
            auto last_dum = std::unique(pts_vec.begin(), pts_vec.end());
            pts_vec.erase(last_dum, pts_vec.end());
        }
        cout << endl;
        
        // Now, go through each set of points and see if there is any overlap between them. If there is, the ARs are neighbors
        for (int i_ar = 0; i_ar < AR_gpts.size(); i_ar++) {
            cout << "\r    Considering AR  i_ar = " << i_ar << "/" << AR_gpts.size() - 1 << "." << std::flush;

            // Add a vector for the AR
            AR_neighbors.push_back( std::vector<int>() );
            
            // Compare the gmsh points tags of i_ar against those of i_ar2
            for (int i_ar2 = 0; i_ar2 < AR_gpts.size(); i_ar2++) {
                // Skip if it is the same AR
                if (i_ar2 == i_ar) { continue; }
                // Combine the points tags of i_ar and i_ar2 into a single vector
                std::vector<int> pts_vec = AR_gpts[i_ar];
                std::vector<int> &pts_vec2 = AR_gpts[i_ar2];
                pts_vec.insert(pts_vec.end(), pts_vec2.begin(), pts_vec2.end());
                // Sort the vector and determine if there are any repeats with std::unique()
                std::sort(pts_vec.begin(), pts_vec.end());
                auto last_dum = std::unique(pts_vec.begin(), pts_vec.end());
                // If there are repeats (i.e., last_dum != pts_vec.end()), add the AR tag
                if (last_dum != pts_vec.end()) { AR_neighbors[i_ar].push_back( AR_tags[i_ar2] ); }
            }

            // Report if the AR does not have any AR neighbors
            if (AR_neighbors[i_ar].size() == 0) { cout << "generate_mesh.cpp: WARNING: AR " << i_ar << " does not have any neighboring AR! The code can continue, but it might be important to note this for your simulation." << endl; }
        }
        cout << endl << "    Complete!" << endl;
        


        /*
        // Old method of finding neighbors
        for (int i_ar = 0; i_ar < AR_tags_toBeMerged.size(); i_ar++) {
            cout << "\r    Considering AR  i_ar = " << i_ar << "/" << AR_tags_toBeMerged.size() << "." << std::flush;

            AR_neighbors.push_back( std::vector<int>() );
            // Go through the tags to be merged for a certain AR
            for (int i_tag = 0; i_tag < AR_tags_toBeMerged[i_ar].size(); i_tag++) {
                
                for (int i_ar2 = 0; i_ar2 < AR_tags_toBeMerged.size(); i_ar2++) {
                    if (i_ar2 == i_ar) { continue; }
                    for (int i_tag2 = 0; i_tag2 < AR_tags_toBeMerged[i_ar2].size(); i_tag2++) {
                        //if (areAdjacentEntities(2 + is3D, AR_tags_toBeMerged[i_ar][i_tag], AR_tags_toBeMerged[i_ar2][i_tag2])) { AR_neighbors[i_ar].push_back( AR_tags[i_ar2] ); break; }  // does not include diagonal neighbors
                        if (shareAnyBoundaryEntities(2 + is3D, AR_tags_toBeMerged[i_ar][i_tag], AR_tags_toBeMerged[i_ar2][i_tag2])) { AR_neighbors[i_ar].push_back( AR_tags[i_ar2] ); break; }  // includes diagonal neighbors
                    }
                }

            }

            // Remove the duplicates in the AR neighbor vector
            std::sort(AR_neighbors[i_ar].begin(), AR_neighbors[i_ar].end());
            auto last_unique = std::unique(AR_neighbors[i_ar].begin(), AR_neighbors[i_ar].end());
            AR_neighbors[i_ar].erase(last_unique, AR_neighbors[i_ar].end());

            // Report if the AR does not have any AR neighbors
            if (AR_neighbors[i_ar].size() == 0) { cout << "generate_mesh.cpp: WARNING: AR " << i_ar << " does not have any neighboring AR! The code can continue, but it might be important to note this for your simulation." << endl; }
        }
        cout << endl << "    Complete!" << endl;
        */


        // ========================================
        // Create inlet, outlet, and no slip physical groups
        // ========================================
        inlet_tag = gmsh::model::addPhysicalGroup(2, left_ln_tags);
        outlet_tag = gmsh::model::addPhysicalGroup(2, right_ln_tags);
        no_slip_tag = gmsh::model::addPhysicalGroup(2, no_slip_ln_tags);
        //unused_tag = gmsh::model::addPhysicalGroup(2, unused_ln_tags);
        
    }
    else
    {
        // Make sure is3D is 0 (i.e., we are looking at a 2D domain for this)
        assert (is3D == 0);

        // Separate the tags of the lines at the domain boundaries (i.e., at +-L_x and +-L_y) from the rest of the line tags.
        // The remaining line tags (in domain_boundary_tags[4]) are the cut geometry boundaries
        cout << "generate_mesh.cpp: Separating the domain boundary tags from cut surface/volume tags..." << endl;
        getDomainBoundaryLineTags(entity_tags_in_AR, domain_boundary_tags, ARs, boundary_neighbor_AR_tags, L_x, L_y, min_elem_size * 0.01); // tolerance can be changed here if needed.
        top_ln_tags = domain_boundary_tags[0];
        right_ln_tags = domain_boundary_tags[1];
        bottom_ln_tags = domain_boundary_tags[2];
        left_ln_tags = domain_boundary_tags[3];
        //unused_ln_tags = domain_boundary_tags[4];
        no_slip_ln_tags = domain_boundary_tags[4];

        
        // Remove duplicates from boundary_neighbor_AR_tags and initialize boundary_neighbor_PG
        for (int i = 0; i < boundary_neighbor_AR_tags.size(); i++) {
            std::sort(boundary_neighbor_AR_tags[i].begin(), boundary_neighbor_AR_tags[i].end());
            auto last = std::unique(boundary_neighbor_AR_tags[i].begin(), boundary_neighbor_AR_tags[i].end());
            boundary_neighbor_AR_tags[i].erase(last, boundary_neighbor_AR_tags[i].end());
            boundary_neighbor_PG.push_back( std::vector<int>() );
        }
        cout << "    Complete!" << endl;

        

        

        // ========================================
        // Merge the small AR pore spaces with others (If things are merged, do not use AR average; use pore-space average in simulations) 
        // ========================================
        cout << "generate_mesh.cpp: Merging small AR with neighboring AR..." << endl;
        std::vector<std::vector<int>> AR_tags_toBeMerged; // The inner list is {merge host tag (i.e., the big AR), merging tag (i.e., the smoll AR)}
        if (merge_ARs == 1) { anything_Merged = createMergedARGroups(ARs, min_area_threshold, max_AR_length, max_AR_length_ratio, AR_tags_toBeMerged, AR_pore_space_areas); }
        if (anything_Merged) { cout << "    ARs were merged." << endl; } else { cout << "    No ARs were merged." << endl; }
        
        if (anything_Merged)
        {
            cout << "    Total number of AR tags: " << AR_tags_toBeMerged.size() << "." << endl;
            // Create a physical groups from the AR merge groups
            cout << "    Creating AR physical groups:" << endl;
            for (int i = 0; i < AR_tags_toBeMerged.size(); i++) {
                cout << "\r        Considering AR  i = " << i << "/" << AR_tags_toBeMerged.size() - 1 << "." << std::flush;
                AR_tags.push_back( gmsh::model::addPhysicalGroup(2 + is3D, AR_tags_toBeMerged[i]) );

                // Collect the physical group tag for the (merged) AR in boundary_neighbor_PG for the (merged) AR tags that touch the boundaries
                for (int i_bn = 0; i_bn < boundary_neighbor_AR_tags.size(); i_bn++) {
                    std::vector<int> combined_vec = boundary_neighbor_AR_tags[i_bn];
                    combined_vec.insert(combined_vec.end(),  AR_tags_toBeMerged[i].begin(), AR_tags_toBeMerged[i].end());
                    std::sort(combined_vec.begin(), combined_vec.end());
                    auto last = std::unique(combined_vec.begin(), combined_vec.end());
                    if (last != combined_vec.end()) { boundary_neighbor_PG[i_bn].push_back( AR_tags[AR_tags.size() - 1] ); }
                }
            }
        }
        else
        {
            cout << "    Total number of AR tags: " << ARs.size() << "." << endl;
            // Clear AR_tags_toBeMerged (if merge_ARs == 1, but nothing was merged, AR_tags_toBeMerged will already be loaded, so clear it)
            AR_tags_toBeMerged.clear();
            AR_pore_space_areas.clear();
            // Create a physical group for each AR
            cout << "    Creating AR physical groups:" << endl;
            for (int i_ar = 0; i_ar < ARs.size(); i_ar++) {
                cout << "\r        Considering AR  i_ar = " << i_ar << "/" << ARs.size() - 1 << "." << std::flush;
                // Get AR information
                int dim = ARs[i_ar].first; int tag = ARs[i_ar].second;
                // Save the AR tags as "tobeMerged" arrays
                AR_tags_toBeMerged.push_back( {tag} );
                // Create the physical groups
                AR_tags.push_back( gmsh::model::addPhysicalGroup(dim, AR_tags_toBeMerged[i_ar]) );
                
                // Compute AR area
                double AR_pore_space_area;
                gmsh::model::occ::getMass(dim, tag, AR_pore_space_area);
                AR_pore_space_areas.push_back( AR_pore_space_area );
                
                // Collect the physical group tag for the (merged) AR in boundary_neighbor_PG for the (merged) AR tags that touch the boundaries
                for (int i_bn = 0; i_bn < boundary_neighbor_AR_tags.size(); i_bn++) {
                    std::vector<int> combined_vec = boundary_neighbor_AR_tags[i_bn];
                    combined_vec.push_back( tag );
                    std::sort(combined_vec.begin(), combined_vec.end());
                    auto last = std::unique(combined_vec.begin(), combined_vec.end());
                    if (last != combined_vec.end()) { boundary_neighbor_PG[i_bn].push_back( AR_tags[AR_tags.size() - 1] ); }
                    //for (int i_tag = 0; i_tag < boundary_neighbor_AR_tags[i_bn].size(); i_tag++) {
                    //    if (boundary_neighbor_AR_tags[i_bn][i_tag] == tag) { boundary_neighbor_PG[i_bn].push_back( AR_tags[AR_tags.size() - 1] ); break; }
                    //}
                }
            }
        }
        cout << "    Complete!" << endl;


        // ========================================
        // Get the neighboring ARs for each AR (this includes diagonal neighbors; the code looks for any AR that share a point/line/surface with other AR)
        // ========================================
        // Go through the AR_tags_toBeMerged, which provides the final number of AR (AR_tags_toBeMerged.size()), and the surface/volume tags to be merged to each that total of ARs
        cout << "generate_mesh.cpp: Obtaining neighboring ARs lists..." << endl;
        for (int i_ar = 0; i_ar < AR_tags_toBeMerged.size(); i_ar++) {
            cout << "\r    Considering AR  i_ar = " << i_ar << "/" << AR_tags_toBeMerged.size() - 1 << "." << std::flush;
            AR_neighbors.push_back( std::vector<int>() );
            // Go through the tags to be merged for a certain AR
            for (int i_tag = 0; i_tag < AR_tags_toBeMerged[i_ar].size(); i_tag++) {
                
                for (int i_ar2 = 0; i_ar2 < AR_tags_toBeMerged.size(); i_ar2++) {
                    if (i_ar2 == i_ar) { continue; }
                    for (int i_tag2 = 0; i_tag2 < AR_tags_toBeMerged[i_ar2].size(); i_tag2++) {
                        //if (areAdjacentEntities(2 + is3D, AR_tags_toBeMerged[i_ar][i_tag], AR_tags_toBeMerged[i_ar2][i_tag2])) { AR_neighbors[i_ar].push_back( AR_tags[i_ar2] ); break; } // does not include diagonal neighbors
                        if (shareAnyBoundaryEntities(2 + is3D, AR_tags_toBeMerged[i_ar][i_tag], AR_tags_toBeMerged[i_ar2][i_tag2])) { AR_neighbors[i_ar].push_back( AR_tags[i_ar2] ); break; } // includes diagonal neighbors
                    }
                }

            }

            // Remove the duplicates in the AR neighbor vector
            std::sort(AR_neighbors[i_ar].begin(), AR_neighbors[i_ar].end());
            auto last_unique = std::unique(AR_neighbors[i_ar].begin(), AR_neighbors[i_ar].end());
            AR_neighbors[i_ar].erase(last_unique, AR_neighbors[i_ar].end());

            // Report if the AR does not have any AR neighbors
            if (AR_neighbors[i_ar].size() == 0) { cout << "generate_mesh.cpp: WARNING: AR " << i_ar << " does not have any neighboring AR! The code can continue, but it might be important to note this for your simulation." << endl; }
        }
        cout << "    Complete!" << endl;
        


        // Moved these below the "Merge the small AR pore spaces" for now to allow AR tags start at 1.
        
        // ========================================
        // Create inlet, outlet, and no slip physical groups
        // ========================================
        vector<int> used_tags; int used_tag;
        
        used_tag = getPhysicalGroup(inlet_tag, inletSide, 1, top_ln_tags, right_ln_tags, bottom_ln_tags, left_ln_tags);
        used_tags.push_back( used_tag );
        inlet_AR_neighbors = boundary_neighbor_PG[used_tag - 1];

        used_tag = getPhysicalGroup(outlet_tag, outletSide, 1, top_ln_tags, right_ln_tags, bottom_ln_tags, left_ln_tags);
        used_tags.push_back( used_tag );
        outlet_AR_neighbors = boundary_neighbor_PG[used_tag - 1];

        // Go through the 4 options (top, right, bottom, left) for what could be used as the inlet/outlet
        for (int i = 1; i <= 4; i++)
        {
            // See if the option was used
            used_tag = 0;
            for (int j = 0; j < used_tags.size(); j++) { if (used_tags[j] == i) { used_tag = 1; } }
            
            // If not, add it to the no slip tags
            if (used_tag == 0)
            {
                if (i == 1) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), top_ln_tags.begin(), top_ln_tags.end()); }
                else if (i == 2) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), right_ln_tags.begin(), right_ln_tags.end()); }
                else if (i == 3) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), bottom_ln_tags.begin(), bottom_ln_tags.end()); }
                else if (i == 4) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), left_ln_tags.begin(), left_ln_tags.end()); }
            }
        }
        
        // Define the no-slip tags
        no_slip_tag = gmsh::model::addPhysicalGroup(1, no_slip_ln_tags);
        //unused_tag = gmsh::model::addPhysicalGroup(1, unused_ln_tags);
    }
    


    // ========================================
    // Prepare the mesh for periodic conditions
    // ========================================
    cout << "generate_mesh.cpp: Preparing mesh for periodic boundary conditions (if needed/in input file)..." << endl;
    if (is3D) {
        if (isPeriodic[0] == 1) {
            for (int i = 0; i < left_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(left_ln_tags[i], (int)(L_y*L_z/min_elem_size)); }
            for (int i = 0; i < right_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(right_ln_tags[i], (int)(L_y*L_z/min_elem_size)); }
        }
        if (isPeriodic[1] == 1) {
            for (int i = 0; i < bottom_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(bottom_ln_tags[i], (int)(L_x*L_z/min_elem_size)); }
            for (int i = 0; i < top_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(top_ln_tags[i], (int)(L_x*L_z/min_elem_size)); }
        }
        if (isPeriodic[2] == 1) {
            for (int i = 0; i < back_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(back_ln_tags[i], (int)(L_x*L_y/min_elem_size)); }
            for (int i = 0; i < front_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(front_ln_tags[i], (int)(L_x*L_y/min_elem_size)); }
        }
        gmsh::model::geo::synchronize();
    }
    else {
        if (isPeriodic[0] == 1) {
            for (int i = 0; i < left_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(left_ln_tags[i], (int)(ell_y/min_elem_size)); } // update: this should probably be the (length of each left_ln_tags) / min_elem_size
            for (int i = 0; i < right_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(right_ln_tags[i], (int)(ell_y/min_elem_size)); } // update: this should probably be the (length of each left_ln_tags) / min_elem_size
        }
        if (isPeriodic[1] == 1) {
            for (int i = 0; i < bottom_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(bottom_ln_tags[i], (int)(L_x/min_elem_size)); }
            for (int i = 0; i < top_ln_tags.size(); i++) { gmsh::model::geo::mesh::setTransfiniteCurve(top_ln_tags[i], (int)(L_x/min_elem_size)); }
        }
        gmsh::model::geo::synchronize();
    }
    cout << "    Complete!" << endl;
    

    
    
    // ========================================
    // Create and save the map structure for providing the physical groups to the solvers
    // ========================================
    cout << "generate_mesh.cpp: Preparing mesh info file..." << endl;
    JSONDict bdr_tag_file_new, AR_dict, stokes_dict, transport_dict, geometry, sim_info;
    
    // Create the AR sub-dictionary
    AR_dict["total_number"] = (int)AR_tags.size();
    AR_dict["tags"] = AR_tags;
    if (anything_Merged) { AR_dict["total_areas"] = AR_pore_space_areas; AR_dict["area_type"] = "pore"; }
    else {
        vector<double> AR_areas;
        if (is3D) { for (int i_ar = 0; i_ar < AR_tags.size(); i_ar++) { AR_areas.push_back( (double)(ell_x*ell_y*ell_z) ); } }
        else { for (int i_ar = 0; i_ar < AR_tags.size(); i_ar++) { AR_areas.push_back( (double)(ell_x*ell_y) ); } }
        AR_dict["total_areas"] = AR_areas;
        AR_dict["area_type"] = "AR";
    }
    AR_dict["pore_areas"] = AR_pore_space_areas; // We now compute and record the pore-space areas in this file, as opposed to the solver file
    AR_dict["neighbors"] = AR_neighbors;
    bdr_tag_file_new["AR"] = &AR_dict;
    
    // Create the Stokes sub-dictionary (containing tags for BCs in the Stokes problem)
    stokes_dict["inlet1"] = inlet_tag;
    stokes_dict["inlet1 AR neighbors"] = inlet_AR_neighbors;
    stokes_dict["noslip1"] = no_slip_tag;
    bdr_tag_file_new["stokes"] = &stokes_dict;

    // Create the scalar closure sub-dictionary (containing tags for BCs in the scalar closure problem)
    transport_dict["inlet1"] = inlet_tag;
    transport_dict["inlet1 AR neighbors"] = inlet_AR_neighbors;
    transport_dict["outlet1"] = outlet_tag;
    transport_dict["outlet1 AR neighbors"] = outlet_AR_neighbors;
    // Create the boundary physical groups defined by the user
    if (!is3D && pg_names.size() != 0) {
        // Add the physical groups, where the belonging lines tags were previous obtained. We add them here so that they can be observed in gmsh.exe (i.e., they overlay the no-slip tags; the lines are tagged as belonging to both physics groups, but in gmsh.exe, this allows for viewing of the user defined boundary physical groups)
        std::vector<int> pg_cut_inds_pgs;
        for (int i_pg = 0; i_pg < pg_ln_tags.size(); i_pg++) { pg_cut_inds_pgs.push_back( gmsh::model::addPhysicalGroup(1, pg_ln_tags[i_pg]) ); }
        for (int i_pg = 0; i_pg < pg_cut_inds_pgs.size(); i_pg++) { transport_dict[pg_names[i_pg]] = pg_cut_inds_pgs[i_pg]; }
    }
    bdr_tag_file_new["scalar_closure"] = &transport_dict;

    // Create the mesh geometry dictionary
    geometry["L"] = {L_x, L_y, L_z};
    geometry["l"] = {ell_x, ell_y, ell_z};
    bdr_tag_file_new["geometry"] = &geometry;

    // Create the simulation information dictionarys
    sim_info["problem_type"] = simulation_type;
    bdr_tag_file_new["simulation_info"] = &sim_info;

    bdr_tag_file_new.saveToFile(mesh_info_file_path);
    cout << "    Complete!" << endl;
    

    // ========================================
    // Generate a 2D mesh
    // ========================================
    cout << "generate_mesh.cpp: Generating and saving mesh..." << endl;
    //generateMesh(output_file_path, is3D, {min_elem_size, max_elem_size});
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", min_elem_size);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", max_elem_size);
    if (is3D) { gmsh::model::mesh::generate(3); } // The "3" is for 3D.
    else { gmsh::model::mesh::generate(2); } // The "2" is for 2D.


    // ========================================
    // Save the mesh
    // ========================================
    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2); // Make the saved version 2.2 so that it is MFEM compatible
    gmsh::write(output_file_path); // Save mesh with defined name
    cout << "    Complete!" << endl;
    

    // ========================================
    // Finalize
    // ========================================
    gmsh::finalize(); // Called to finish using the Gmsh C++ API
    

    return 0;
}
