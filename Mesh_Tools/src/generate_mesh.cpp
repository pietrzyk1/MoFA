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
    }
    else
    {
        cerr << "generate_mesh.cpp: main: Error in loading config file." << endl;
        return 1;
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
    std::vector<std::pair<int, int>> tool_sfs;
    std::vector<std::pair<int, int>> tool_vols;
    std::vector<std::vector<std::vector<std::vector<double>>>> ARs_uc_geo;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> ARs_uc_geo_3D;

    if (is3D) { MakeARMesh_Uniform3DRectangular(N_AR_x, ell_x, L_x, N_AR_y, ell_y, L_y, N_AR_z, ell_z, L_z, ARs_uc, ARs_uc_geo_3D); }
    else { MakeARMesh_Uniform2DRectangular(N_AR_x, ell_x, L_x, N_AR_y, ell_y, L_y, ARs_uc, ARs_uc_geo); }
    

    // ========================================
    // Define the geometry to be cut into the averaging regions
    // ========================================
    std::vector<std::vector<std::vector<std::vector<double>>>> tool_geo; // This has indices [surface, lines for those surfaces, points in those lines, coords of those points]
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> tool_geo_3D; // This has indices [volumes, surfaces for those volumes, lines for those surfaces, points in those lines, coords of those points]
    
    std::vector<std::pair<int, int>> ARs;
    std::vector<std::vector<std::pair<int, int>>> mapOut;

    if (geometryToggle == 1)
    {
        if (is3D)
        {
            int cut_geo1 = gmsh::model::occ::addBox(-ell_x/4, -ell_y/4, -ell_z/4, ell_x/2, ell_y/2, ell_z/2);
            gmsh::model::occ::synchronize();
            tool_vols.push_back({3, cut_geo1});
            tool_geo_3D.push_back( getVolumeCoords(cut_geo1) );
            
            int cut_geo2 = gmsh::model::occ::addBox(ell_x*2.1, -ell_y/3, -ell_z * 2/3, ell_x*1.8, ell_y/3, ell_z/2);
            gmsh::model::occ::synchronize();
            tool_vols.push_back({3, cut_geo2});
            tool_geo_3D.push_back( getVolumeCoords(cut_geo2) );
            
            // Cut the averaging regions with the geometry
            gmsh::model::occ::cut(ARs_uc, tool_vols, ARs, mapOut);
        }
        else
        {
            //int cut_geo = createCircle(-ell_x/2, 0.0, 0.0, 0.3 * geo_scale, 0.05 * geo_scale);
            std::vector<int> cut_geo = getCutSurfaceFromFile(geometry_file_path, geo_scale);
            gmsh::model::occ::synchronize();
            for (int i_surf = 0; i_surf < cut_geo.size(); i_surf++)
            {
                tool_sfs.push_back({2, cut_geo[i_surf]});
                tool_geo.push_back( getSurfaceCoords(cut_geo[i_surf]) );
            }
            
            // Cut the averaging regions with the geometry
            gmsh::model::occ::cut(ARs_uc, tool_sfs, ARs, mapOut);
        }
    }
    else
    {
        if (is3D)
        {
            gmsh::model::occ::fragment(ARs_uc, tool_vols, ARs, mapOut);
        }
        else
        {
            gmsh::model::occ::fragment(ARs_uc, tool_sfs, ARs, mapOut);
        }
    }
    

    // ========================================
    // Synchronize the model
    // ========================================
    gmsh::model::occ::synchronize();
    
    
    // ========================================
    // Re-obtain the geometry curves
    // ========================================
    std::vector<std::vector<std::vector<double>>> AR_c_geo;
    std::vector<std::vector<std::pair<int, int>>> tool_sfs_lns(tool_geo.size());
    std::vector<std::vector<std::pair<int, int>>> ARs_lns(ARs.size());
    std::vector<std::pair<int, int>> AR_all_lns;
    bool isGeoLine;
    int sf_ind;
    
    std::vector<std::vector<std::vector<std::vector<double>>>> AR_c_geo_3D;
    std::vector<std::vector<std::pair<int, int>>> tool_vol_srfcs(tool_geo_3D.size());
    std::vector<std::vector<std::pair<int, int>>> ARs_srfcs(ARs.size());
    std::vector<std::pair<int, int>> AR_all_srfcs;
    bool isGeoSurface;
    int vol_ind;
    

    if (is3D)
    {
        // Get the boundaries of the geometry in each AR after the cut
        for (int i_ar = 0; i_ar < ARs.size(); i_ar++)
        {
            // Get the surface tags for the AR after the cut
            assert (ARs[i_ar].first == 3);
            gmsh::model::getBoundary({ARs[i_ar]}, AR_all_srfcs, false, false, false);
            
            // Get the geometry for the AR after the cut
            AR_c_geo_3D = getVolumeCoords(ARs[i_ar].second);

            // Look through and see if any of the surfaces match the previous geo surfaces
            for (int i_srfc_c = 0; i_srfc_c < AR_c_geo_3D.size(); i_srfc_c++)
            {
                isGeoSurface = false;
                for (int i_vol_uc = 0; i_vol_uc < tool_geo_3D.size(); i_vol_uc++) // Looking through the uncut geo volumes
                {
                    for (int i_srfc_uc = 0; i_srfc_uc < tool_geo_3D[i_vol_uc].size(); i_srfc_uc++) // Looking through the uncut geo surfaces
                    {
                        // See if the surfaces from the cut geometry and uncut geometry are parallel and share a corner point. If so, it is likely that they are/from the same surface
                        if (areParallel(AR_c_geo_3D[i_srfc_c], tool_geo_3D[i_vol_uc][i_srfc_uc]) && compareSurfaceEndPoints(AR_c_geo_3D[i_srfc_c], tool_geo_3D[i_vol_uc][i_srfc_uc]))
                        {
                            isGeoSurface = true;
                            vol_ind = i_vol_uc;
                            break;
                        }
                    }
                    if (isGeoSurface) { break; }
                }

                // If the cut geo surface did compare to the uncut geo surface, save it
                if (isGeoSurface)
                {
                    tool_vol_srfcs[vol_ind].push_back(AR_all_srfcs[i_srfc_c]); // NOTE: tool_vol_srfcs retains the same volumes as tool_geo. There may be more surfaces in each volume though
                }
                else
                {
                    ARs_srfcs[i_ar].push_back(AR_all_srfcs[i_srfc_c]); // Assuming that the remainder surfaces are the AR boundaries
                }
            }
        }
    }
    else
    {
        // Get the boundaries of the geometry in each AR after the cut
        for (int i_ar = 0; i_ar < ARs.size(); i_ar++)
        {
            // Get the line tags for the AR after the cut
            assert (ARs[i_ar].first == 2);
            gmsh::model::getBoundary({ARs[i_ar]}, AR_all_lns, false, false, false);
            
            // Get the geometry for the AR after the cut
            AR_c_geo = getSurfaceCoords(ARs[i_ar].second);

            // Look through and see if any of the lines match the previous geo lines
            for (int i_ln_c = 0; i_ln_c < AR_c_geo.size(); i_ln_c++)
            {
                isGeoLine = false;
                for (int i_sf_uc = 0; i_sf_uc < tool_geo.size(); i_sf_uc++) // Looking through the uncut geo surfaces
                {
                    for (int i_ln_uc = 0; i_ln_uc < tool_geo[i_sf_uc].size(); i_ln_uc++) // Looking through the uncut geo lines
                    {
                        // See if the lines from the cut geometry and uncut geometry are parallel and share an end point. If so, it is likely that they are/from the same line
                        if (areParallel(AR_c_geo[i_ln_c], tool_geo[i_sf_uc][i_ln_uc]) && compareLineEndPoints(AR_c_geo[i_ln_c], tool_geo[i_sf_uc][i_ln_uc]))
                        {
                            isGeoLine = true;
                            sf_ind = i_sf_uc;
                            break;
                        }
                    }
                    if (isGeoLine) { break; }
                }

                // If the cut geo line did compare to the uncut geo line, save it
                if (isGeoLine)
                {
                    tool_sfs_lns[sf_ind].push_back(AR_all_lns[i_ln_c]); // NOTE: tool_sfs_lns retains the same surfaces as tool_geo. There may be more lines in each surface though
                }
                else
                {
                    ARs_lns[i_ar].push_back(AR_all_lns[i_ln_c]); // Assuming that the remainder lines are the AR boundaries
                }
            }
        }
    }
    

    // ========================================
    // Tag the AR domains and decided if the
    // top, bottom, left, and right tags are
    // inlet, outlet, no-slip, or periodic
    // ========================================
    std::vector<int> AR_tags, inlet_ln_tags, outlet_ln_tags, no_slip_ln_tags, top_ln_tags, bottom_ln_tags, left_ln_tags, right_ln_tags, unused_ln_tags;
    int inlet_tag, outlet_tag, no_slip_tag, unused_tag;


    
    if (is3D)
    {

        for (int i_ar = 0; i_ar < ARs.size(); i_ar++)
        {
            int dim = ARs[i_ar].first;
            int tag = ARs[i_ar].second;
            
            assert (dim > 1);
            if (verbose) { std::cout << "occ::cut created entity: dim = " << dim << ", tag = " << tag << std::endl; }
            
            // Create a physical group from the AR
            AR_tags.push_back( gmsh::model::addPhysicalGroup(dim, {tag}) );
            
            for (int i_srfc = 0; i_srfc < ARs_srfcs[i_ar].size(); i_srfc++)
            {
                assert (ARs_srfcs[i_ar][i_srfc].first == 2); // Ensure that only surfaces have made it into ARs_srfcs
                
                // Check if the surface lies on the inlet (i.e., the left-most side of the domain).
                if (isCoincidesSurface(ARs_srfcs[i_ar][i_srfc].second, 0, -L_x/2))
                {
                    inlet_ln_tags.push_back( ARs_srfcs[i_ar][i_srfc].second );
                    if (verbose) std::cout << "Added surface " << ARs_srfcs[i_ar][i_srfc].second << " to the inlet Physical Group." << std::endl;
                }
                // Check if the surface lies on the outlet (i.e., the right-most side of the domain).
                else if (isCoincidesSurface(ARs_srfcs[i_ar][i_srfc].second, 0, L_x/2))
                {
                    outlet_ln_tags.push_back( ARs_srfcs[i_ar][i_srfc].second );
                    if (verbose) std::cout << "Added surface " << ARs_srfcs[i_ar][i_srfc].second << " to the outlet Physical Group." << std::endl;
                }
                // Check if the surface lies on the top boundary
                else if (isCoincidesSurface(ARs_srfcs[i_ar][i_srfc].second, 1, L_y/2))
                {
                    no_slip_ln_tags.push_back( ARs_srfcs[i_ar][i_srfc].second );
                    if (verbose) std::cout << "Added surface " << ARs_srfcs[i_ar][i_srfc].second << " to the no-slip Physical Group." << std::endl;
                }
                // Check if the surface lies on the bottom boundary
                else if (isCoincidesSurface(ARs_srfcs[i_ar][i_srfc].second, 1, -L_y/2))
                {
                    no_slip_ln_tags.push_back( ARs_srfcs[i_ar][i_srfc].second );
                    if (verbose) std::cout << "Added surface " << ARs_srfcs[i_ar][i_srfc].second << " to the no-slip Physical Group." << std::endl;
                }
                // Check if the surface lies on the front boundary
                else if (isCoincidesSurface(ARs_srfcs[i_ar][i_srfc].second, 2, L_z/2))
                {
                    no_slip_ln_tags.push_back( ARs_srfcs[i_ar][i_srfc].second );
                    if (verbose) std::cout << "Added surface " << ARs_srfcs[i_ar][i_srfc].second << " to the no-slip Physical Group." << std::endl;
                }
                // Check if the surface lies on the back boundary
                else if (isCoincidesSurface(ARs_srfcs[i_ar][i_srfc].second, 2, -L_z/2))
                {
                    no_slip_ln_tags.push_back( ARs_srfcs[i_ar][i_srfc].second );
                    if (verbose) std::cout << "Added surface " << ARs_srfcs[i_ar][i_srfc].second << " to the no-slip Physical Group." << std::endl;
                }
                // If none of the previous, put it in the unused tags
                else
                {
                    unused_ln_tags.push_back( ARs_srfcs[i_ar][i_srfc].second );
                    if (verbose) std::cout << "Added surface " << ARs_srfcs[i_ar][i_srfc].second << " to the unused Physical Group." << std::endl;
                }
            }
        }

        // ========================================
        // Add tool volumes to the no-slip tags physical group
        // ========================================
        for (int i_vol = 0; i_vol < tool_vol_srfcs.size(); i_vol++)
        {
            for (int i_srfc = 0; i_srfc < tool_vol_srfcs[i_vol].size(); i_srfc++)
            {
                no_slip_ln_tags.push_back(tool_vol_srfcs[i_vol][i_srfc].second);
            }
        }

        // ========================================
        // Create inlet, outlet, and no slip physical groups
        // ========================================
        inlet_tag = gmsh::model::addPhysicalGroup(2, inlet_ln_tags);
        outlet_tag = gmsh::model::addPhysicalGroup(2, outlet_ln_tags);
        no_slip_tag = gmsh::model::addPhysicalGroup(2, no_slip_ln_tags);
        unused_tag = gmsh::model::addPhysicalGroup(2, unused_ln_tags);

    }
    else
    {

        for (int i_ar = 0; i_ar < ARs.size(); i_ar++)
        {
            int dim = ARs[i_ar].first;
            int tag = ARs[i_ar].second;
            
            assert (dim > 1);
            if (verbose) { std::cout << "occ::cut created entity: dim = " << dim << ", tag = " << tag << std::endl; }
            
            // Create a physical group from the AR
            AR_tags.push_back( gmsh::model::addPhysicalGroup(dim, {tag}) );
            
            for (int i_ln = 0; i_ln < ARs_lns[i_ar].size(); i_ln++)
            {
                assert (ARs_lns[i_ar][i_ln].first == 1); // Ensure that only lines have made it into ARs_lns
                
                // Check if the line lies on the inlet (i.e., the left-most side of the domain).
                if (isCoincides(ARs_lns[i_ar][i_ln].second, 0, -L_x/2))
                {
                    left_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    //inlet_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    if (verbose) std::cout << "Added curve " << ARs_lns[i_ar][i_ln].second << " to the inlet Physical Group." << std::endl;
                }
                // Check if the line lies on the outlet (i.e., the right-most side of the domain).
                else if (isCoincides(ARs_lns[i_ar][i_ln].second, 0, L_x/2))
                {
                    right_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    //outlet_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    if (verbose) std::cout << "Added curve " << ARs_lns[i_ar][i_ln].second << " to the outlet Physical Group." << std::endl;
                }
                // Check if the line lies on the top boundary
                else if (isCoincides(ARs_lns[i_ar][i_ln].second, 1, L_y/2))
                {
                    top_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    //no_slip_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    if (verbose) std::cout << "Added curve " << ARs_lns[i_ar][i_ln].second << " to the no-slip Physical Group." << std::endl;
                }
                // Check if the line lies on the bottom boundary
                else if (isCoincides(ARs_lns[i_ar][i_ln].second, 1, -L_y/2))
                {
                    bottom_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    //no_slip_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    if (verbose) std::cout << "Added curve " << ARs_lns[i_ar][i_ln].second << " to the no-slip Physical Group." << std::endl;
                }
                // If none of the previous, put it in the unused tags
                else
                {
                    unused_ln_tags.push_back( ARs_lns[i_ar][i_ln].second );
                    if (verbose) std::cout << "Added curve " << ARs_lns[i_ar][i_ln].second << " to the unused Physical Group." << std::endl;
                }
            }
        }

        // ========================================
        // Add tool surfaces to the no-slip tags physical group
        // ========================================
        for (int i_srfc = 0; i_srfc < tool_sfs_lns.size(); i_srfc++)
        {
            for (int i_ln = 0; i_ln < tool_sfs_lns[i_srfc].size(); i_ln++)
            {
                no_slip_ln_tags.push_back(tool_sfs_lns[i_srfc][i_ln].second);
            }
        }

        // ========================================
        // Create inlet, outlet, and no slip physical groups
        // ========================================
        vector<int> used_tags; int used_tag;
        
        used_tag = getPhysicalGroup(inlet_tag, inletSide, 1, top_ln_tags, right_ln_tags, bottom_ln_tags, left_ln_tags);
        used_tags.push_back( used_tag );
        
        used_tag = getPhysicalGroup(outlet_tag, outletSide, 1, top_ln_tags, right_ln_tags, bottom_ln_tags, left_ln_tags);
        used_tags.push_back( used_tag );

        //if (isPeriodic[0] == 1)
        //{
        //    for (int i = 0; i < used_tags.size(); i++) { assert (used_tags[i] != 2 && used_tags[i] != 4); }
        //    used_tags.push_back( 2 ); used_tags.push_back( 4 );
        //}
        //if (isPeriodic[1] == 1)
        //{
        //    for (int i = 0; i < used_tags.size(); i++) { assert (used_tags[i] != 1 && used_tags[i] != 3); }
        //    used_tags.push_back( 1 ); used_tags.push_back( 3 );
        //}

        for (int i = 1; i <= 4; i++)
        {
            used_tag = 0;
            for (int j = 0; j < used_tags.size(); j++)
            {
                if (used_tags[j] == i) { used_tag = 1; }
            } 
            if (used_tag == 0)
            {
                if (i == 1) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), top_ln_tags.begin(), top_ln_tags.end()); }
                else if (i == 2) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), right_ln_tags.begin(), right_ln_tags.end()); }
                else if (i == 3) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), bottom_ln_tags.begin(), bottom_ln_tags.end()); }
                else if (i == 4) { no_slip_ln_tags.insert(no_slip_ln_tags.end(), left_ln_tags.begin(), left_ln_tags.end()); }
            }
        }
        
        no_slip_tag = gmsh::model::addPhysicalGroup(1, no_slip_ln_tags);
        unused_tag = gmsh::model::addPhysicalGroup(1, unused_ln_tags);

        /*
        if (false) // Classic, horizontal configuration
        {
            inlet_tag = gmsh::model::addPhysicalGroup(1, left_ln_tags);
            outlet_tag = gmsh::model::addPhysicalGroup(1, right_ln_tags);

            no_slip_ln_tags.insert(no_slip_ln_tags.end(), top_ln_tags.begin(), top_ln_tags.end());
            no_slip_ln_tags.insert(no_slip_ln_tags.end(), bottom_ln_tags.begin(), bottom_ln_tags.end());
            no_slip_tag = gmsh::model::addPhysicalGroup(1, no_slip_ln_tags);
            
            unused_tag = gmsh::model::addPhysicalGroup(1, unused_ln_tags);
        }
        else // Testing vertical configuration
        {
            inlet_tag = gmsh::model::addPhysicalGroup(1, bottom_ln_tags);
            outlet_tag = gmsh::model::addPhysicalGroup(1, top_ln_tags);

            no_slip_ln_tags.insert(no_slip_ln_tags.end(), left_ln_tags.begin(), left_ln_tags.end());
            no_slip_ln_tags.insert(no_slip_ln_tags.end(), right_ln_tags.begin(), right_ln_tags.end());
            no_slip_tag = gmsh::model::addPhysicalGroup(1, no_slip_ln_tags);
            
            unused_tag = gmsh::model::addPhysicalGroup(1, unused_ln_tags);
        }
        */
        //inlet_tag = gmsh::model::addPhysicalGroup(1, inlet_ln_tags);
        //outlet_tag = gmsh::model::addPhysicalGroup(1, outlet_ln_tags);
        //no_slip_tag = gmsh::model::addPhysicalGroup(1, no_slip_ln_tags);
        //unused_tag = gmsh::model::addPhysicalGroup(1, unused_ln_tags);

    }
    




    // ========================================
    // Prepare the mesh for periodic conditions
    // ========================================
    if (is3D)
    {
        if (isPeriodic[0] == 1 || isPeriodic[1] == 1 || isPeriodic[2] == 1)
        {
            cout << "CRITICAL ERROR: generate_mesh.cpp: main(): Periodic mesh generation for 3D is not implemented yet. Please make isPeriodic = {0, 0, 0}." << endl;
            return 1;
        }
    }
    else
    {
        if (isPeriodic[0] == 1)
        {
            for (int i = 0; i < left_ln_tags.size(); i++)
            {
                gmsh::model::geo::mesh::setTransfiniteCurve(left_ln_tags[i], (int)(ell_y/min_elem_size));
            }
            for (int i = 0; i < right_ln_tags.size(); i++)
            {
                gmsh::model::geo::mesh::setTransfiniteCurve(right_ln_tags[i], (int)(ell_y/min_elem_size));
            }
        }
        if (isPeriodic[1] == 1)
        {
            for (int i = 0; i < bottom_ln_tags.size(); i++)
            {
                gmsh::model::geo::mesh::setTransfiniteCurve(bottom_ln_tags[i], (int)(L_x/min_elem_size));
            }
            for (int i = 0; i < top_ln_tags.size(); i++)
            {
                gmsh::model::geo::mesh::setTransfiniteCurve(top_ln_tags[i], (int)(L_x/min_elem_size));
            }
        }
        gmsh::model::geo::synchronize();
    }
    

    
    
    
    // ========================================
    // Create and save the map structure for providing the physical groups to the solvers
    // ========================================
    JSONDict bdr_tag_file_new, AR_dict, stokes_dict, transport_dict, geometry, sim_info;
    
    // Create the AR sub-dictionary
    AR_dict["total_number"] = (int)AR_tags.size();
    AR_dict["tags"] = AR_tags; // The AR tags are not actually used yet. Maybe in the future.
    vector<double> AR_areas;
    if (is3D) { for (int i_ar = 0; i_ar < AR_tags.size(); i_ar++) { AR_areas.push_back( (double)(ell_x*ell_y*ell_z) ); } }
    else { for (int i_ar = 0; i_ar < AR_tags.size(); i_ar++) { AR_areas.push_back( (double)(ell_x*ell_y) ); } }
    AR_dict["total_areas"] = AR_areas;
    bdr_tag_file_new["AR"] = &AR_dict;
    
    // Create the Stokes sub-dictionary (containing tags for BCs in the Stokes problem)
    stokes_dict["inlet1"] = inlet_tag;
    stokes_dict["noslip1"] = no_slip_tag;
    bdr_tag_file_new["stokes"] = &stokes_dict;

    // Create the scalar closure sub-dictionary (containing tags for BCs in the scalar closure problem)
    transport_dict["inlet1"] = inlet_tag;
    transport_dict["outlet1"] = outlet_tag;
    bdr_tag_file_new["scalar_closure"] = &transport_dict;

    // Create the mesh geometry dictionary
    geometry["L"] = {L_x, L_y, L_z};
    geometry["l"] = {ell_x, ell_y, ell_z};
    bdr_tag_file_new["geometry"] = &geometry;

    // Create the simulation information dictionarys
    sim_info["problem_type"] = simulation_type;
    bdr_tag_file_new["simulation_info"] = &sim_info;

    bdr_tag_file_new.saveToFile(mesh_info_file_path);
    

    // ========================================
    // Generate a 2D mesh
    // ========================================
    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", min_elem_size);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", max_elem_size);
    if (is3D) { gmsh::model::mesh::generate(3); } // The "3" is for 3D.
    else { gmsh::model::mesh::generate(2); } // The "2" is for 2D.


    // ========================================
    // Save the mesh
    // ========================================
    gmsh::option::setNumber("Mesh.MshFileVersion", 2.2); // Make the saved version 2.2 so that it is MFEM compatible
    gmsh::write(output_file_path); // Save mesh with defined name


    // ========================================
    // Finalize
    // ========================================
    gmsh::finalize(); // Called to finish using the Gmsh C++ API
    

    return 0;
}
