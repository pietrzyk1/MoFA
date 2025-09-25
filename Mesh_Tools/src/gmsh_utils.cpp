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

#include <vector>
#include <string>
#include <gmsh.h>
#include <cassert>
#include "gmsh_utils.h"



// ========================================
// Functions for getting point coordinates of different geometric entities
// ========================================  

// Define a function that takes in a point tag and returns the point coordinates.
std::vector<double> getPointCoords(int pt_tag)
{
    std::vector<double> pt_coords;
    gmsh::model::getValue(0, pt_tag, {0.0}, pt_coords); // Get the coordinates of the points.
    return pt_coords;
}

// Define a function that takes in a line tag and returns the coordinates of the end points.
std::vector<std::vector<double>> getLineCoords(int ln_tag)
{
    std::vector<std::vector<double>> ln_pts_coords;
    std::vector<std::pair<int, int>> pts_dimtag;
    
    gmsh::model::getBoundary({{1, ln_tag}}, pts_dimtag, false, false, false);

    for (int i_pts = 0; i_pts < pts_dimtag.size(); i_pts++)
    {
        ln_pts_coords.push_back(getPointCoords(pts_dimtag[i_pts].second));
    }

    return ln_pts_coords;
}

// Define a function that takes in a surface tag and returns the coordinates of the end points of the boundary lines.
std::vector<std::vector<std::vector<double>>> getSurfaceCoords(int srfc_tag)
{
    std::vector<std::vector<std::vector<double>>> srfc_ln_pts_coords;
    std::vector<std::pair<int, int>> lns_dimtag;

    gmsh::model::getBoundary({{2, srfc_tag}}, lns_dimtag, false, false, false);          

    for (int i_ln = 0; i_ln < lns_dimtag.size(); i_ln++)
    {
        srfc_ln_pts_coords.push_back(getLineCoords(lns_dimtag[i_ln].second));
    }

    return srfc_ln_pts_coords;
}

// Define a function that takes in a volume tag and returns the coordinates of the points creating the boundary surfaces.
std::vector<std::vector<std::vector<std::vector<double>>>> getVolumeCoords(int vol_tag)
{
    std::vector<std::vector<std::vector<std::vector<double>>>> vol_srfc_ln_pts_coords;
    std::vector<std::pair<int, int>> srfc_dimtag;

    gmsh::model::getBoundary({{3, vol_tag}}, srfc_dimtag, false, false, false);          

    for (int i_srfc = 0; i_srfc < srfc_dimtag.size(); i_srfc++)
    {
        vol_srfc_ln_pts_coords.push_back(getSurfaceCoords(srfc_dimtag[i_srfc].second));
    }

    return vol_srfc_ln_pts_coords;
}





// ========================================
// Functions for mathematical calculations
// ========================================  

// Compute the cross product of two 3D vectors
std::vector<double> cross(const std::vector<double>& a, const std::vector<double> &b)
{
    // Make sure each vector has 3 components
    assert (a.size() == 3);
    assert (b.size() == 3);

    std::vector<double> sol;
    sol.push_back( a[1]*b[2] - a[2]*b[1] );
    sol.push_back( a[2]*b[0] - a[0]*b[2] );
    sol.push_back( a[0]*b[1] - a[1]*b[0] );
    return sol;
}

// Compute the magnitude/norm of a 3D vector
double norm(const std::vector<double>& a)
{
    // Make sure the vector has 3 components
    assert (a.size() == 3);
    return std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

// Define a function that computes the direction vector between 2 points
std::vector<double> computeDirectionVector(const std::vector<double> &a, const std::vector<double> &b)
{
    // Make sure there are 3 coordinates to each point provided
    assert (a.size() == 3);
    assert (b.size() == 3);

    std::vector<double> direction_vec;
    for (int i = 0; i < a.size(); i++) { direction_vec.push_back(a[i] - b[i]); }
    return direction_vec;
}

// Define a function that computes the direction vectors for a flat 3 point surface, and then return the normal
std::vector<double> computeSurfaceNormal_3Point(const std::vector<std::vector<double>> &a)
{
    assert (a.size() == 3);

    std::vector<double> direction_vec1 = computeDirectionVector(a[1], a[0]);
    std::vector<double> direction_vec2 = computeDirectionVector(a[2], a[0]);

    return cross(direction_vec1, direction_vec2);
}





// ========================================
// Functions for comparing different geometric entities
// ========================================  

// Define a function for comparing two doubles to within a given tolerance
bool compareValues(const double a, const double b, const double tol)
{
    return std::abs(a - b) < tol;
}

// Define a function for comparing two vectors of doubles. The components are compared to within a given tolerance
bool compareVectors(const std::vector<double> a, const std::vector<double> b, const double tol)
{
    assert (a.size() == b.size());
    bool sol = true;
    for (int i = 0; i < a.size(); i++)
    {
        if (!compareValues(a[i], b[i], tol))
        {
            sol = false;
            break;
        }
    }
    return sol;
}

// Define a function for comparing two vectors of doubles. The components are compared to within a given tolerance
bool compareLines(const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b, const double tol)
{
    // Make sure each line has 2 points
    assert (a.size() == 2);
    assert (b.size() == 2);
    bool sol = false;
    if (compareVectors(a[0], b[0], tol) && compareVectors(a[1], b[1], tol)) { sol = true; }
    else if (compareVectors(a[0], b[1], tol) && compareVectors(a[1], b[0], tol)) { sol = true; }
    return sol;
}

// Define a function for getting three corner points of a flat surface (to eventually evaluate the surface's normal vector)
std::vector<std::vector<double>> getThreeSurfacePoints(const std::vector<std::vector<std::vector<double>>> &a)
{
    // Make sure the surface has at least 3 lines
    assert (a.size() >= 3);
    
    std::vector<std::vector<double>> a_points;
    a_points.push_back( a[0][0] ); // Use the points of the first line, and one of the points from the second line. 
    a_points.push_back( a[0][1] );
    if ( compareVectors(a[0][0], a[1][0]) || compareVectors(a[0][1], a[1][0]) ) { a_points.push_back( a[1][1] ); }
    else { a_points.push_back( a[1][0] ); }

    return a_points;
}

// Define a function for determining if 2 lines are parallel
bool areParallel(const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b, const double tol)
{
    // Make sure each line has 2 points
    assert (a.size() == 2);
    assert (b.size() == 2);

    // Compute the direction vectors
    std::vector<double> a_direction_vec = computeDirectionVector(a[1], a[0]);
    std::vector<double> b_direction_vec = computeDirectionVector(b[1], b[0]);
    
    // Take the cross of them
    std::vector<double> a_cross_b = cross(a_direction_vec, b_direction_vec);
    
    bool sol = (std::abs(a_cross_b[0]) < tol && std::abs(a_cross_b[1]) < tol && std::abs(a_cross_b[2]) < tol);
    
    return sol;
}

// Define a function for determining if 2 surfaces are parallel
bool areParallel(const std::vector<std::vector<std::vector<double>>> a, const std::vector<std::vector<std::vector<double>>> b, const double tol)
{
    // Make sure each surface has at least 3 lines
    assert (a.size() >= 3);
    assert (b.size() >= 3);
    
    // Get 3 points from each surface (to compute the normal vectors)
    std::vector<std::vector<double>> a_points = getThreeSurfacePoints(a);
    std::vector<std::vector<double>> b_points = getThreeSurfacePoints(b);
    
    // Compute the normal vectors
    std::vector<double> normal1 = computeSurfaceNormal_3Point(a_points);
    std::vector<double> normal2 = computeSurfaceNormal_3Point(b_points);
    
    bool sol = (norm(cross(normal1, normal2)) < tol);
    
    return sol;
}

// Define a function for determining if the normal vector of a surface and a vector are parallel
bool areParallel(const std::vector<std::vector<std::vector<double>>> a, const std::vector<double> b, const double tol)
{
    // Make sure each surface has at least 3 lines
    assert (a.size() >= 3);
    assert (b.size() == 3);

    // Get 3 points from the surface (to compute the normal vector)
    std::vector<std::vector<double>> a_points = getThreeSurfacePoints(a);
    
    // Compute the normal vector of the surface
    std::vector<double> normal1 = computeSurfaceNormal_3Point(a_points);
    
    bool sol = (norm(cross(normal1, b)) < tol);
    
    return sol;
}

// Define a function that takes in a line tag, a spatial dimension, and a coordinate, and determines if the line tag coincides with the axis
bool isCoincides(const int ln_tag, const int dim, const double test_val, const double tol)
{
    // Get the coordinates of the end points of the provided line    
    std::vector<std::vector<double>> ln_geo = getLineCoords(ln_tag);

    // Create a dummy line using the test values and the origin to see if the provided line is parallel to it
    std::vector<std::vector<double>> dummy_ln_geo;
    std::vector<double> helper_pt1, helper_pt2;
    helper_pt1 = {0.0, 0.0, 0.0};
    helper_pt1[dim] = test_val;
    helper_pt2 = {1.0, 1.0, 0.0};
    helper_pt2[dim] = test_val;
    dummy_ln_geo.push_back( helper_pt1 );
    dummy_ln_geo.push_back( helper_pt2 );

    // Determine if the given line is parallel to the dummy line
    bool isParallel = areParallel(ln_geo, dummy_ln_geo, tol);

    // Determine if the given line intersects the dummy line
    bool isIntersect = ( compareValues(ln_geo[0][dim], test_val, tol) && compareValues(ln_geo[1][dim], test_val, tol) );
    
    // Return true or false depending on the results
    if (isParallel && isIntersect) { return true; }
    return false;
}

// Define a function that takes in a surface tag, a spatial dimensions (that is parallel with the surface normal), and a coordinate (in that spatial dimension), and determines if the surface tag coincides with the surface perpendicular to the specified dimension at the given coordinate
bool isCoincidesSurface(const int srfc_tag, const int dim, const double test_val, const double tol)
{
    // Get the coordinates of the corner points of the provided surface
    std::vector<std::vector<std::vector<double>>> srfc_geo = getSurfaceCoords(srfc_tag);

    // Create a dummy vector to see if the normal of the provided surface is parallel to it
    std::vector<double> helper_vec;
    helper_vec = {0.0, 0.0, 0.0};
    helper_vec[dim] = 1.0;
    
    // Determine if the normal of the given surface is parallel to the dummy vector
    bool isParallel = areParallel(srfc_geo, helper_vec, tol);
    if (!isParallel) { return false; }
    
    // Determine if the given line intersects the dummy line
    for (size_t i = 0; i < srfc_geo.size(); i++)
    {
        if ( !compareValues(srfc_geo[i][0][dim], test_val, tol) || !compareValues(srfc_geo[i][1][dim], test_val, tol) )
        {
            return false;
        }
    }
    
    return true;
}

// Define a function for comparing the end points of lines
bool compareLineEndPoints(const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b, const double tol)
{
    // Make sure each line has 2 points
    assert (a.size() == 2);
    assert (b.size() == 2);
    bool sol = (compareVectors(a[0], b[0], tol) || compareVectors(a[0], b[1], tol) ||
                compareVectors(a[1], b[0], tol) || compareVectors(a[1], b[1], tol));
    return sol;
}

// Define a function for comparing the corner points of a surface
bool compareSurfaceEndPoints(const std::vector<std::vector<std::vector<double>>> a, const std::vector<std::vector<std::vector<double>>> b, const double tol)
{
    // Make sure each line has 2 points
    for (size_t i = 0; i < a.size(); i++)
    {
        for (size_t j = 0; j < b.size(); j++)
        {
            if (compareLineEndPoints(a[i], b[j])) { return true; }
        }
    }
    return false;
}





// ========================================
// Functions for creating physical groups from input/config file
// ========================================  

// Create the physical groups that combine the left, right, top, bottom boundaries of the domain 
int getPhysicalGroup(int &physGroupTag, const string side, const int dim, std::vector<int> &top_tags, std::vector<int> &right_tags, std::vector<int> &bottom_tags, std::vector<int> &left_tags)
{
    if (side == "none") { physGroupTag = -1; return 0; }
    else if (side == "top") { physGroupTag = gmsh::model::addPhysicalGroup(dim, top_tags); return 1; }
    else if (side == "right") { physGroupTag = gmsh::model::addPhysicalGroup(dim, right_tags); return 2; }
    else if (side == "bottom") { physGroupTag = gmsh::model::addPhysicalGroup(dim, bottom_tags); return 3; }
    else if (side == "left") { physGroupTag = gmsh::model::addPhysicalGroup(dim, left_tags); return 4; }
    else
    {
        std::cout << "CRITICAL ERROR: gmsh_utils.cpp: getPhysicalGroup: Provided side, " << side << ", unidentified." << std::endl;
        return -100;
    }
}




// ========================================
// Functions for simplistic geometry creation
// ========================================  

// Define a function for making a circle from lines
int createCircle(const double x_center, const double y_center, const double z_center, const double r, const double dx_input)
{
    // Restrict the largest dx to be as large as the circle radius (i.e., which draws a hexagon)
    const double dx = (dx_input < r) ? dx_input : r;

    // Define the number of segments on the circle circumferance, and the angle between them
    int N_segments = static_cast<int>(std::ceil(2 * pi * r / dx));
    double arc_ang = 2 * pi / N_segments;

    // Create points and lines on the circumference
    double theta, x, y, z;
    std::vector<int> circ_pts, circ_lns, circ_cl;
    for (int i = 0; i < N_segments; i++)
    {
        theta = i * arc_ang; // In radians
        x = x_center + r * std::cos(theta);
        y = y_center + r * std::sin(theta);
        z = z_center;
        circ_pts.push_back( gmsh::model::occ::addPoint(x, y, z, 1.0) );
        if (i > 0) { circ_lns.push_back( gmsh::model::occ::addLine(circ_pts[i], circ_pts[i - 1]) ); }
    }
    
    // Add the final line
    circ_lns.push_back( gmsh::model::occ::addLine(circ_pts[0], circ_pts.back()) );
    
    // Create the curved loop and surface
    circ_cl.push_back( gmsh::model::occ::addCurveLoop(circ_lns) );
    int circ_srfc = gmsh::model::occ::addPlaneSurface(circ_cl);

    return circ_srfc;
}





// ========================================
// Functions for importing geometries from SVG/text files
// ========================================  

std::vector<int> getCutSurfaceFromFile(const string &geometry_file_path, const double scale)
{
    JSONDict cut_geo_dict, geo_dict;
    cut_geo_dict.loadFromFile(geometry_file_path);
    geo_dict = *cut_geo_dict["geometry"];
    
    int N_surf = (int)(*cut_geo_dict["other"])["N_surfaces"]; // Get the number of surfaces in the geo file

    std::vector<double> x_coords, y_coords;
    std::vector<int> lines, curve_loops, surfs;
    for (int i_surf = 0; i_surf < N_surf; i_surf++)
    {
        // Get coordinate data
        x_coords = (*geo_dict["surface_" + to_string(i_surf)])["x_coords"];
        y_coords = (*geo_dict["surface_" + to_string(i_surf)])["y_coords"];
        assert(x_coords.size() == y_coords.size());

        // Scale coordinate data
        for (int i = 0; i < x_coords.size(); i++) { x_coords[i] *= scale; y_coords[i] *= scale; }
        
        // Make lines, curve loops, and surfaces from coordinate data
        lines = makeLines(x_coords, y_coords);
        curve_loops.push_back( gmsh::model::occ::addCurveLoop(lines) );
        surfs.push_back( gmsh::model::occ::addPlaneSurface(curve_loops) );
    }
    return surfs;
}

std::vector<int> makeLines(const std::vector<double> &x, const std::vector<double> &y)
{
    assert(x.size() == y.size());
    //assert(x.size() == z.size());

    int point1, point2;
    std::vector<int> lines;
    for (int i_pnt = 0; i_pnt < x.size(); i_pnt++)
    {
        point1 = gmsh::model::occ::addPoint(x[i_pnt], y[i_pnt], 0.0, 1.0);
        if (i_pnt < x.size() - 1) { point2 = gmsh::model::occ::addPoint(x[i_pnt + 1], y[i_pnt + 1], 0.0, 1.0); }
        else { point2 = gmsh::model::occ::addPoint(x[0], y[0], 0.0, 1.0); }
        lines.push_back( gmsh::model::occ::addLine(point1, point2) );
    }
    return lines;
}





// ========================================
// Functions for making averaging region arrangments
// ========================================  

void MakeARMesh_Uniform2DRectangular(const int N_AR_x, const double &ell_x, const double &L_x,
    const int N_AR_y, const double &ell_y, const double &L_y, std::vector<std::pair<int, int>> &ARs_uc,
    std::vector<std::vector<std::vector<std::vector<double>>>> &ARs_uc_geo)
{
    for (int j = 0; j < N_AR_y; ++j)
    {
        for (int i = 0; i < N_AR_x; ++i)
        {
            // Define the bottom left corner of the AR
            double x = i * ell_x - L_x/2;
            double y = j * ell_y - L_y/2;

            // Create the AR
            int AR = gmsh::model::occ::addRectangle(x, y, 0.0, ell_x, ell_y);
            
            // Synchronize the model
            gmsh::model::occ::synchronize();
            
            // Store the AR so that the geometry may be cut into it later
            ARs_uc.push_back({2, AR}); // (dim = 2 for surface)

            // Store the coordinates of the AR's points so that we can distinguish the AR interfaces from 
            // geometry surfaces after cutting the AR. ARs_uc_geo is:
            //      2D: ARs_uc_geo[AR, lines, points, coords]
            ARs_uc_geo.push_back(getSurfaceCoords(AR));
        }
    }
}

void MakeARMesh_Uniform3DRectangular(const int N_AR_x, const double &ell_x, const double &L_x,
    const int N_AR_y, const double &ell_y, const double &L_y,
    const int N_AR_z, const double &ell_z, const double &L_z,
    std::vector<std::pair<int, int>> &ARs_uc,
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &ARs_uc_geo)
{
    for (int k = 0; k < N_AR_z; ++k)
    {
        for (int j = 0; j < N_AR_y; ++j)
        {
            for (int i = 0; i < N_AR_x; ++i)
            {
                // Define the bottom left corner of the AR
                double x = i * ell_x - L_x/2;
                double y = j * ell_y - L_y/2;
                double z = k * ell_z - L_z/2;

                // Create the AR
                int AR = gmsh::model::occ::addBox(x, y, z, ell_x, ell_y, ell_z);
                
                // Synchronize the model
                gmsh::model::occ::synchronize();
                
                // Store the AR so that the geometry may be cut into it later
                ARs_uc.push_back({3, AR}); // (dim = 2 for surface)

                // Store the coordinates of the AR's points so that we can distinguish the AR interfaces from 
                // geometry surfaces after cutting the AR. ARs_uc_geo is:
                //      3D: ARs_uc_geo[AR, surface, lines, points, coords]
                ARs_uc_geo.push_back(getVolumeCoords(AR));
            }
        }
    }
}
    
