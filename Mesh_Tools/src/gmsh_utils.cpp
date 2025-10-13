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

// Define a function that takes in a line tag and returns the minimum/maximum x/y coordinates of the line's boundary points.
std::vector<double> getMinMaxLineCoords(int tag)
{
    std::vector<double> min_max_x_y = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // min x, max x, min y, max y, min z, max z ---> min_max_x_y[0], min_max_x_y[1], min_max_x_y[2], min_max_x_y[3], min_max_x_y[4], min_max_x_y[5] 
    std::vector<std::vector<double>> coords = getLineCoords(tag);
    for (int i_pt = 0; i_pt < coords.size(); i_pt++)
    {
        if (i_pt == 0)
        {
            min_max_x_y[0] = coords[i_pt][0]; min_max_x_y[1] = coords[i_pt][0];
            min_max_x_y[2] = coords[i_pt][1]; min_max_x_y[3] = coords[i_pt][1];
            min_max_x_y[4] = coords[i_pt][2]; min_max_x_y[5] = coords[i_pt][2];
        }
        if (coords[i_pt][0] < min_max_x_y[0]) { min_max_x_y[0] = coords[i_pt][0]; }
        else if (coords[i_pt][0] > min_max_x_y[1]) { min_max_x_y[1] = coords[i_pt][0]; }
        if (coords[i_pt][1] < min_max_x_y[2]) { min_max_x_y[2] = coords[i_pt][1]; }
        else if (coords[i_pt][1] > min_max_x_y[3]) { min_max_x_y[3] = coords[i_pt][1]; }
        if (coords[i_pt][2] < min_max_x_y[4]) { min_max_x_y[4] = coords[i_pt][2]; }
        else if (coords[i_pt][2] > min_max_x_y[5]) { min_max_x_y[5] = coords[i_pt][2]; }
    }
    return min_max_x_y;
}

// Define a function that takes in a surface tag and returns the minimum/maximum x/y coordinates of the surface's boundary points.
std::vector<double> getMinMaxSurfaceCoords(int tag)
{
    std::vector<double> min_max_x_y = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // min x, max x, min y, max y, min z, max z ---> min_max_x_y[0], min_max_x_y[1], min_max_x_y[2], min_max_x_y[3], min_max_x_y[4], min_max_x_y[5] 
    std::vector<std::vector<std::vector<double>>> coords = getSurfaceCoords(tag);
    for (int i_ln = 0; i_ln < coords.size(); i_ln++)
    {
        for (int i_pt = 0; i_pt < coords[i_ln].size(); i_pt++)
        {
            if (i_ln == 0 && i_pt == 0)
            {
                min_max_x_y[0] = coords[i_ln][i_pt][0]; min_max_x_y[1] = coords[i_ln][i_pt][0];
                min_max_x_y[2] = coords[i_ln][i_pt][1]; min_max_x_y[3] = coords[i_ln][i_pt][1];
                min_max_x_y[4] = coords[i_ln][i_pt][2]; min_max_x_y[5] = coords[i_ln][i_pt][2];
            }
            if (coords[i_ln][i_pt][0] < min_max_x_y[0]) { min_max_x_y[0] = coords[i_ln][i_pt][0]; }
            else if (coords[i_ln][i_pt][0] > min_max_x_y[1]) { min_max_x_y[1] = coords[i_ln][i_pt][0]; }
            if (coords[i_ln][i_pt][1] < min_max_x_y[2]) { min_max_x_y[2] = coords[i_ln][i_pt][1]; }
            else if (coords[i_ln][i_pt][1] > min_max_x_y[3]) { min_max_x_y[3] = coords[i_ln][i_pt][1]; }
            if (coords[i_ln][i_pt][2] < min_max_x_y[4]) { min_max_x_y[4] = coords[i_ln][i_pt][2]; }
            else if (coords[i_ln][i_pt][2] > min_max_x_y[5]) { min_max_x_y[5] = coords[i_ln][i_pt][2]; }
        }
    }
    return min_max_x_y;
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

bool areAdjacentEntities(const int entity_dim, const int A_tag, const int B_tag)
{
    // Get the boundary entities for each provided entity
    std::vector<std::pair<int, int>> A_boundaries, B_boundaries;
    gmsh::model::getBoundary({{entity_dim, A_tag}}, A_boundaries, false, false, false);
    gmsh::model::getBoundary({{entity_dim, B_tag}}, B_boundaries, false, false, false);
    
    // Get the boundary tags
    std::vector<int> A_boundary_tags; for (int i = 0; i < A_boundaries.size(); i++) { A_boundary_tags.push_back( A_boundaries[i].second ); }
    std::vector<int> B_boundary_tags; for (int i = 0; i < B_boundaries.size(); i++) { B_boundary_tags.push_back( B_boundaries[i].second ); }

    // Compare the boundary tags
    for (int i = 0; i < A_boundary_tags.size(); i++)
    {
        for (int j = 0; j < B_boundary_tags.size(); j++)
        {
            if (A_boundary_tags[i] == B_boundary_tags[j]) { return true; }
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
// Functions for creating simple geometries
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

void load2DToolGeo(std::vector<std::pair<int, int>> &tool_sfs, std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo, const double geo_scale, const string geometry_file_path)
{
    std::vector<int> cut_geo = getCutSurfaceFromFile(geometry_file_path, geo_scale);

    gmsh::model::occ::synchronize();
    
    for (int i_surf = 0; i_surf < cut_geo.size(); i_surf++)
    {
        tool_sfs.push_back({2, cut_geo[i_surf]});
        tool_geo.push_back( getSurfaceCoords(cut_geo[i_surf]) );
    }
}
void load2DToolGeo(std::vector<std::pair<int, int>> &tool_sfs, std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo, const double geo_scale, const double x, const double y, const double z)
{
    std::vector<int> cut_geo;
    cut_geo.push_back( createCircle(x, y, z, 0.3 * geo_scale, 0.05 * geo_scale) );
    
    gmsh::model::occ::synchronize();
    
    for (int i_surf = 0; i_surf < cut_geo.size(); i_surf++)
    {
        tool_sfs.push_back({2, cut_geo[i_surf]});
        tool_geo.push_back( getSurfaceCoords(cut_geo[i_surf]) );
    }
}





// ========================================
// Functions for making/editting averaging region arrangments
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

bool createMergedARGroups(const std::vector<std::pair<int, int>> &ARs,
    double min_area_threshold, double max_AR_length, double max_AR_length_ratio,
    std::vector<std::vector<int>> &AR_tags_toBeMerged, std::vector<double> &AR_pore_space_areas)
{
    // In AR_tags_toBeMerged, the inner list is {merge host tag (i.e., the big AR), merging tag (i.e., the smoll AR)}
    bool anything_Merged = false;

    for (int i_ar = 0; i_ar < ARs.size(); i_ar++)
    {
        // Get AR information
        int dim = ARs[i_ar].first;
        int tag = ARs[i_ar].second;

        // Compute AR area
        double AR_pore_space_area;
        gmsh::model::occ::getMass(dim, tag, AR_pore_space_area);

        // Decide if the AR must be merged (i.e., it is smaller than the threshold)
        if (AR_pore_space_area < min_area_threshold)
        {
            // Get the max/min x/y coordinates of the AR
            std::vector<double> min_max_x_y = getMinMaxSurfaceCoords(tag); // min x, max x, min y, max y, min z, max z ---> min_max_x_y[0], min_max_x_y[1], min_max_x_y[2], min_max_x_y[3], min_max_x_y[4], min_max_x_y[5] 
            
            // See if a "host" AR can be found to merge the small AR to
            bool found_host = false;
            for (int j_ar = 0; j_ar < ARs.size(); j_ar++)
            {
                if (i_ar == j_ar) { continue; } // AR cannot host itself
                
                // Get the tag of the potential host AR
                int tag2 = ARs[j_ar].second;

                // If host AR, AR[j_ar], is not adjacent to small AR, AR[i_ar], try the next AR
                if (!areAdjacentEntities(dim, tag, tag2)) { continue; }
                
                // For now, we require the host AR to satisfy the min_area_threshold itself; otherwise, try the next AR
                double ARj_area;
                gmsh::model::occ::getMass(dim, tag2, ARj_area);
                if (ARj_area < min_area_threshold) { continue; }


                // TODO: check AR_tags_toBeMerged to see if host AR has already been merged. Then, consider merged host AR
                //       as a potential host, not just AR[j_ar]


                // If potential host AR is adjacent and has large enough area, get the max/min x/y coordinates of the host AR
                std::vector<double> new_min_max_x_y = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                std::vector<double> min_max_x_y2 = getMinMaxSurfaceCoords(tag2);
                new_min_max_x_y[0] = min_max_x_y[0] < min_max_x_y2[0] ? min_max_x_y[0] : min_max_x_y2[0];
                new_min_max_x_y[1] = min_max_x_y[1] > min_max_x_y2[1] ? min_max_x_y[1] : min_max_x_y2[1];
                new_min_max_x_y[2] = min_max_x_y[2] < min_max_x_y2[2] ? min_max_x_y[2] : min_max_x_y2[2];
                new_min_max_x_y[3] = min_max_x_y[3] > min_max_x_y2[3] ? min_max_x_y[3] : min_max_x_y2[3];
                new_min_max_x_y[4] = min_max_x_y[4] < min_max_x_y2[4] ? min_max_x_y[4] : min_max_x_y2[4];
                new_min_max_x_y[5] = min_max_x_y[5] > min_max_x_y2[5] ? min_max_x_y[5] : min_max_x_y2[5];
                

                // TODO: again, host AR might be a merged AR. Get the min max x y for the merged AR (i.e., AR_tags_toBeMerged)


                // If merging the host and small ARs would make the new AR too big, try the next AR
                if ( new_min_max_x_y[1] - new_min_max_x_y[0] > max_AR_length ||
                     new_min_max_x_y[3] - new_min_max_x_y[2] > max_AR_length ||
                     new_min_max_x_y[5] - new_min_max_x_y[4] > max_AR_length ) { continue; }
                
                // If merging the host and small ARs would create an AR with too large a length scale ratio, try the next AR
                if (new_min_max_x_y[4] == 0.0 && new_min_max_x_y[5] == 0.0)
                {
                    if ( (new_min_max_x_y[1] - new_min_max_x_y[0])/(new_min_max_x_y[3] - new_min_max_x_y[2]) >= max_AR_length_ratio ||
                         (new_min_max_x_y[3] - new_min_max_x_y[2])/(new_min_max_x_y[1] - new_min_max_x_y[0]) >= max_AR_length_ratio ) { continue; }
                }
                else
                {
                    if ( (new_min_max_x_y[1] - new_min_max_x_y[0])/(new_min_max_x_y[3] - new_min_max_x_y[2]) >= max_AR_length_ratio ||
                         (new_min_max_x_y[3] - new_min_max_x_y[2])/(new_min_max_x_y[1] - new_min_max_x_y[0]) >= max_AR_length_ratio ||
                         (new_min_max_x_y[1] - new_min_max_x_y[0])/(new_min_max_x_y[5] - new_min_max_x_y[4]) >= max_AR_length_ratio ||
                         (new_min_max_x_y[5] - new_min_max_x_y[4])/(new_min_max_x_y[1] - new_min_max_x_y[0]) >= max_AR_length_ratio ||
                         (new_min_max_x_y[3] - new_min_max_x_y[2])/(new_min_max_x_y[5] - new_min_max_x_y[4]) >= max_AR_length_ratio ||
                         (new_min_max_x_y[5] - new_min_max_x_y[4])/(new_min_max_x_y[3] - new_min_max_x_y[2]) >= max_AR_length_ratio ) { continue; }
                }
                
                
                // If everything is satisfied, merge the host and small ARs
                bool new_entry = true;
                for (int i = 0; i < AR_tags_toBeMerged.size(); i++)
                {
                    if (AR_tags_toBeMerged[i][0] == tag2)
                    {
                        AR_tags_toBeMerged[i].push_back( tag );
                        AR_pore_space_areas[i] += AR_pore_space_area;
                        new_entry = false;
                        break;
                    }
                }
                if (new_entry) { AR_tags_toBeMerged.push_back( {tag2, tag} ); AR_pore_space_areas.push_back( ARj_area + AR_pore_space_area ); }
                
                found_host = true;
                anything_Merged = true;
                break;
            }
            
            if (!found_host) { AR_tags_toBeMerged.push_back( {tag} ); AR_pore_space_areas.push_back( AR_pore_space_area ); }

        }
        else
        {
            bool new_entry = true;
            for (int i = 0; i < AR_tags_toBeMerged.size(); i++)
            {
                if (AR_tags_toBeMerged[i][0] == tag) { new_entry = false; break; }
            }
            if (new_entry) { AR_tags_toBeMerged.push_back( {tag} ); AR_pore_space_areas.push_back( AR_pore_space_area ); }
        }
    }

    return anything_Merged;
}







// ========================================
// Functions for separating averaging region lines from cut geometry lines 
// ========================================  

void separateARandGeometrySurfaceLines(const std::vector<std::pair<int, int>> &ARs,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo,
    std::vector<std::vector<std::pair<int, int>>> &ARs_lns,
    std::vector<std::vector<std::pair<int, int>>> &tool_sfs_lns)
{
    // Define necessary variables
    std::vector<std::pair<int, int>> AR_all_lns;
    std::vector<std::vector<std::vector<double>>> AR_c_geo;
    int sf_ind;

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
            bool isGeoLine = false;
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
            if (isGeoLine) { tool_sfs_lns[sf_ind].push_back(AR_all_lns[i_ln_c]); } // NOTE: tool_sfs_lns retains the same surfaces as tool_geo. There may be more lines in each surface though
            else { ARs_lns[i_ar].push_back(AR_all_lns[i_ln_c]); } // Assuming that the remainder lines are the AR boundaries
        }
    }
}

void separateARandGeometryVolumeSurfaces(const std::vector<std::pair<int, int>> &ARs,
    const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &tool_geo_3D,
    std::vector<std::vector<std::pair<int, int>>> &ARs_srfcs,
    std::vector<std::vector<std::pair<int, int>>> &tool_vol_srfcs)
{
    // Define necessary variables
    std::vector<std::pair<int, int>> AR_all_srfcs;
    std::vector<std::vector<std::vector<std::vector<double>>>> AR_c_geo_3D;
    int vol_ind;

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
            bool isGeoSurface = false;
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
            if (isGeoSurface) { tool_vol_srfcs[vol_ind].push_back(AR_all_srfcs[i_srfc_c]); } // NOTE: tool_vol_srfcs retains the same volumes as tool_geo. There may be more surfaces in each volume though
            else { ARs_srfcs[i_ar].push_back(AR_all_srfcs[i_srfc_c]); } // Assuming that the remainder surfaces are the AR boundaries
        }
    }
}


