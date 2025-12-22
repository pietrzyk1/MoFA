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
        ln_pts_coords.push_back( getPointCoords(pts_dimtag[i_pts].second) );
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
        srfc_ln_pts_coords.push_back( getLineCoords(lns_dimtag[i_ln].second) );
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
        vol_srfc_ln_pts_coords.push_back( getSurfaceCoords(srfc_dimtag[i_srfc].second) );
    }

    return vol_srfc_ln_pts_coords;
}

/*
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
*/

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

// Define a function that takes in an entity tag and returns the minimum/maximum x/y/z coordinates of the entity's boundary points.
std::vector<double> getEntityMinMaxCoords(const std::pair<int, int> &dimtag)
{
    // Initialize variables
    std::vector<double> min_max_x_y_z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // min x, max x, min y, max y, min z, max z ---> min_max_x_y[0], min_max_x_y[1], min_max_x_y[2], min_max_x_y[3], min_max_x_y[4], min_max_x_y[5] 
    std::vector<std::vector<double>> coords;
    
    // Get the point coordinates from the provided entity
    if (dimtag.first == 0 || dimtag.first == 1) {
        cerr << "gmsh_utils.cpp: getEntityMinMaxCoords(): CRITICAL ERROR: Entity dimension is 0 or 1, but the function is currently written for dimensions of 2 or 3. More code is needed for these cases." << endl; exit(1); }
    else if (dimtag.first == 2) {
        std::vector<std::vector<std::vector<double>>> entity_coords = getSurfaceCoords(dimtag.second);
        for (int i_ln = 0; i_ln < entity_coords.size(); i_ln++) {
            for (int i_pt = 0; i_pt < entity_coords[i_ln].size(); i_pt++) {
                coords.push_back( entity_coords[i_ln][i_pt] ); } } }
    else if (dimtag.first == 3) {
        std::vector<std::vector<std::vector<std::vector<double>>>> entity_coords = getVolumeCoords(dimtag.second);
        for (int i_surf = 0; i_surf < entity_coords.size(); i_surf++) {
            for (int i_ln = 0; i_ln < entity_coords[i_surf].size(); i_ln++) {
                for (int i_pt = 0; i_pt < entity_coords[i_surf][i_ln].size(); i_pt++) {
                    coords.push_back( entity_coords[i_surf][i_ln][i_pt] ); } } } }
    else {
        cerr << "gmsh_utils.cpp: getEntityMinMaxCoords(): CRITICAL ERROR: Entity dimension unrecognized." << endl; exit(1); }

    // Determine the minimum and maximum coordinate values of the entity points
    for (int i_pt = 0; i_pt < coords.size(); i_pt++)
    {
        if (i_pt == 0) {
            min_max_x_y_z[0] = coords[i_pt][0]; min_max_x_y_z[1] = coords[i_pt][0];
            min_max_x_y_z[2] = coords[i_pt][1]; min_max_x_y_z[3] = coords[i_pt][1];
            min_max_x_y_z[4] = coords[i_pt][2]; min_max_x_y_z[5] = coords[i_pt][2];
        }
        if (coords[i_pt][0] < min_max_x_y_z[0]) { min_max_x_y_z[0] = coords[i_pt][0]; }
        else if (coords[i_pt][0] > min_max_x_y_z[1]) { min_max_x_y_z[1] = coords[i_pt][0]; }
        if (coords[i_pt][1] < min_max_x_y_z[2]) { min_max_x_y_z[2] = coords[i_pt][1]; }
        else if (coords[i_pt][1] > min_max_x_y_z[3]) { min_max_x_y_z[3] = coords[i_pt][1]; }
        if (coords[i_pt][2] < min_max_x_y_z[4]) { min_max_x_y_z[4] = coords[i_pt][2]; }
        else if (coords[i_pt][2] > min_max_x_y_z[5]) { min_max_x_y_z[5] = coords[i_pt][2]; }
    }

    return min_max_x_y_z;
}





// ========================================
// Functions for mathematical calculations
// ========================================  

// Compute the cross product of two 3D vectors
std::vector<double> cross(const std::vector<double> &a, const std::vector<double> &b)
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
bool compareVectors(const std::vector<double> &a, const std::vector<double> &b, const double tol)
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
    //for (int i = 0; i < A_boundary_tags.size(); i++) {
    //    for (int j = 0; j < B_boundary_tags.size(); j++) {
    //        if (A_boundary_tags[i] == B_boundary_tags[j]) { return true; }
    //    }
    //}
    //return false;
    
    // Another, potentially faster way:
    A_boundary_tags.insert(A_boundary_tags.end(), B_boundary_tags.begin(), B_boundary_tags.end());
    std::sort(A_boundary_tags.begin(), A_boundary_tags.end());
    for (int i = 0; i < A_boundary_tags.size() - 1; i++) { if (A_boundary_tags[i] == A_boundary_tags[i + 1]) { return true; } }
    return false;
}
bool areAdjacentEntities(const std::vector<std::pair<int, int>> &A_in_dimtag, const std::vector<std::pair<int, int>> &B_in_dimtag,
    std::vector<std::pair<int, int>> &A_out_dimtag, std::vector<std::pair<int, int>> &B_out_dimtag)
{
    // Get the boundary entities for each provided dim-tag pair
    gmsh::model::getBoundary(A_in_dimtag, A_out_dimtag, false, false, false);
    gmsh::model::getBoundary(B_in_dimtag, B_out_dimtag, false, false, false);
    
    // Remove duplicates from A_out_dimtag and B_out_dimtag before continuing
    std::sort(A_out_dimtag.begin(), A_out_dimtag.end(),
        [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.second < b.second; });
    auto A_last_unique = std::unique(A_out_dimtag.begin(), A_out_dimtag.end());
    A_out_dimtag.erase(A_last_unique, A_out_dimtag.end());
    std::sort(B_out_dimtag.begin(), B_out_dimtag.end(),
        [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.second < b.second; });
    auto B_last_unique = std::unique(B_out_dimtag.begin(), B_out_dimtag.end());
    B_out_dimtag.erase(B_last_unique, B_out_dimtag.end());

    // Get the boundary tags
    std::vector<int> A_boundary_tags; for (int i = 0; i < A_out_dimtag.size(); i++) { A_boundary_tags.push_back( A_out_dimtag[i].second ); }
    std::vector<int> B_boundary_tags; for (int i = 0; i < B_out_dimtag.size(); i++) { B_boundary_tags.push_back( B_out_dimtag[i].second ); }
    
    // Append boundary tags, sort them, and see if there are any duplicates. If so, then entities are shared
    A_boundary_tags.insert(A_boundary_tags.end(), B_boundary_tags.begin(), B_boundary_tags.end());
    std::sort(A_boundary_tags.begin(), A_boundary_tags.end());
    auto last_unique = std::unique(A_boundary_tags.begin(), A_boundary_tags.end());
    if (last_unique == A_boundary_tags.end()) { return false; } else { return true; }
    
    //for (int i = 0; i < A_boundary_tags.size() - 1; i++) { if (A_boundary_tags[i] == A_boundary_tags[i + 1]) { return true; } }
    //return false;
}

bool shareAnyBoundaryEntities(const int entity_dim, const int A_tag, const int B_tag)
{
    // Make sure the tags provided are lines, surfaces, or volumes
    assert (entity_dim > 0);

    // Initiate vectors for the input and output dim-tag pairs
    std::vector<std::pair<int, int>> A_input_entities = {{entity_dim, A_tag}}, B_input_entities = {{entity_dim, B_tag}};
    std::vector<std::pair<int, int>> A_output_entities, B_output_entities;

    // See if any of the boundary entities of input entity(s) A are shared with input entity(s) B. If so, return true
    if (areAdjacentEntities(A_input_entities, B_input_entities, A_output_entities, B_output_entities)) { return true; }
    
    // If there are no shared entities, ensure that a lower dimension of boundary entities can be compared and continue. If there isn't a lower dimension of boundary entities, return false 
    if (A_output_entities[0].first == 0) { return false; }
    A_input_entities.clear(); B_input_entities.clear();

    // See if any of the boundary entities of input entity(s) A are shared with input entity(s) B. If so, return true
    if (areAdjacentEntities(A_output_entities, B_output_entities, A_input_entities, B_input_entities)) { return true; }
    
    // If there are no shared entities, ensure that a lower dimension of boundary entities can be compared and continue. If there isn't a lower dimension of boundary entities, return false 
    if (A_input_entities[0].first == 0) { return false; }
    A_output_entities.clear(); B_output_entities.clear();
    
    // See if any of the boundary entities of input entity(s) A are shared with input entity(s) B. If so, return true
    if (areAdjacentEntities(A_input_entities, B_input_entities, A_output_entities, B_output_entities)) { return true; }

    // If there are no shared entities, ensure that a lower dimension of boundary entities can be compared and continue. If there isn't a lower dimension of boundary entities, return false 
    if (A_output_entities[0].first == 0) { return false; }
    
    cerr << "gmsh_utils.cpp: shareAnyBoundaryEntities(): CRITICAL ERROR: The code should not get to this. If it does, it means an entity dimension greater than 3 was initially provided." << endl;
    exit(1);
}

bool linesIntersect(const std::vector<std::vector<double>>& ln1, const std::vector<std::vector<double>>& ln2, double tol)
{

    double c1 = cross({ln1[1][0] - ln1[0][0], ln1[1][1] - ln1[0][1], 0.0}, {ln2[0][0] - ln1[0][0], ln2[0][1] - ln1[0][1], 0.0})[2];
    double c2 = cross({ln1[1][0] - ln1[0][0], ln1[1][1] - ln1[0][1], 0.0}, {ln2[1][0] - ln1[0][0], ln2[1][1] - ln1[0][1], 0.0})[2];
    double c3 = cross({ln2[1][0] - ln2[0][0], ln2[1][1] - ln2[0][1], 0.0}, {ln1[0][0] - ln2[0][0], ln1[0][1] - ln2[0][1], 0.0})[2];
    double c4 = cross({ln2[1][0] - ln2[0][0], ln2[1][1] - ln2[0][1], 0.0}, {ln1[1][0] - ln2[0][0], ln1[1][1] - ln2[0][1], 0.0})[2];
    
    // Only looks at full intersections (i.e., if lines "straddle"; if only the end points of line A touch line B, or if lines partially overlap, it does not count)
    if ( (c1 > tol && c2 < -tol) || (c1 < -tol && c2 > tol) ) {
        if ( (c3 > tol && c4 < -tol) || (c3 < -tol && c4 > tol) ) {
            return true; } }
    return false;
}

bool isCollinear(const std::vector<std::vector<double>> &ln_pts, const std::vector<double> &pt, double tol)
{
    assert (ln_pts.size() == 2);
    std::vector<std::vector<double>> dist;
    for (int i = 0; i < ln_pts.size(); i++) {
        // Compute the displacement vector between the line point and the point being checked
        dist.push_back( {pt[0] - ln_pts[i][0], pt[1] - ln_pts[i][1], pt[2] - ln_pts[i][2]} );
        
        // Get the magnitude of the displacement vector
        double dist_mag = norm(dist[i]);
        // If pt is one of the points in ln_pts, it is considered collinear
        if (compareValues(0.0, dist_mag, tol)) { return true; }
        
        // Normalize the displacement vector
        for (int j = 0; j < dist[i].size(); j++) { dist[i][j] /= dist_mag; }
    }
    
    // See if the normalized displacement vectors are the same
    if (compareVectors(dist[0], dist[1], tol)) { return true; }
    
    // If not, try multiplying one of the displacement vectors by -1
    for (int j = 0; j < dist[0].size(); j++) { dist[0][j] *= -1.0; }
    if (compareVectors(dist[0], dist[1], tol)) { return true; }
    return false;
    
    // Old stuff
    //if (!compareValues(pt[1], (ln_pts[1][1] - ln_pts[0][1]) / (ln_pts[1][0] - ln_pts[0][0]) * (pt[0] - ln_pts[0][0]) + ln_pts[0][1], tol)) { return false; }
    //if (!compareValues(pt[2], (ln_pts[1][2] - ln_pts[0][2]) / (ln_pts[1][0] - ln_pts[0][0]) * (pt[0] - ln_pts[0][0]) + ln_pts[0][2], tol)) { return false; }
    //return true;
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
        std::cout << "CRITICAL ERROR: gmsh_utils.cpp: getPhysicalGroup(): Provided side, " << side << ", unidentified." << std::endl;
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
        tool_sfs.push_back( {2, cut_geo[i_surf]} );
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














/*
std::vector<int> assignARIndex(const std::vector<double> &coord, const std::vector<std::vector<double>> &AR_planes, const double tol)
{
    // Initiate variables
    std::vector<int> AR_inds = {0};
    std::vector<int> N_boxes; for (int i = 0; i < AR_planes.size(); i++) { N_boxes.push_back( (int)AR_planes[i].size() - 1 ); }
    int box_multiplier = 1;

    for (int i_d = 0; i_d < AR_planes.size(); i_d++)
    {
        // Update the box multiplier (i.e., multiply by N_boxes[0] if considering the i_d = 1, N_boxes[0] * N_boxes[1] if considering the i_d = 2)
        if (i_d > 0) { box_multiplier *= N_boxes[i_d - 1]; }

        bool foundAR = false;
        for (int i_ARp = 0; i_ARp < AR_planes[i_d].size(); i_ARp++) {
            // Consider the case of the point lying on the AR cut planes independently
            if (compareValues(coord[i_d], AR_planes[i_d][i_ARp], tol)) {
                if (i_ARp == 0) { foundAR = true; }
                else if (i_ARp == AR_planes[i_d].size() - 1) { for (int i = 0; i < AR_inds.size(); i++) { AR_inds[i] += (i_ARp - 1) * box_multiplier; } foundAR = true; }
                else {
                    int N_AR_old = AR_inds.size();
                    for (int i = 0; i < N_AR_old; i++) {
                        AR_inds.push_back(AR_inds[i]);
                        AR_inds[i] += (i_ARp - 1) * box_multiplier;
                        AR_inds[AR_inds.size() - 1] += (i_ARp) * box_multiplier; }
                    foundAR = true; }
                break; }
            // Consider the case of the point lying inside an AR
            else if (coord[i_d] < AR_planes[i_d][i_ARp]) {
                assert (i_ARp > 0);
                for (int i = 0; i < AR_inds.size(); i++) { AR_inds[i] += (i_ARp - 1) * box_multiplier; }
                foundAR = true;
                break; }
        }
        assert (foundAR);
    }

    return AR_inds;
}


template <typename T>
class DisplacedVector
{
private:
    int min_ind = 0;
    std::vector<T> data;
public:
    // Class constructors
    DisplacedVector() {}
    DisplacedVector(int min_ind_) : min_ind(min_ind_) {}

    T &operator[](const int i) {
        if (i < min_ind) { cerr << "DisplacedVector: provided index " << i << " is too small for DisplacedVector with minimum index " << min_ind << "." << endl; exit(1); }
        if (i > min_ind + data.size() - 1) { cerr << "DisplacedVector: provided index " << i << " is too large for DisplacedVector with maximum index " << min_ind + data.size() - 1 << "." << endl; exit(1); }
        return data[i - min_ind]; }
    void SetMinInd(const int i) { min_ind = i; }
    void push_back(const T info) { data.push_back(info); }
    int size() { return data.size(); }
    int get_min_ind() { return min_ind; }
    int get_max_ind() { return min_ind + data.size() - 1; }
};

template <typename T, typename lamfunc>
void removeDuplicates(std::vector<T> &vec, lamfunc func)
{
    std::sort(vec.begin(), vec.end(), func);
    auto last = std::unique(vec.begin(), vec.end());
    vec.erase(last, vec.end());
}

template <typename T, typename lamfunc>
void getUniqueEntries(std::vector<T> &vec, lamfunc func)
{
    std::vector<T> new_vec;
    std::sort(vec.begin(), vec.end(), func);
    bool skip = false;
    for (int i = 0; i < vec.size(); i++) {
        if (skip) { skip = false; continue; }
        if (i == vec.size() - 1) { new_vec.push_back( vec[i] ); continue; }
        if (vec[i] == vec[i + 1]) { skip = true; continue; }
        new_vec.push_back( vec[i] );
    }
    vec = new_vec;
}

class BoundingSurface
{
public:
    int dim; // The dimension parallel to the surface's normal vector
    double coord; // The coordinate value (on direction dim) which defines the surface

    std::vector<int> alt_dim;
    std::vector<vector<double>> alt_dim_min_max;


    // Class constructors
    BoundingSurface() {}
    BoundingSurface(const int &dim_, const double &coord_, const std::vector<int> &alt_dim_, const std::vector<vector<double>> &alt_dim_min_max_) :
    dim(dim_), coord(coord_), alt_dim(alt_dim_), alt_dim_min_max(alt_dim_min_max_) {}

    void Set(const int &dim_, const double &coord_, const std::vector<int> &alt_dim_, const std::vector<vector<double>> &alt_dim_min_max_)
    { dim = dim_; coord = coord_; alt_dim = alt_dim_; alt_dim_min_max = alt_dim_min_max_; }
};


class ARGridManager
{
public:
    double tol;
    std::vector<std::vector<double>> AR_planes;
    std::vector<std::tuple<std::vector<double>, bool, int>> AR_corner_pts; // in the tuple: get<0> = coords, get<1> = is already created, get<0> = gmsh tag (if it is already created)
    std::vector<int> isInside; // 0 = outside, 1 = inside, 2 = on the border (2 means the point has also already been created)

    std::vector<std::vector<std::vector<std::vector<double>>>> pts_on_xy_AR_planes, pts_on_xz_AR_planes, pts_on_yz_AR_planes;
    

    // Class constructors
    ARGridManager(const std::vector<std::vector<double>> &AR_planes_, const double tol_) : AR_planes(AR_planes_), tol(tol_) {
        // Initialize AR_corner_pts
        for (int i_z = 0; i_z < AR_planes[2].size(); i_z++) {
            for (int i_y = 0; i_y < AR_planes[1].size(); i_y++) {
                for (int i_x = 0; i_x < AR_planes[0].size(); i_x++) {
                    AR_corner_pts.push_back( std::tuple<std::vector<double>, bool, int>() );
                    std::get<0>(AR_corner_pts[AR_corner_pts.size() - 1]) = {AR_planes[0][i_x], AR_planes[1][i_y], AR_planes[2][i_z]};
                    std::get<1>(AR_corner_pts[AR_corner_pts.size() - 1]) = false;
                    std::get<2>(AR_corner_pts[AR_corner_pts.size() - 1]) = -1;
                }
            }
        }

        // Initialize pts_on_xy_AR_planes
        for (int i_x = 0; i_x < AR_planes[0].size(); i_x++) {
            pts_on_xy_AR_planes.push_back( std::vector<std::vector<std::vector<double>>>() );
            for (int i_y = 0; i_y < AR_planes[1].size(); i_y++) {
                pts_on_xy_AR_planes[i_x].push_back( std::vector<std::vector<double>>() );
            }
        }
        // Initialize pts_on_xz_AR_planes
        for (int i_x = 0; i_x < AR_planes[0].size(); i_x++) {
            pts_on_xz_AR_planes.push_back( std::vector<std::vector<std::vector<double>>>() );
            for (int i_z = 0; i_z < AR_planes[2].size(); i_z++) {
                pts_on_xz_AR_planes[i_x].push_back( std::vector<std::vector<double>>() );
            }
        }
        // Initialize pts_on_yz_AR_planes
        for (int i_y = 0; i_y < AR_planes[1].size(); i_y++) {
            pts_on_yz_AR_planes.push_back( std::vector<std::vector<std::vector<double>>>() );
            for (int i_z = 0; i_z < AR_planes[2].size(); i_z++) {
                pts_on_yz_AR_planes[i_y].push_back( std::vector<std::vector<double>>() );
            }
        }
    }
    

    void checkAndAddPoint(const std::pair<std::vector<double>, int> &coords_gpt) {
        checkAndAddForARCornerPoint(coords_gpt);
        checkAndAddForOnRayCast(coords_gpt);
        //if (compareVectors(coords_gpt.first, {-5, -5, -3})) {
        //    cout << "asdfasdfasdf" << endl;
        //    exit(1);
        //}
    }

    void checkAndAddForARCornerPoint(const std::pair<std::vector<double>, int> &coords_gpt) {
        for (int i_pt = 0; i_pt < AR_corner_pts.size(); i_pt++) {
            if (compareVectors(coords_gpt.first, std::get<0>(AR_corner_pts[i_pt]), tol)) {
                std::get<1>(AR_corner_pts[i_pt]) = true;
                std::get<2>(AR_corner_pts[i_pt]) = coords_gpt.second;
                //cout << coords_gpt.first[0] << " " << coords_gpt.first[1] << " " << coords_gpt.first[2] << endl;
                break;
            }
        }
    }

    void checkAndAddForOnRayCast(const std::pair<std::vector<double>, int> &coords_gpt) {
        for (int i_x = 0; i_x < AR_planes[0].size(); i_x++) {
            if (compareValues(coords_gpt.first[0], AR_planes[0][i_x], tol)) {
                for (int i_y = 0; i_y < AR_planes[1].size(); i_y++) {
                    if (compareValues(coords_gpt.first[1], AR_planes[1][i_y], tol)) {
                        pts_on_xy_AR_planes[i_x][i_y].push_back( coords_gpt.first );
                    }
                }
                for (int i_z = 0; i_z < AR_planes[2].size(); i_z++) {
                    if (compareValues(coords_gpt.first[2], AR_planes[2][i_z], tol)) {
                        pts_on_xz_AR_planes[i_x][i_z].push_back( coords_gpt.first );
                    }
                }
            }
        }
        for (int i_y = 0; i_y < AR_planes[1].size(); i_y++) {
            if (compareValues(coords_gpt.first[1], AR_planes[1][i_y], tol)) {
                for (int i_z = 0; i_z < AR_planes[2].size(); i_z++) {
                    if (compareValues(coords_gpt.first[2], AR_planes[2][i_z], tol)) {
                        pts_on_yz_AR_planes[i_y][i_z].push_back( coords_gpt.first );
                    }
                }
            }
        }
    }

    void updateIsInside() {
        // To figure out if an AR corner point is inside the pore-space, we ray cast in the x,y,z directions and accept the majority vote of inside vs. outside 
        isInside.clear();
        int count = 0;
        for (int i_z = 0; i_z < AR_planes[2].size(); i_z++) {
            for (int i_y = 0; i_y < AR_planes[1].size(); i_y++) {
                for (int i_x = 0; i_x < AR_planes[0].size(); i_x++) {
                    std::vector<int> inside_outside_votes(2);

                    // First check if the corner point has been created into a gmsh point already. This would indicate that it is on the border.
                    if (std::get<1>(AR_corner_pts[count])) { isInside.push_back(2); count += 1; continue; }
                    
                    // If not, check the above and below points in each direction
                    std::vector<double> corner_coord = {AR_planes[0][i_x], AR_planes[1][i_y], AR_planes[2][i_z]};

                    int N_above = 0, N_below = 0;
                    for (int i = 0; i < pts_on_xy_AR_planes[i_x][i_y].size(); i++) { if (corner_coord[2] < pts_on_xy_AR_planes[i_x][i_y][i][2]) { N_above += 1; } else { N_below += 1; } }
                    if (N_above == 0 || N_below == 0) { isInside.push_back(0); count += 1; continue; }
                    if (N_above % 2 == 1) { inside_outside_votes[0] += 1; } else { inside_outside_votes[1] += 1; }
                    if (N_below % 2 == 1) { inside_outside_votes[0] += 1; } else { inside_outside_votes[1] += 1; }
                    
                    N_above = 0; N_below = 0;
                    for (int i = 0; i < pts_on_xz_AR_planes[i_x][i_z].size(); i++) { if (corner_coord[1] < pts_on_xz_AR_planes[i_x][i_z][i][1]) { N_above += 1; } else { N_below += 1; } }
                    if (N_above == 0 || N_below == 0) { isInside.push_back(0); count += 1; continue; }
                    if (N_above % 2 == 1) { inside_outside_votes[0] += 1; } else { inside_outside_votes[1] += 1; }
                    if (N_below % 2 == 1) { inside_outside_votes[0] += 1; } else { inside_outside_votes[1] += 1; }
                    
                    N_above = 0; N_below = 0;
                    for (int i = 0; i < pts_on_yz_AR_planes[i_y][i_z].size(); i++) { if (corner_coord[0] < pts_on_yz_AR_planes[i_y][i_z][i][0]) { N_above += 1; } else { N_below += 1; } }
                    if (N_above == 0 || N_below == 0) { isInside.push_back(0); count += 1; continue; }
                    if (N_above % 2 == 1) { inside_outside_votes[0] += 1; } else { inside_outside_votes[1] += 1; }
                    if (N_below % 2 == 1) { inside_outside_votes[0] += 1; } else { inside_outside_votes[1] += 1; }
                    
                    if (inside_outside_votes[0] == inside_outside_votes[1]) { 
                        cerr << "gmsh_utils.cpp: ARGridManager::updateIsInside: CRITICAL ERROR: Could not determine if AR corner point is inside, outside, or on the boundary of the domain. ";
                        cerr << "Point is on: {" << corner_coord[0] << ", " << corner_coord[1] << ", " << corner_coord[2] << "}. ";
                        cerr << "There were " << inside_outside_votes[0] << " inside votes and " << inside_outside_votes[1] << " outside votes." << endl;
                        exit(1); }
                    if (inside_outside_votes[0] == inside_outside_votes[1] + 1 || inside_outside_votes[0] + 1 == inside_outside_votes[1] + 1) {
                        cerr << "gmsh_utils.cpp: ARGridManager::updateIsInside: Warning: The ray-casting used to determining whether an AR corner point is inside, outside, or on the boundary of the domain may be inaccurate. Be sure to check the mesh before use." << endl;
                    }
                    if (inside_outside_votes[0] > inside_outside_votes[1]) { isInside.push_back(1); } else { isInside.push_back(0); }
                    count += 1;
                }
            }
        }

    }
};
*/



/*
void followCurveLoop(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns, const std::vector<std::vector<std::pair<int, int>>> &bndry_pts,
    std::vector<int> &curveLoop, const string &add_mode, bool verbose = false);
void followCurveLoop(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns, const std::vector<std::vector<std::pair<int, int>>> &bndry_pts,
    std::vector<int> &curveLoop, const string &add_mode, bool verbose)
{
    // Figure out the orientation of the other surface lines
    for (int i_ln = 0; i_ln < bndry_pts.size(); i_ln++) {
        for (int i_ln2 = 0; i_ln2 < unused_inds.size(); i_ln2++) {
                        
            // If the line was already added to curveLoop, skip it. Each line can only be used once right now
            //bool wasUsed = true;
            //int unused_ind;
            //for (int i_u = 0; i_u < unused_inds.size(); i_u++) { if (i_ln2 == unused_inds[i_u]) { unused_ind = i_u; wasUsed = false; break; } }
            //if (wasUsed) { continue; }
            
            //cout << "unused_inds" << endl;
            //for (int i = 0; i < unused_inds.size(); i++) { cout << unused_inds[i] << endl; }
 
            if (verbose) {
                cout << "gmsh_utils.cpp: followCurveLoop(): Considering line " << unused_inds[i_ln2];
                cout << " with points " << bndry_pts[unused_inds[i_ln2]][0].second << " " << bndry_pts[unused_inds[i_ln2]][1].second << endl;
                std::vector<double> pt_para_coords, dummy_coords;
                gmsh::model::getValue(0, bndry_pts[unused_inds[i_ln2]][0].second, pt_para_coords, dummy_coords);
                cout << "gmsh_utils.cpp: followCurveLoop(): They have endpoints: ";
                cout << "{" << dummy_coords[0] << ", " << dummy_coords[1] << ", " << dummy_coords[2] << "}    ";
                gmsh::model::getValue(0, bndry_pts[unused_inds[i_ln2]][1].second, pt_para_coords, dummy_coords);
                cout << "{" << dummy_coords[0] << ", " << dummy_coords[1] << ", " << dummy_coords[2] << "}" << endl;
            }
            
            // Determine if one of the end points is equivalent to the target point, and adjust/add to the variables accordingly
            if (bndry_pts[unused_inds[i_ln2]][0].second == target_pt) {
                if (verbose) { cout << "gmsh_utils.cpp: followCurveLoop(): First point matches target!" << endl; }
                if (add_mode == "back") { curveLoop.push_back( glns[unused_inds[i_ln2]] ); }
                else if (add_mode == "front") {
                    std::vector<int> dummy_curveLoop = {-glns[unused_inds[i_ln2]]};
                    dummy_curveLoop.insert(dummy_curveLoop.end(), curveLoop.begin(), curveLoop.end());
                    curveLoop = dummy_curveLoop;
                }
                
                target_pt = bndry_pts[unused_inds[i_ln2]][1].second;
                unused_inds.erase(unused_inds.begin() + i_ln2);
                if (verbose) { cout << "gmsh_utils.cpp: followCurveLoop(): New target point tag: " << target_pt << endl; }
                break;
            }
            if (bndry_pts[unused_inds[i_ln2]][1].second == target_pt) {
                if (verbose) { cout << "gmsh_utils.cpp: followCurveLoop(): Second point matches target!" << endl; }
                if (add_mode == "back") { curveLoop.push_back( -glns[unused_inds[i_ln2]] ); }
                else if (add_mode == "front") {
                    std::vector<int> dummy_curveLoop = {glns[unused_inds[i_ln2]]};
                    dummy_curveLoop.insert(dummy_curveLoop.end(), curveLoop.begin(), curveLoop.end());
                    curveLoop = dummy_curveLoop;
                }
                
                target_pt = bndry_pts[unused_inds[i_ln2]][0].second;
                unused_inds.erase(unused_inds.begin() + i_ln2);
                if (verbose) { cout << "gmsh_utils.cpp: followCurveLoop(): New target point tag: " << target_pt << endl; }
                break;
            }
        }
    }
}

void obtainFirstLineForCurveLoopFinder(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns,
    const std::vector<std::vector<std::pair<int, int>>> &bndry_pts, std::vector<int> &curveLoop, bool verbose = false);
void obtainFirstLineForCurveLoopFinder(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns,
    const std::vector<std::vector<std::pair<int, int>>> &bndry_pts, std::vector<int> &curveLoop, bool verbose)
{
    // If the line is available, use it as the starting line
    target_pt = bndry_pts[unused_inds[0]][1].second; // Get the pt tag of the [1] point of the [i_ln] line
    curveLoop.push_back( glns[unused_inds[0]] );
    unused_inds.erase(unused_inds.begin());
    if (verbose) {
        std::vector<double> pt_para_coords, dummy_coords;
        gmsh::model::getValue(0, target_pt, pt_para_coords, dummy_coords);
        cout << "gmsh_utils.cpp: obtainFirstLineForCurveLoopFinder(): Found starting target point: tag = " << target_pt << "    , coord =  {" << dummy_coords[0] << ", " << dummy_coords[1] << ", " << dummy_coords[2] << "}" << endl;
    }
}

int getCurveLoopEndPointCoords(const std::vector<int> &curveLoop, std::vector<std::vector<double>> &end_coords)
{
    // Clear and intialize the provided end coords variable
    end_coords.clear();
    for (int i = 0; i < 2; i++) { end_coords.push_back( std::vector<double>() ); }

    // Get the coordinates of the end points to the curve loop
    std::vector<double> pt_para_coords; std::vector<std::pair<int, int>> pt_tags;
    gmsh::model::getBoundary({{1, curveLoop[0]}}, pt_tags, false, false, false);
    int target_pt = pt_tags[0].second; // just in case we need to try add_mode = "front"
    gmsh::model::getValue(0, pt_tags[0].second, pt_para_coords, end_coords[0]);
    gmsh::model::getBoundary({{1, curveLoop[curveLoop.size() - 1]}}, pt_tags, false, false, false);
    gmsh::model::getValue(0, pt_tags[1].second, pt_para_coords, end_coords[1]);
    return target_pt;
}

void findCurveLoopsRecursive(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns,
    const std::vector<std::vector<std::pair<int, int>>> &bndry_pts, std::vector<std::vector<int>> &curveLoops, const double tol, bool verbose = false);
void findCurveLoopsRecursive(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns,
    const std::vector<std::vector<std::pair<int, int>>> &bndry_pts, std::vector<std::vector<int>> &curveLoops, const double tol, bool verbose)
{
    // Initiate a curve loop vector
    std::vector<int> curveLoop;

    // Obtain an initial line for the curveloop
    obtainFirstLineForCurveLoopFinder(target_pt, unused_inds, glns, bndry_pts, curveLoop, verbose);

    if (verbose) { cout << "gmsh_utils.cpp: findCurveLoopsRecursive(): There are " << glns.size() << " lines." << endl; }

    // Perform the initial attempt at finding the curve loop
    string add_mode = "back";
    followCurveLoop(target_pt, unused_inds, glns, bndry_pts, curveLoop, add_mode, verbose);
    
    if (verbose) { cout << "gmsh_utils.cpp: findCurveLoopsRecursive(): The curve loop has been followed to the end." << endl; }

    // Get the coordinates of the end points to the curve loop (this is to determine which of the following cases to pursue)
    std::vector<std::vector<double>> end_coords;
    target_pt = getCurveLoopEndPointCoords(curveLoop, end_coords);


    // If end points do not match, and not all lines were used, try going the other way.
    //      If end points do not match, and not all lines were used, redo process with remaining lines to see if there is a loop in there.
    //      If end points match, but not all lines were used, define the loop/surface and redo process with remaining lines to see if there is a loop in there
    //      If end points do not match, but all lines were used, there is no surface
    //      If end points match, and all lines were used, ERROR. This should not happen
    // If end points match, but not all lines were used, define the loop/surface and redo process with remaining lines to see if there is a loop in there
    // If end points do not match, but all lines were used, there is no surface
    // If end points match, and all lines were used, define the surface

    // If end points do not match, and not all lines were used, try going the other way.
    if (!compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() != 0) {
        if (verbose) { cout << "gmsh_utils.cpp: findCurveLoopsRecursive(): The end points do not match, and not all lines were used. Going to try the other way." << endl; }
        if (verbose) { cout << "gmsh_utils.cpp: findCurveLoopsRecursive(): New target point tag: " << target_pt << endl; }
        if (verbose) { cout << "gmsh_utils.cpp: findCurveLoopsRecursive(): There are " << unused_inds.size() << " lines left." << endl; }
        
        string add_mode = "front";
        followCurveLoop(target_pt, unused_inds, glns, bndry_pts, curveLoop, add_mode, verbose);
        target_pt = getCurveLoopEndPointCoords(curveLoop, end_coords);
        
        // If end points do not match, and not all lines were used, redo process with remaining lines to see if there is a loop in there
        if (!compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() != 0) {
            findCurveLoopsRecursive(target_pt, unused_inds, glns, bndry_pts, curveLoops, tol, verbose);
            return;
        }

        // If end points match, but not all lines were used, define the loop/surface and redo process with remaining lines to see if there is a loop in there
        if (compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() != 0) {
            curveLoops.push_back( curveLoop );
            findCurveLoopsRecursive(target_pt, unused_inds, glns, bndry_pts, curveLoops, tol, verbose);
            return;
        }
    
        // If end points do not match, but all lines were used, there is no surface
        if (!compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() != 0) { return; }
    
        // If end points match, and all lines were used, ERROR. This should not happen
        if (compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() != 0) {
            cerr << "gmsh_utils.cpp: findCurveLoopsRecursive(): CRITICAL ERROR: End points match and all lines are used after going through the lines backwards. This should never happen." << endl;
            exit(1);
        }
    }


    // If end points match, but not all lines were used, define the loop/surface and redo process with remaining lines to see if there is a loop in there
    if (compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() != 0) {
        curveLoops.push_back( curveLoop );
        findCurveLoopsRecursive(target_pt, unused_inds, glns, bndry_pts, curveLoops, tol, verbose);
        return;
    }

    // If end points do not match, but all lines were used, there is no surface
    if (!compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() == 0) { return; }

    // If end point coordinates match, and all lines were used, define the surface
    if (compareVectors(end_coords[0], end_coords[1], tol) && unused_inds.size() == 0) {
        curveLoops.push_back( curveLoop );
        return;
    }

    // One of the previous options should have been chosen
    cerr << "gmsh_utils.cpp: findCurveLoopsRecursive(): CRITICAL ERROR: No option chosen. This should never happen." << endl;
    exit(1);
}

std::vector<std::vector<int>> findCurveLoops(std::vector<int> &glns, double tol, bool verbose = false);
std::vector<std::vector<int>> findCurveLoops(std::vector<int> &glns, double tol, bool verbose)
{
    // From the gmsh line tags (glns), get the gmsh point tags of the end points in gmsh format
    std::vector<std::vector<std::pair<int, int>>> bndry_pts;
    for (int i_ln = 0; i_ln < glns.size(); i_ln++) {
        bndry_pts.push_back( std::vector<std::pair<int, int>>() );
        gmsh::model::getBoundary({{1, glns[i_ln]}}, bndry_pts[bndry_pts.size() - 1], false, false, false);
    }

    // Define a vector to store the unused point indices, a vector for the curve loop, and an integer for the target point
    std::vector<int> unused_inds; for (int i = 0; i < bndry_pts.size(); i++) { unused_inds.push_back(i); }
    std::vector<std::vector<int>> curveLoops;
    int target_pt;
    
    // Start the recursive curve-loop-finding process 
    findCurveLoopsRecursive(target_pt, unused_inds, glns, bndry_pts, curveLoops, tol, verbose);

    return curveLoops;
}
*/
/*
void followSurfaceLoop(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<int> &surfaceLoop, bool verbose = false);
void followSurfaceLoop(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<int> &surfaceLoop, bool verbose)
{
    for (int i_s = 0; i_s < bndry_glns.size(); i_s++) {
        for (int i_s2 = 0; i_s2 < unused_inds.size(); i_s2++) {
            // For debugging
            if (verbose) {
                cout << "gmsh_utils.cpp: followSurfaceLoop(): Considering surface index " << unused_inds[i_s2] << " with lines";
                for (int i = 0; i < bndry_glns[unused_inds[i_s2]].size(); i++) { cout << " " << bndry_glns[unused_inds[i_s2]][i]; }
                cout << endl;
            }
            
            // Determine if one of the lines in surface unused_inds[i_s2] is equivalent to a target line. Do this by
            // creating and sorting a list with the target glns and the glns from surface unused_inds[i_s2]
            std::vector<int> glns_list = bndry_glns[unused_inds[i_s2]];
            glns_list.insert(glns_list.end(), target_lns.begin(), target_lns.end());
            std::sort(glns_list.begin(), glns_list.end(), [](const int &a, const int &b){ return a < b; });
            
            // For debugging
            if (verbose) { cout << "Gmsh line tags in combined list:"; for (int i = 0; i < glns_list.size(); i++) { cout << " " << glns_list[i]; } cout << endl; }
                
            // Now, go through glns_list and look for doubles. If a double is found, the surface can be removed from
            // glns_list and added to the surface loop. Continue through the whole list to remove all doubles.
            // There should never be triples of a gln in glns_list
            bool skip = false;
            std::vector<int> glns_list2;
            for (int i_db = 0; i_db < glns_list.size(); i_db++) {
                if (skip) { if (i_db < glns_list.size() - 1) { assert (glns_list[i_db] != glns_list[i_db + 1]); } skip = false; continue; }
                if (i_db == glns_list.size() - 1) { glns_list2.push_back( glns_list[i_db] ); continue; }
                if (glns_list[i_db] == glns_list[i_db + 1]) { skip = true; continue; }
                glns_list2.push_back( glns_list[i_db] );
            }

            // If any lines matched, the beginning and ending sizes will be different. Add the surface and adjust
            // the other variables if the surface is added
            if (glns_list2.size() != glns_list.size()) {
                if (verbose) { cout << "gmsh_utils.cpp: followSurfaceLoop(): Some lines matched the target lines!" << endl; }
                
                surfaceLoop.push_back( gsurfs[unused_inds[i_s2]] );
                target_lns = glns_list2;
                unused_inds.erase(unused_inds.begin() + i_s2);
                
                if (verbose) {
                    cout << "gmsh_utils.cpp: followSurfaceLoop(): New target line tags:";
                    for (int i = 0; i < target_lns.size(); i++) { cout << " " << target_lns[i]; }
                    cout << endl;
                }
                
                break;
            }
        }
    }
}

void obtainFirstSurfaceForSurfaceLoopFinder(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<int> &surfaceLoop, bool verbose = false);
void obtainFirstSurfaceForSurfaceLoopFinder(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<int> &surfaceLoop, bool verbose)
{
    // Use the first available surface to initiate the surface loop, and use its lines as the target lines
    target_lns = bndry_glns[unused_inds[0]];
    surfaceLoop.push_back( gsurfs[unused_inds[0]] );
    unused_inds.erase(unused_inds.begin());
    if (verbose) {
        cout << "gmsh_utils.cpp: obtainFirstSurfaceForSurfaceLoopFinder(): Found starting surface " << surfaceLoop[surfaceLoop.size() - 1] << " with line tags ";
        for (int i = 0; i < target_lns.size(); i++) { cout << " " << target_lns[i]; }
        cout << endl;
    }
}

void findSurfaceLoopsRecursive(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<std::vector<int>> &surfaceLoops, bool verbose = false);
void findSurfaceLoopsRecursive(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<std::vector<int>> &surfaceLoops, bool verbose)
{
    // Initiate a surface loop vector
    std::vector<int> surfaceLoop;

    if (verbose) { cout << "gmsh_utils.cpp: findSurfaceLoopsRecursive(): Starting finding process. There are " << unused_inds.size() << " surfaces." << endl; }

    // Obtain an initial surface for the surface loop
    obtainFirstSurfaceForSurfaceLoopFinder(target_lns, unused_inds, gsurfs, bndry_glns, surfaceLoop, verbose);

    // Perform the initial attempt at finding the surface loop
    followSurfaceLoop(target_lns, unused_inds, gsurfs, bndry_glns, surfaceLoop, verbose);
    
    if (verbose) { cout << "gmsh_utils.cpp: findSurfaceLoopsRecursive(): The surface loop has been followed to the end." << endl; }

    // After following the surface loop to the end, there should not be any target lines left, regardless of how many
    // surface loops there are in the provided surfaces
    assert (target_lns.size() == 0);
    
    // Add the loop to the surface loops
    surfaceLoops.push_back( surfaceLoop );

    // If all surfaces were used, exit.
    if (unused_inds.size() == 0) { return; }
    // If not all surfaces were used, call findSurfaceLoopsRecursive again.
    else {
        if (verbose) {
            cout << "gmsh_utils.cpp: findSurfaceLoopsRecursive(): Not all surfaces were used. ";
            cout << "There are " << unused_inds.size() << " surfaces left. ";
            cout << "Beginning search for additional surface loops." << endl; }
        findSurfaceLoopsRecursive(target_lns, unused_inds, gsurfs, bndry_glns, surfaceLoops, verbose);
        return;
    }
}

std::vector<std::vector<int>> findSurfaceLoops(const std::vector<int> &gsurfs, bool verbose = false);
std::vector<std::vector<int>> findSurfaceLoops(const std::vector<int> &gsurfs, bool verbose)
{
    // From the gmsh surface tags (gsurfs), get the gmsh line tags of the lines bounding each surface
    std::vector<std::vector<int>> bndry_glns;
    for (int i_surf = 0; i_surf < gsurfs.size(); i_surf++) {
        std::vector<std::pair<int, int>> glns_GF;
        gmsh::model::getBoundary({{2, gsurfs[i_surf]}}, glns_GF, false, false, false);
        bndry_glns.push_back( std::vector<int>() );
        for (int i_ln = 0; i_ln < glns_GF.size(); i_ln++) { bndry_glns[bndry_glns.size() - 1].push_back( glns_GF[i_ln].second ); }
    }

    // Define a vector to store the unused surface indices, a vector for the surface loop, and a vector for the target lines
    std::vector<int> unused_inds; for (int i = 0; i < gsurfs.size(); i++) { unused_inds.push_back(i); }
    std::vector<std::vector<int>> surfaceLoops;
    std::vector<int> target_lns;

    // Start the recursive surface-loop-finding process 
    findSurfaceLoopsRecursive(target_lns, unused_inds, gsurfs, bndry_glns, surfaceLoops, verbose);

    return surfaceLoops;
}




class LineObj
{
public:
    double tol;
    const std::pair<std::vector<double>, std::vector<double>> &pt_coords;
    
    int plane_dim = -1, plane_val;
    std::vector<int> alt_dim;

    bool isVertical = false, isHorizontal = false;
    double vert_val, hori_val;
    double m, b; // the slope and "y-intercept" for y = m*x + b

    bool returned_half = false;

    
    // Class constructors
    LineObj(const std::pair<std::vector<double>, std::vector<double>> &pt_coords_, int plane_dim_, const double tol_) : pt_coords(pt_coords_), plane_dim(plane_dim_), tol(tol_)
    {
        assert (compareValues(pt_coords.first[plane_dim], pt_coords.second[plane_dim], tol));
        plane_val = pt_coords.first[plane_dim];

        if (plane_dim == 0) { alt_dim = {1, 2}; }
        if (plane_dim == 1) { alt_dim = {0, 2}; }
        if (plane_dim == 2) { alt_dim = {0, 1}; }

        // Check if the line is vertical (w.r.t. the chosen alt_dim). If so, store the alt_dim[0] value. If not, compute the slope and "y-intercept" for y = m*x + b
        if (compareValues(pt_coords.first[alt_dim[0]], pt_coords.second[alt_dim[0]], tol)) { isVertical = true; vert_val = pt_coords.first[alt_dim[0]]; }
        else {
            m = (pt_coords.second[alt_dim[1]] - pt_coords.first[alt_dim[1]]) / (pt_coords.second[alt_dim[0]] - pt_coords.first[alt_dim[0]]);
            b = -m*pt_coords.first[alt_dim[0]] + pt_coords.first[alt_dim[1]];
        }
        if (compareValues(pt_coords.first[alt_dim[1]], pt_coords.second[alt_dim[1]], tol)) { isHorizontal = true; hori_val = pt_coords.first[alt_dim[1]]; }
    }

    int isSegmentIntersected(const std::pair<std::vector<double>, std::vector<double>> &in_pts, const string check_dir = "pos") {
        // Make sure the input line points are on the same plane
        assert (compareValues(in_pts.first[plane_dim], plane_val, tol) && compareValues(in_pts.second[plane_dim], plane_val, tol));
        
        // Initiate variables
        double alt_dim_0_int, alt_dim_1_int;

        cout << "m = " << m << ", b = " << b << endl;
        
        //bool inputIsHorizontal = false;
        //if (compareValues(in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]], tol)) { inputIsHorizontal = true; }
        
        // Check if the input line is vertical (w.r.t. alt_dim[0])
        bool inputIsVertical = false;
        bool inputIsHorizontal = false;
        if (compareValues(in_pts.first[alt_dim[0]], in_pts.second[alt_dim[0]], tol)) {
            inputIsVertical = true;
            // If the class line is also vertical...
            if (isVertical) {
                // ... see if in_pts.first[alt_dim[0]] == vert_val. This means lines overlap and we treat the analysis as "inconclusive"
                if (compareValues(in_pts.first[alt_dim[0]], vert_val, tol)) { return -1; }
                // If in_pts.first[alt_dim[0]] != vert_val, the lines do not intersect.
                return 0; // Let's treat both as 0; the lines attached to the input line will together return 1 due to returned_half.
            }

            // If the input line is vertical, but the class line is not, compute the intersection using the class line slope and "y-intercept"
            alt_dim_0_int = in_pts.first[alt_dim[0]];
            if (isHorizontal) { alt_dim_1_int = hori_val; }
            else { alt_dim_1_int = m*alt_dim_0_int + b; }

            cout << "INTERSECTION COORDS: alt_dim_0_int = " << alt_dim_0_int << ", alt_dim_1_int = " << alt_dim_1_int << endl;

            if (check_dir == "pos") {
                // If alt_dim_0_int is not on the positive side of the class point, it does not intersect on the correct side
                if (alt_dim_0_int < std::max({pt_coords.first[alt_dim[0]], pt_coords.second[alt_dim[0]]})) { return 0; }
                double large_1 = std::max({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                double small_1 = std::min({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                // If alt_dim_1_int is smaller than the larger of the alt_dim[1] values of in_pts, but larger than the smaller one, then it intersects the segment
                if (alt_dim_1_int < large_1 && alt_dim_1_int > small_1) { return 1; }
                // If alt_dim_1_int is larger than the larger alt_dim[1] values of in_pts, or smaller than the smaller one, then it does not intersect the segment
                if (alt_dim_1_int > large_1 || alt_dim_1_int < small_1) { return 0; }
                // If alt_dim_1_int is equal to the "y" values of the input points, then it intersects at the input points, and the analysis is inconclusive
                if (compareValues(alt_dim_1_int, in_pts.first[alt_dim[1]], tol) || compareValues(alt_dim_1_int, in_pts.second[alt_dim[1]], tol)) { return -1; }
                cerr << "gmsh_utils.cpp: LineObj().isSegmentIntersected(): CRITICAL ERROR: Cannot be determined whether input line is intersected by class line." << endl;
                exit(1);
            }
            if (check_dir == "neg") {
                // If alt_dim_0_int is not on the negative side of the class point, it does not intersect on the correct side
                if (alt_dim_0_int > std::min({pt_coords.first[alt_dim[0]], pt_coords.second[alt_dim[0]]})) { return 0; }
                double large_1 = std::max({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                double small_1 = std::min({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                // If alt_dim_1_int is smaller than the larger of the alt_dim[1] values of in_pts, but larger than the smaller one, then it intersects the segment
                if (alt_dim_1_int < large_1 && alt_dim_1_int > small_1) { return 1; }
                // If alt_dim_1_int is larger than the larger alt_dim[1] values of in_pts, or smaller than the smaller one, then it does not intersect the segment
                if (alt_dim_1_int > large_1 || alt_dim_1_int < small_1) { return 0; }
                // If alt_dim_1_int is equal to the "y" values of the input points, then it intersects at the input points, and the analysis is inconclusive
                if (compareValues(alt_dim_1_int, in_pts.first[alt_dim[1]], tol) || compareValues(alt_dim_1_int, in_pts.second[alt_dim[1]], tol)) { return -1; }
                cerr << "gmsh_utils.cpp: LineObj().isSegmentIntersected(): CRITICAL ERROR: Cannot be determined whether input line is intersected by class line." << endl;
                exit(1);
            }
        }
        // Check if the input line is horizontal (w.r.t. alt_dim[1])
        else if (compareValues(in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]], tol)) {
            inputIsHorizontal = true;
            // If the class line is also horizontal...
            if (isHorizontal) {
                // ... see if in_pts.first[alt_dim[1]] == hori_val. This means lines overlap and we treat the analysis as "inconclusive"
                if (compareValues(in_pts.first[alt_dim[1]], hori_val, tol)) { return -1; }
                // If in_pts.first[alt_dim[1]] != hori_val, the lines do not intersect.
                return 0;
            }

            // If the input line is horizontal, but the class line is not, compute the intersection using the class line slope and "y-intercept"
            alt_dim_1_int = in_pts.first[alt_dim[1]];
            if (isVertical) { alt_dim_0_int = vert_val; }
            else { alt_dim_0_int = (alt_dim_1_int - b) / m; }

            cout << "INTERSECTION COORDS: alt_dim_0_int = " << alt_dim_0_int << ", alt_dim_1_int = " << alt_dim_1_int << endl;

            if (check_dir == "pos") {
                // If alt_dim_0_int is not on the positive side of the class point, it does not intersect on the correct side
                if (alt_dim_0_int < std::max({pt_coords.first[alt_dim[0]], pt_coords.second[alt_dim[0]]})) { return 0; }
                double large_0 = std::max({in_pts.first[alt_dim[0]], in_pts.second[alt_dim[0]]});
                double small_0 = std::min({in_pts.first[alt_dim[0]], in_pts.second[alt_dim[0]]});
                // If alt_dim_0_int is smaller than the larger of the alt_dim[0] values of in_pts, but larger than the smaller one, then it intersects the segment
                if (alt_dim_0_int < large_0 && alt_dim_0_int > small_0) { return 1; }
                // If alt_dim_0_int is larger than the larger alt_dim[0] values of in_pts, or smaller than the smaller one, then it does not intersect the segment
                if (alt_dim_0_int > large_0 || alt_dim_0_int < small_0) { return 0; }
                // If alt_dim_0_int is equal to the "x" values of the input points, then it intersects at the input points, and the analysis is inconclusive
                if (compareValues(alt_dim_0_int, in_pts.first[alt_dim[0]], tol) || compareValues(alt_dim_0_int, in_pts.second[alt_dim[0]], tol)) { return -1; }
                cerr << "gmsh_utils.cpp: LineObj().isSegmentIntersected(): CRITICAL ERROR: Cannot be determined whether input line is intersected by class line." << endl;
                exit(1);
            }
            if (check_dir == "neg") {
                // If alt_dim_0_int is not on the negative side of the class point, it does not intersect on the correct side
                if (alt_dim_0_int > std::min({pt_coords.first[alt_dim[0]], pt_coords.second[alt_dim[0]]})) { return 0; }
                double large_0 = std::max({in_pts.first[alt_dim[0]], in_pts.second[alt_dim[0]]});
                double small_0 = std::min({in_pts.first[alt_dim[0]], in_pts.second[alt_dim[0]]});
                // If alt_dim_0_int is smaller than the larger of the alt_dim[0] values of in_pts, but larger than the smaller one, then it intersects the segment
                if (alt_dim_0_int < large_0 && alt_dim_0_int > small_0) { return 1; }
                // If alt_dim_0_int is larger than the larger alt_dim[0] values of in_pts, or smaller than the smaller one, then it does not intersect the segment
                if (alt_dim_0_int > large_0 || alt_dim_0_int < small_0) { return 0; }
                // If alt_dim_0_int is equal to the "x" values of the input points, then it intersects at the input points, and the analysis is inconclusive
                if (compareValues(alt_dim_0_int, in_pts.first[alt_dim[0]], tol) || compareValues(alt_dim_0_int, in_pts.second[alt_dim[0]], tol)) { return -1; }
                cerr << "gmsh_utils.cpp: LineObj().isSegmentIntersected(): CRITICAL ERROR: Cannot be determined whether input line is intersected by class line." << endl;
                exit(1);
            }
        }
        else {
            // If the input line is not vertical, compute the slope and "y-intercept" for y = m*x + b for the input line
            double m2 = (in_pts.second[alt_dim[1]] - in_pts.first[alt_dim[1]]) / (in_pts.second[alt_dim[0]] - in_pts.first[alt_dim[0]]);
            double b2 = -m2*in_pts.first[alt_dim[0]] + in_pts.first[alt_dim[1]];

            cout << "m2 = " << m2 << ", b2 = " << b2 << endl;
        
            // If the class line is vertical, compute the intersection using the input line slope and "y-intercept"
            if (isVertical) { alt_dim_0_int = vert_val; alt_dim_1_int = m2*alt_dim_0_int + b2; }
            else if (isHorizontal) { alt_dim_1_int = hori_val ; alt_dim_0_int = (alt_dim_1_int - b2) / m2; }
            // If the class line is not vertical either, compute the intersecting "x" and "y" (could be other dimensions; x and y are used for descriptive purposes) by using the y = m*x + b equations  
            else { alt_dim_0_int = (b2 - b) / (m - m2); alt_dim_1_int = m2*alt_dim_0_int + b2; }
            
            cout << "INTERSECTION COORDS: alt_dim_0_int = " << alt_dim_0_int << ", alt_dim_1_int = " << alt_dim_1_int << endl;

            if (check_dir == "pos") {
                // If alt_dim_0_int is not on the positive side of the class point, it does not intersect on the correct side
                if (alt_dim_0_int < std::max({pt_coords.first[alt_dim[0]], pt_coords.second[alt_dim[0]]})) { return 0; }
                double large_1 = std::max({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                double small_1 = std::min({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                // If alt_dim_1_int is smaller than the larger of the alt_dim[1] values of in_pts, but larger than the smaller one, then it intersects the segment
                if (alt_dim_1_int < large_1 && alt_dim_1_int > small_1) { return 1; }
                // If alt_dim_1_int is larger than the larger alt_dim[1] values of in_pts, or smaller than the smaller one, then it does not intersect the segment
                if (alt_dim_1_int > large_1 || alt_dim_1_int < small_1) { return 0; }
                // If alt_dim_1_int is equal to the "y" values of the input points, then it intersects at the input points, and the analysis is inconclusive
                if (compareValues(alt_dim_1_int, in_pts.first[alt_dim[1]], tol) || compareValues(alt_dim_1_int, in_pts.second[alt_dim[1]], tol)) { return -1; }
                cerr << "gmsh_utils.cpp: LineObj().isSegmentIntersected(): CRITICAL ERROR: Cannot be determined whether input line is intersected by class line." << endl;
                exit(1);
            }
            if (check_dir == "neg") {
                // If alt_dim_0_int is not on the negative side of the class point, it does not intersect on the correct side
                if (alt_dim_0_int > std::min({pt_coords.first[alt_dim[0]], pt_coords.second[alt_dim[0]]})) { return 0; }
                double large_1 = std::max({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                double small_1 = std::min({in_pts.first[alt_dim[1]], in_pts.second[alt_dim[1]]});
                // If alt_dim_1_int is smaller than the larger of the alt_dim[1] values of in_pts, but larger than the smaller one, then it intersects the segment
                if (alt_dim_1_int < large_1 && alt_dim_1_int > small_1) { return 1; }
                // If alt_dim_1_int is larger than the larger alt_dim[1] values of in_pts, or smaller than the smaller one, then it does not intersect the segment
                if (alt_dim_1_int > large_1 || alt_dim_1_int < small_1) { return 0; }
                // If alt_dim_1_int is equal to the "y" values of the input points, then it intersects at the input points, and the analysis is inconclusive
                if (compareValues(alt_dim_1_int, in_pts.first[alt_dim[1]], tol) || compareValues(alt_dim_1_int, in_pts.second[alt_dim[1]], tol)) { return -1; }
                cerr << "gmsh_utils.cpp: LineObj().isSegmentIntersected(): CRITICAL ERROR: Cannot be determined whether input line is intersected by class line." << endl;
                exit(1);
            }
        }
        cerr << "gmsh_utils.cpp: LineObj().isSegmentIntersected(): CRITICAL ERROR: Should never get here." << endl;
        exit(1);
        return -1;
    }

    int computeCurveLoopIntersections(const std::vector<std::pair<std::vector<double>, std::vector<double>>> &in_cl, const string check_dir = "pos") {
        // Find how many times the line of curve loop i_cl intersects curve loop i_cl2
        int intersection, total_intersections = 0;
        //cout << "in_cl.size(): " << in_cl.size() << endl;
        for (int i_ln2 = 0; i_ln2 < in_cl.size(); i_ln2++) {
            
            cout << "i_ln2: " << i_ln2 << endl;
            cout << "pt_coords.first:  " << pt_coords.first[0] << " " << pt_coords.first[1] << " " << pt_coords.first[2] << endl;
            cout << "pt_coords.second:  " << pt_coords.second[0] << " " << pt_coords.second[1] << " " << pt_coords.second[2] << endl;
            cout << "in_pts.first:  " << in_cl[i_ln2].first[0] << " " << in_cl[i_ln2].first[1] << " " << in_cl[i_ln2].first[2] << endl;
            cout << "in_pts.second:  " << in_cl[i_ln2].second[0] << " " << in_cl[i_ln2].second[1] << " " << in_cl[i_ln2].second[2] << endl;
            
            intersection = isSegmentIntersected(in_cl[i_ln2], check_dir);
            
            cout << "intersection: " << intersection << endl;
            
            if (intersection == -1) { return -1; }
            else { total_intersections += intersection; }
            //cout << ", total_intersections: " << total_intersections;
            //if (returned_half) { cout << ", returned_half: true" << endl; }
            //else { cout << ", returned_half: false" << endl; }
        }
        // If returned_half is true at the end of the previous loop, it means the class line is likely tangential to the curve loop.
        // In this case, return an inconclusive result
        if (returned_half) { return -1; }
        return total_intersections;
    }
};



std::vector<std::vector<std::vector<int>>> arrangeCurveLoopsForSurfaces(const std::vector<std::vector<int>> &curveLoops, int dim, const double tol, bool verbose = false);
std::vector<std::vector<std::vector<int>>> arrangeCurveLoopsForSurfaces(const std::vector<std::vector<int>> &curveLoops, int dim, const double tol, bool verbose)
{
    int N_cl = curveLoops.size();

    // Get the coordinates of the end points of each line in the curve loops
    std::vector<std::vector<std::pair<std::vector<double>, std::vector<double>>>> curveLoops_coords;
    for (int i_cl = 0; i_cl < curveLoops.size(); i_cl++) {
        curveLoops_coords.push_back( std::vector<std::pair<std::vector<double>, std::vector<double>>>() );
        // Get the unique gmsh point tags of the lines in the curveloop
        std::vector<int> glns;
        for (int i_ln = 0; i_ln < curveLoops[i_cl].size(); i_ln++) {
            int gln = curveLoops[i_cl][i_ln]; if (gln < 0) { gln *= -1.0; }
            glns.push_back( gln ); }
        removeDuplicates(glns, [](const int &a, const int &b){ return a < b; });
        
        // Get the coordinates of the points in the curve loop
        for (int i_ln = 0; i_ln < glns.size(); i_ln++) {
            curveLoops_coords[i_cl].push_back( std::pair<std::vector<double>, std::vector<double>>() );
            // Get the gmsh point tags of the gmsh line end points
            std::vector<std::pair<int, int>> gpts_endpoints;
            gmsh::model::getBoundary({{1, glns[i_ln]}}, gpts_endpoints, false, false, false);
            // Get the coordinates of the end points
            std::vector<double> pt_para_coords;
            gmsh::model::getValue(0, gpts_endpoints[0].second, pt_para_coords, curveLoops_coords[i_cl][i_ln].first);
            gmsh::model::getValue(0, gpts_endpoints[1].second, pt_para_coords, curveLoops_coords[i_cl][i_ln].second);
        }
    }

    std::vector<std::vector<int>> curveLoops_groups; for (int i = 0; i < N_cl; i++) { curveLoops_groups.push_back( {i} ); }
    
    for (int i_cl = 0; i_cl < N_cl; i_cl++) {
        for (int i_cl2 = 0; i_cl2 < N_cl; i_cl2++) {
            if (i_cl2 == i_cl) { continue; }

            bool isInside = false, foundAnswer = false;
            std::vector<int> votesInsideOutside(2);
            int majorityDiff = 3;
            for (int i_ln = 0; i_ln < curveLoops_coords[i_cl].size(); i_ln++) {
                cout << "i_ln: " << i_ln << endl;
                // Make a line object with the line's end points
                LineObj line(curveLoops_coords[i_cl][i_ln], dim, tol);
                
                // Find how many times the line of curve loop i_cl intersects curve loop i_cl2. First check the positive direction, then the negative
                std::vector<int> total_intersections_pos_neg(2);
                total_intersections_pos_neg[0] = line.computeCurveLoopIntersections(curveLoops_coords[i_cl2], "pos");
                // Easiest case is if there are no intersections; this means curve loop i_cl is definitely not in i_cl2
                if (total_intersections_pos_neg[0] == 0) { isInside = false; foundAnswer = true; break; }
                // Check the negative direction
                total_intersections_pos_neg[1] = line.computeCurveLoopIntersections(curveLoops_coords[i_cl2], "neg");
                // Easiest case is if there are no intersections; this means curve loop i_cl is definitely not in i_cl2
                if (total_intersections_pos_neg[1] == 0) { isInside = false; foundAnswer = true; break; }
                // Because of numerical error, we let the intersection test come back inconclusive. In the case that both pos and neg come back inconclusive, move to the next i_ln 
                if (total_intersections_pos_neg[0] == -1 && total_intersections_pos_neg[1] == -1) { continue; }
                for (int i_tipn = 0; i_tipn < total_intersections_pos_neg.size(); i_tipn++) {
                    if (total_intersections_pos_neg[i_tipn] != -1) {
                        // Add vote according to total_intersections
                        if (total_intersections_pos_neg[i_tipn] % 2 == 0) { votesInsideOutside[1] += 1; } else { votesInsideOutside[0] += 1; }
                    }
                }
                
                // If one criteria receives a majority (set by majorityDiff), then exit with the answer
                cout << "VOTES IN AND OUT: " << votesInsideOutside[0] << ", " << votesInsideOutside[1] << endl;
                if (votesInsideOutside[0] >= votesInsideOutside[1] + majorityDiff) { isInside = true; foundAnswer = true; break; }
                if (votesInsideOutside[1] >= votesInsideOutside[0] + majorityDiff) { isInside = false; foundAnswer = true; break; }
            }
            if (!foundAnswer) { cerr << "gmsh_utils.cpp: arrangeCurveLoopsForSurfaces(): CRITICAL ERROR: Could not figure out if curve loop " << i_cl << " is inside or outside curve loop " << i_cl2 << "." << endl; exit(1); }
            
            if (isInside) {
                // Go through the groups and look for i_cl. Make sure its group does not have any one else in it. Then, copy that group to another variable
                std::vector<int> i_cl_group;
                int i_cl_ind = -1;
                for (int i_clg = 0; i_clg < curveLoops_groups.size(); i_clg++) {
                    for (int i_clg2 = 0; i_clg2 < curveLoops_groups[i_clg].size(); i_clg2++) {
                        if (curveLoops_groups[i_clg][i_clg2] == i_cl) {
                            assert (curveLoops_groups[i_clg].size() == 1);
                            i_cl_group = curveLoops_groups[i_clg];
                            i_cl_ind = i_clg;
                            break;
                        }
                    }
                    if (i_cl_ind != -1) { break; };
                }
                if (i_cl_ind == -1) { cerr << "gmsh_utils.cpp: arrangeCurveLoopsForSurfaces(): CRITICAL ERROR: Could not find an i_cl group in curveLoops_groups." << endl; exit(1); }
                
                // Erase the curve loop group of i_cl from curveLoops_groups
                curveLoops_groups.erase(curveLoops_groups.begin() + i_cl_ind);
                
                // Add the curve loop group of i_cl to the group of i_cl2 in curveLoops_groups. Make sure the previous size of i_cl2's group was 1
                bool flag = false;
                for (int i_clg = 0; i_clg < curveLoops_groups.size(); i_clg++) {
                    for (int i_clg2 = 0; i_clg2 < curveLoops_groups[i_clg].size(); i_clg2++) {
                        if (curveLoops_groups[i_clg][i_clg2] == i_cl2) {
                            assert (curveLoops_groups[i_clg].size() == 1);
                            curveLoops_groups[i_clg].insert(curveLoops_groups[i_clg].end(), i_cl_group.begin(), i_cl_group.end());
                            flag = true;
                            break;
                        }
                    }
                    if (flag) { break; };
                }
                if (!flag) { cerr << "gmsh_utils.cpp: arrangeCurveLoopsForSurfaces(): CRITICAL ERROR: Could not find curve loop group for i_cl2 in curveLoops_groups." << endl; exit(1); }
                break;
            }
        }
    }

    std::vector<std::vector<std::vector<int>>> curveLoops_arranged;
    for (int i_clg = 0; i_clg < curveLoops_groups.size(); i_clg++) {
        curveLoops_arranged.push_back( std::vector<std::vector<int>>() );
        for (int i_clg2 = 0; i_clg2 < curveLoops_groups[i_clg].size(); i_clg2++) {
            curveLoops_arranged[i_clg].push_back( curveLoops[curveLoops_groups[i_clg][i_clg2]] );
        }
    }

    return curveLoops_arranged;
}

























std::vector<std::pair<int, int>> importAndCutSTLVolume(const string &STL_file_path, const double &geo_scale, const std::vector<std::vector<double>> &AR_planes)
{
    gmsh::option::setNumber("Geometry.OCCAutoFix", 0); // We set this to zero, because if we do not, addSurfaceLoop may create duplicate surfaces.

    // Define the function name
    string FILE_NAME = "gmsh_utils.cpp";
    string FUNCTION_NAME = "importAndCutSTLVolume()";
    double tol = 1e-10;

    // Get the number of AR volumes defined according to the AR cut planes provided
    int N_AR = 1;
    for (int i_d = 0; i_d < AR_planes.size(); i_d++) { N_AR *= AR_planes[i_d].size() - 1; }
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Detected " << N_AR << " ARs defined by the AR cut planes provided." << endl;
    cout << "    NOTE: this is only the INITIAL number of ARs. The final number will depend on the geometry and if ARs are combined to minimize tiny AR." << endl;
    
    // Get the surfaces for each AR that are going to be analyzed
    //std::vector<std::vector<std::pair<int, double>>> AR_bounding_surfs; // Indices are [AR][surface], and then .first is the dim, .second is the value
    std::vector<std::vector<BoundingSurface>> AR_bounding_surfs; // Indices are [AR][surface], and then .first is the dim, .second is the value
    for (int i_z = 0; i_z < AR_planes[2].size() - 1; i_z++) {
        for (int i_y = 0; i_y < AR_planes[1].size() - 1; i_y++) {
            for (int i_x = 0; i_x < AR_planes[0].size() - 1; i_x++) {
                AR_bounding_surfs.push_back( std::vector<BoundingSurface>() );
                BoundingSurface surf(0, AR_planes[0][i_x], {1, 2}, {{AR_planes[1][i_y], AR_planes[1][i_y + 1]}, {AR_planes[2][i_z], AR_planes[2][i_z + 1]}});
                AR_bounding_surfs[AR_bounding_surfs.size() - 1].push_back( surf );
                surf.Set(1, AR_planes[1][i_y], {0, 2}, {{AR_planes[0][i_x], AR_planes[0][i_x + 1]}, {AR_planes[2][i_z], AR_planes[2][i_z + 1]}});
                AR_bounding_surfs[AR_bounding_surfs.size() - 1].push_back( surf );
                surf.Set(2, AR_planes[2][i_z], {0, 1}, {{AR_planes[0][i_x], AR_planes[0][i_x + 1]}, {AR_planes[1][i_y], AR_planes[1][i_y + 1]}});
                AR_bounding_surfs[AR_bounding_surfs.size() - 1].push_back( surf );

                if (i_x == AR_planes[0].size() - 2) {
                    surf.Set(0, AR_planes[0][i_x + 1], {1, 2}, {{AR_planes[1][i_y], AR_planes[1][i_y + 1]}, {AR_planes[2][i_z], AR_planes[2][i_z + 1]}});
                    AR_bounding_surfs[AR_bounding_surfs.size() - 1].push_back( surf ); }
                if (i_y == AR_planes[1].size() - 2) {
                    surf.Set(1, AR_planes[1][i_y + 1], {0, 2}, {{AR_planes[0][i_x], AR_planes[0][i_x + 1]}, {AR_planes[2][i_z], AR_planes[2][i_z + 1]}});
                    AR_bounding_surfs[AR_bounding_surfs.size() - 1].push_back( surf ); }
                if (i_z == AR_planes[2].size() - 2) {
                    surf.Set(2, AR_planes[2][i_z + 1], {0, 1}, {{AR_planes[0][i_x], AR_planes[0][i_x + 1]}, {AR_planes[1][i_y], AR_planes[1][i_y + 1]}});
                    AR_bounding_surfs[AR_bounding_surfs.size() - 1].push_back( surf ); }
            }
        }
    }

    // Initiate the ARGridManager
    ARGridManager AR_corner_pts(AR_planes, tol);




    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Getting STL points from file..." << endl;
    // Get the point cooridnates from the STL file. Assign preliminary line and surface/triangle information
    int tri_count;
    std::vector<std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>> coord_ln_surf_pairs;
    std::vector<std::vector<double>> surf_normals;
    getTriPointsFromSTL(STL_file_path, geo_scale, coord_ln_surf_pairs, surf_normals, tri_count);
    cout << "gmsh_utils.cpp: " << FUNCTION_NAME << ": Complete!" << endl;
    


    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Creating gmsh points from STL points..." << endl;
    // Sort coord_ln_pairs based on the point coordinates (so that duplicates are next to each other)
    std::sort(coord_ln_surf_pairs.begin(), coord_ln_surf_pairs.end(),
        [](const std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> &a,
           const std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> &b)
        {
            if (std::get<0>(a)[0] == std::get<0>(b)[0]) {
                if (std::get<0>(a)[1] == std::get<0>(b)[1]) {
                    return std::get<0>(a)[2] < std::get<0>(b)[2];
                }
                else { return std::get<0>(a)[1] < std::get<0>(b)[1]; }
            }
            else { return std::get<0>(a)[0] < std::get<0>(b)[0]; }
        });
    
    // Initialize ln_pt_surf_pairs
    std::vector<std::tuple<std::vector<int>, std::vector<int>>> ln_pt_surf_pairs; // given a line number, this contains (1) the corresponding pt tags, and (2) the fictitious surface tags
    for (int i = 0; i < coord_ln_surf_pairs.size(); i++) { ln_pt_surf_pairs.push_back( std::tuple<std::vector<int>, std::vector<int>>() ); }

    // Initialize surf_ln_pt_pairs
    std::vector<std::tuple<std::vector<int>, std::vector<int>>> surf_ln_pt_pairs; // given a fictitious surface tag, this contains (1) the corresponding pt tags, and (2) the fictitious surface tags
    for (int i = 0; i < tri_count; i++) { surf_ln_pt_pairs.push_back( std::tuple<std::vector<int>, std::vector<int>>() ); }

    // Create a vector where given the gmsh point tag, it gives the AR volumes that the point belongs to. This should be updated everytime a gmsh point is defined
    DisplacedVector<std::vector<int>> gpt_vol_pairs;

    // Create gmsh points from coord_ln_pairs and collect the tags (associated with each line) in ln_pt_surf_pairs. Also store the surfaces of the corresponding lines in ln_pt_surf_pairs
    int pt_gID, pt_gID_first;
    for (int i_coord = 0; i_coord < coord_ln_surf_pairs.size(); i_coord++) {
        cout << "\r    Considering point  i_coord = " << i_coord << "/" << coord_ln_surf_pairs.size() - 1 << "." << std::flush;
        std::vector<double> &coords_current = std::get<0>(coord_ln_surf_pairs[i_coord]);
        if (i_coord == 0) {
            pt_gID = gmsh::model::occ::addPoint(coords_current[0], coords_current[1], coords_current[2], 1.0); // Create the point
            gpt_vol_pairs.SetMinInd(pt_gID);
            gpt_vol_pairs.push_back( assignARIndex(coords_current, AR_planes, tol) );
            AR_corner_pts.checkAndAddPoint({coords_current, pt_gID});
            //isARCornerPoint({coords_current, pt_gID}, AR_corner_pts, tol);
            pt_gID_first = pt_gID;
        }
        else {
            std::vector<double> &coords_prev = std::get<0>(coord_ln_surf_pairs[i_coord - 1]);
            if (coords_current[0] == coords_prev[0] && coords_current[1] == coords_prev[1] && coords_current[2] == coords_prev[2]) { } // pass; do not create a new point if the coordinates match the previous point
            else {
                pt_gID = gmsh::model::occ::addPoint(coords_current[0], coords_current[1], coords_current[2], 1.0); // Create the point
                gpt_vol_pairs.push_back( assignARIndex(coords_current, AR_planes, tol) );
                AR_corner_pts.checkAndAddPoint({coords_current, pt_gID});
                //isARCornerPoint({coords_current, pt_gID}, AR_corner_pts, tol);
            } 
        }
        
        std::vector<int> &ln_inds = std::get<1>(coord_ln_surf_pairs[i_coord]), &surf_inds = std::get<2>(coord_ln_surf_pairs[i_coord]);
        
        for (int i_ln = 0; i_ln < ln_inds.size(); i_ln++) {
            // Add the point to corresponding lines
            std::get<0>( ln_pt_surf_pairs[ln_inds[i_ln]] ).push_back( pt_gID );
            // Add the point's surfaces to corresponding lines
            std::vector<int> &surf_list = std::get<1>( ln_pt_surf_pairs[ln_inds[i_ln]] );
            surf_list.insert( surf_list.end(), surf_inds.begin(), surf_inds.end() );
            //  ==================================================
            // Consider adding the coordinates to ln_pt_surf_pairs so that they do not have to be continuously retrieved in the next step while handling the lines.
            //  ==================================================
        }

        for (int i_surf = 0; i_surf < surf_inds.size(); i_surf++) {
            // Add the point to corresponding surfaces
            std::get<1>( surf_ln_pt_pairs[surf_inds[i_surf]] ).push_back( pt_gID );
        }
    }

    // Because of 1.) the way in which points in coord_ln_surf_pairs are initially assigned line and surface numbers (in getTriPointsFromSTL()), and
    //            2.) the method of filling ln_pt_surf_pairs with surface information,
    // The resulting lines in ln_pt_surf_pairs will be associated with 2 surface IDs that are identical. We remove this duplication here
    for (int i_ln = 0; i_ln < ln_pt_surf_pairs.size(); i_ln++) {
        std::vector<int> &surf_list = std::get<1>(ln_pt_surf_pairs[i_ln]);
        removeDuplicates(surf_list, [](const int &a, const int &b){ return a < b; });
        assert (surf_list.size() == 1);
    }
    cout << endl << "    Created " << pt_gID - pt_gID_first << " unique gmsh points." << endl;
    cout << "    Complete!" << endl;
    


    // Synchronize so that the points can be used in gmsh::model::getValue() later
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Synchronizing..." << endl;
    gmsh::model::occ::synchronize();
    cout << "    Complete!" << endl;



    // Remove duplicate point tags from surf_ln_pt_pairs
    for (int i_surf = 0; i_surf < surf_ln_pt_pairs.size(); i_surf++) {
        std::vector<int> &pt_list = std::get<1>(surf_ln_pt_pairs[i_surf]);
        removeDuplicates(pt_list, [](const int &a, const int &b){ return a < b; });
    }
    
    

    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Removing duplicate lines..." << endl;
    // Sort each component of ln_pt_surf_pairs. Put the smaller point tag first
    for (int i_ln = 0; i_ln < ln_pt_surf_pairs.size(); i_ln++) {
        std::vector<int> &pt_tags_dummy = std::get<0>(ln_pt_surf_pairs[i_ln]);
        assert (pt_tags_dummy.size() == 2); // Make sure there are only 2 points per line
        if (pt_tags_dummy[0] > pt_tags_dummy[1]) { pt_tags_dummy = {pt_tags_dummy[1], pt_tags_dummy[0]}; }
    }
    
    // Sort ln_pt_surf_pairs based on the first line tag, and then the second line tag if the first ones are the same
    std::sort(ln_pt_surf_pairs.begin(), ln_pt_surf_pairs.end(),
        [](const std::tuple<std::vector<int>, std::vector<int>> &a, const std::tuple<std::vector<int>, std::vector<int>> &b) {
            if (std::get<0>(a)[0] == std::get<0>(b)[0]) { return std::get<0>(a)[1] < std::get<0>(b)[1]; }
            else { return std::get<0>(a)[0] < std::get<0>(b)[0]; }
        });
    
    // Remove duplicate lines, but compile their fictitious surface tags with the kept line (of the duplicate pair). Also provide surf_ln_pt_pairs with the lines for each surface
    std::vector<std::tuple<std::vector<int>, std::vector<int>>> ln_pt_surf_pairs2; // given a line number, this contains (1) the corresponding pt tags, and (2) the fictitious surface tags
    for (int i_ln = 0; i_ln < ln_pt_surf_pairs.size(); i_ln++) {
        // If it is the first line...
        if (i_ln == 0) {
            // ... add it to the new lines
            ln_pt_surf_pairs2.push_back( ln_pt_surf_pairs[i_ln] );
            // and add it to the corresponding surfaces
            std::vector<int> &surf_entry = std::get<1>(ln_pt_surf_pairs[i_ln]);
            for (int i_surf = 0; i_surf < surf_entry.size(); i_surf++) { std::get<0>( surf_ln_pt_pairs[ surf_entry[i_surf] ] ).push_back( ln_pt_surf_pairs2.size() - 1 ); }
            continue; }
        // If the line is a copy of the previous line, just add the fictitious surface tags to the last entry in ln_pt_surf_pairs2
        //if (std::get<0>(ln_pt_surf_pairs[i_ln])[0] == std::get<0>(ln_pt_surf_pairs[i_ln - 1])[0] && std::get<0>(ln_pt_surf_pairs[i_ln])[1] == std::get<0>(ln_pt_surf_pairs[i_ln - 1])[1]) {
        std::vector<int> &pt_tags = std::get<0>(ln_pt_surf_pairs[i_ln]); std::vector<int> &prev_pt_tags = std::get<0>(ln_pt_surf_pairs[i_ln - 1]);
        if (pt_tags[0] == prev_pt_tags[0] && pt_tags[1] == prev_pt_tags[1]) {
            // Add the fictitious surface tags to the last entry
            std::vector<int> &last_surf_entry = std::get<1>(ln_pt_surf_pairs2[ln_pt_surf_pairs2.size() - 1]);
            last_surf_entry.insert(last_surf_entry.end(), std::get<1>(ln_pt_surf_pairs[i_ln]).begin(), std::get<1>(ln_pt_surf_pairs[i_ln]).end());
            // Add the previous line index to surf_ln_pt_pairs
            std::vector<int> &surf_entry = std::get<1>(ln_pt_surf_pairs[i_ln]);
            for (int i_surf = 0; i_surf < surf_entry.size(); i_surf++) { std::get<0>( surf_ln_pt_pairs[ surf_entry[i_surf] ] ).push_back( ln_pt_surf_pairs2.size() - 1 ); }
            continue;
        }
        // If the line is not the first entry and not a copy, add it to the new lines and the corresponding surfaces
        ln_pt_surf_pairs2.push_back( ln_pt_surf_pairs[i_ln] );
        std::vector<int> &surf_entry = std::get<1>(ln_pt_surf_pairs[i_ln]);
        for (int i_surf = 0; i_surf < surf_entry.size(); i_surf++) { std::get<0>( surf_ln_pt_pairs[ surf_entry[i_surf] ] ).push_back( ln_pt_surf_pairs2.size() - 1 ); }
    }
    ln_pt_surf_pairs.clear();
    cout << "    Complete!" << endl;
    






    

    // Split STL triangles/surfaces that lie across AR cut planes
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Cutting STL triangles with AR interface planes...";
    assert (AR_planes.size() <= 3);
    for (int i_d = 0; i_d < AR_planes.size(); i_d++) {
        // Print the current cut plane
        cout << endl << "    Cutting with ";
        if (i_d == 0) { cout << "x = "; } else if (i_d == 1) { cout << "y = "; } else { cout << "z = "; }

        for (int i_ARp = 0; i_ARp < AR_planes[i_d].size(); i_ARp++) {
            // Print the current cut plane
            cout << AR_planes[i_d][i_ARp] << ", " << std::flush;
            
            // Go through the lines
            std::vector<std::vector<int>> surf_cut_pts; for (int i = 0; i < surf_ln_pt_pairs.size(); i++) { surf_cut_pts.push_back( std::vector<int>() ); }
            int N_lns = ln_pt_surf_pairs2.size();
            std::vector<bool> above_cut(N_lns);
            for (int i_ln = 0; i_ln < N_lns; i_ln++) {
                // Check that each line is only associated with 2 surfaces
                assert (std::get<1>(ln_pt_surf_pairs2[i_ln]).size() == 2);

                // Get the tags of the endpoints of the line. Make sure there are only 2 points for the line
                std::vector<int> &pt_tags = std::get<0>(ln_pt_surf_pairs2[i_ln]);
                assert (pt_tags.size() == 2);
                
                // Get the coordinates of the line endpoints
                std::vector<std::vector<double>> pt_coords; for (int i = 0; i < pt_tags.size(); i++) { pt_coords.push_back( std::vector<double>() ); }
                {
                    std::vector<double> pt_para_coords;
                    for (int i = 0; i < pt_tags.size(); i++) { gmsh::model::getValue(0, pt_tags[i], pt_para_coords, pt_coords[i]); } // I do not know what the parametric coords are. They might not be relevant to points; only lines, surfaces, or volumes?
                }

                // If the endpoints lie on---or on the same side of---the line, the AR cut plane does not intersect it; go to the next line
                if ( (pt_coords[0][i_d] >= AR_planes[i_d][i_ARp] - tol && pt_coords[1][i_d] >= AR_planes[i_d][i_ARp] - tol) ) { above_cut[i_ln] = true; continue; }
                if ( (pt_coords[0][i_d] <= AR_planes[i_d][i_ARp] + tol && pt_coords[1][i_d] <= AR_planes[i_d][i_ARp] + tol) ) { above_cut[i_ln] = false; continue; }
                
                // If the line is split by the AR cut plane, get the coordinates of the intersection point between the line and the AR cut plane, and create the point
                std::vector<double> new_pt_coord = {0.0, 0.0, 0.0};
                for (int i_d2 = 0; i_d2 < AR_planes.size(); i_d2++) {
                    if (i_d2 == i_d) { new_pt_coord[i_d2] = AR_planes[i_d][i_ARp]; }
                    else { new_pt_coord[i_d2] = (pt_coords[1][i_d2] - pt_coords[0][i_d2]) / (pt_coords[1][i_d] - pt_coords[0][i_d]) * (AR_planes[i_d][i_ARp] - pt_coords[0][i_d]) + pt_coords[0][i_d2]; }
                }
                pt_gID = gmsh::model::occ::addPoint(new_pt_coord[0], new_pt_coord[1], new_pt_coord[2], 1.0); // Create the point
                gpt_vol_pairs.push_back( assignARIndex(new_pt_coord, AR_planes, tol) );
                AR_corner_pts.checkAndAddPoint({new_pt_coord, pt_gID});
                //isARCornerPoint({new_pt_coord, pt_gID}, AR_corner_pts, tol);
                
                // Create the new lines (which come from splitting the original line by the AR cut plane)
                std::vector<std::tuple<std::vector<int>, std::vector<int>>> new_lines;
                for (int i_addln = 0; i_addln < 2; i_addln++) {
                    new_lines.push_back( std::tuple<std::vector<int>, std::vector<int>>() );
                    std::get<0>(new_lines[i_addln]) = {pt_tags[i_addln], pt_gID};
                    std::get<1>(new_lines[i_addln]) = std::get<1>(ln_pt_surf_pairs2[i_ln]);
                }
                
                // Add one to the current position in the line stack, and one to the end of the line stack
                ln_pt_surf_pairs2[i_ln] = new_lines[0];
                ln_pt_surf_pairs2.push_back(new_lines[1]);
                
                // Determine whether each line was above or below the AR cut plane
                if (pt_coords[0][i_d] >= AR_planes[i_d][i_ARp] - tol) { above_cut[i_ln] = true; above_cut.push_back( false ); }
                else { above_cut[i_ln] = false; above_cut.push_back( true ); }
                
                // Add the new line and point to the surface stack
                std::vector<int> &surf_list_dummy = std::get<1>(ln_pt_surf_pairs2[i_ln]);
                assert (surf_list_dummy.size() == 2);
                for (int i_surf = 0; i_surf < surf_list_dummy.size(); i_surf++) {
                    std::get<0>(surf_ln_pt_pairs[surf_list_dummy[i_surf]]).push_back( ln_pt_surf_pairs2.size() - 1 );
                    std::get<1>(surf_ln_pt_pairs[surf_list_dummy[i_surf]]).push_back( pt_gID );
                    surf_cut_pts[surf_list_dummy[i_surf]].push_back( pt_gID );
                }
            }
            

            // Synchronize so that the points can be used in gmsh::model::getValue() later
            gmsh::model::occ::synchronize();

            
            for (int i_surf = 0; i_surf < surf_cut_pts.size(); i_surf++) {
                // If surface was not cut, continue to the next surface (if the surface is along the plane (i.e., all points are on the plane), it is not considered cut)
                if (surf_cut_pts[i_surf].size() == 0) { continue; }
                    
                // If there is only 1 cut point, look through the other points and find the other point that is on the cut plane
                if (surf_cut_pts[i_surf].size() == 1) {
                    std::vector<int> &pt_tags = std::get<1>(surf_ln_pt_pairs[i_surf]);                    
                    std::vector<int> pt_tags_additional;
                    for (int i_pt = 0; i_pt < pt_tags.size(); i_pt++) {
                        if (pt_tags[i_pt] == surf_cut_pts[i_surf][0]) { continue; } // Do not consider the point already in surf_cut_pts
                        std::vector<double> pt_coords, pt_para_coords;
                        gmsh::model::getValue(0, pt_tags[i_pt], pt_para_coords, pt_coords); // I do not know what the parametric coords are. They might not be relevant to points; only lines, surfaces, or volumes?
                        if (compareValues(pt_coords[i_d], AR_planes[i_d][i_ARp], tol)) { pt_tags_additional.push_back( pt_tags[i_pt] ); } }
                    assert (pt_tags_additional.size() == 1);
                    surf_cut_pts[i_surf].push_back( pt_tags_additional[0] );
                }
                
                // Each cut should now be defined with exactly 2 points, regardless if the cut plane goes through a previously existing point
                assert (surf_cut_pts[i_surf].size() == 2);
                
                // Create a new surface entry in surf_ln_pt_pairs and surf_normals
                surf_ln_pt_pairs.push_back( std::tuple<std::vector<int>, std::vector<int>>() );
                surf_normals.push_back( std::vector<double>() );
                int new_surf_ind = surf_ln_pt_pairs.size() - 1;

                // Update the new entry in surf_normals
                surf_normals[new_surf_ind] = surf_normals[i_surf];

                // Get the line tags for the surface that was cut
                std::vector<int> &ln_inds = std::get<0>(surf_ln_pt_pairs[i_surf]);
                
                // Divide the lines according to if they were above or below the cut plane. If they were above, change their associated surface index in ln_pt_surf_pairs2 to the new surface index
                std::vector<int> lns_above, lns_below;
                for (int i_ln = 0; i_ln < ln_inds.size(); i_ln++) {
                    assert (ln_inds[i_ln] < above_cut.size());
                    assert (above_cut[ln_inds[i_ln]] == true || above_cut[ln_inds[i_ln]] == false);
                    if (above_cut[ln_inds[i_ln]]) {
                        lns_above.push_back( ln_inds[i_ln] );
                        std::vector<int> &surf_inds_dummy = std::get<1>(ln_pt_surf_pairs2[ln_inds[i_ln]]);
                        for (int i = 0; i < surf_inds_dummy.size(); i++) { if (surf_inds_dummy[i] == i_surf) { surf_inds_dummy[i] = new_surf_ind; } }
                    }
                    else { lns_below.push_back( ln_inds[i_ln] ); }
                }

                // Create a fictitious line from the cut points and add it to lns_above, lns_below, and ln_pt_surf_pairs2
                std::tuple<std::vector<int>, std::vector<int>> new_line;
                std::get<0>(new_line) = surf_cut_pts[i_surf];
                std::get<1>(new_line) = {i_surf, new_surf_ind};
                ln_pt_surf_pairs2.push_back( new_line );
                lns_above.push_back( ln_pt_surf_pairs2.size() - 1 );
                lns_below.push_back( ln_pt_surf_pairs2.size() - 1 );

                // Update surf_ln_pt_pairs with the new information
                // For surface above
                std::get<0>(surf_ln_pt_pairs[new_surf_ind]) = lns_above;
                std::vector<int> pt_tags_above;
                for (int i_ln = 0; i_ln < lns_above.size(); i_ln++) {
                    std::vector<int> &pt_tags_dummy = std::get<0>(ln_pt_surf_pairs2[lns_above[i_ln]]);
                    pt_tags_above.insert(pt_tags_above.end(), pt_tags_dummy.begin(), pt_tags_dummy.end()); }
                // Delete duplicate points
                removeDuplicates(pt_tags_above, [](const int &a, const int &b){ return a < b;});
                //std::sort(pt_tags_above.begin(), pt_tags_above.end()); auto last_above = std::unique(pt_tags_above.begin(), pt_tags_above.end()); pt_tags_above.erase(last_above, pt_tags_above.end());
                std::get<1>(surf_ln_pt_pairs[new_surf_ind]) = pt_tags_above;
                
                // For surface below
                std::get<0>(surf_ln_pt_pairs[i_surf]) = lns_below;
                std::vector<int> pt_tags_below;
                for (int i_ln = 0; i_ln < lns_below.size(); i_ln++) {
                    std::vector<int> &pt_tags_dummy = std::get<0>(ln_pt_surf_pairs2[lns_below[i_ln]]);
                    pt_tags_below.insert(pt_tags_below.end(), pt_tags_dummy.begin(), pt_tags_dummy.end()); }
                // Delete duplicate points
                removeDuplicates(pt_tags_below, [](const int &a, const int &b){ return a < b;});
                //std::sort(pt_tags_below.begin(), pt_tags_below.end()); auto last_below = std::unique(pt_tags_below.begin(), pt_tags_below.end()); pt_tags_below.erase(last_below, pt_tags_below.end());
                std::get<1>(surf_ln_pt_pairs[i_surf]) = pt_tags_below;
            }
        }
    }
    cout << endl << "    Complete!" << endl;
    




    // Initialize surf_gln_pairs
    std::vector<std::vector<std::pair<int, int>>> surf_gln_pairs;
    for (int i = 0; i < surf_ln_pt_pairs.size(); i++) { surf_gln_pairs.push_back( std::vector<std::pair<int, int>>() ); }

    DisplacedVector<std::vector<int>> gln_surf_pairs, gln_gsurf_pairs;

    // Create gmsh lines from ln_pt_surf_pairs2 and collect the tags (associated with each surface) in surf_gln_pairs
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Creating gmsh lines..." << endl;
    int ln_gID;
    for (int i_ln = 0; i_ln < ln_pt_surf_pairs2.size(); i_ln++) {
        cout << "\r    Considering line  i_ln = " << i_ln << "/" << ln_pt_surf_pairs2.size() - 1 << "." << std::flush;

        // Create the line
        ln_gID = gmsh::model::occ::addLine(std::get<0>(ln_pt_surf_pairs2[i_ln])[0], std::get<0>(ln_pt_surf_pairs2[i_ln])[1]);            
        // Initiate gln_surf_pairs if i_ln = 0
        if (i_ln == 0) { gln_surf_pairs.SetMinInd(ln_gID); gln_gsurf_pairs.SetMinInd(ln_gID); }
        // Add the surface indices entry to gln_surf_pairs for the corresponding gline
        gln_surf_pairs.push_back( std::vector<int>() ); gln_gsurf_pairs.push_back( std::vector<int>() );
        gln_surf_pairs[ln_gID] = std::get<1>(ln_pt_surf_pairs2[i_ln]);
        // Add the gline to the corresponding surfaces in surf_gln_pairs
        for (int i_surf = 0; i_surf < std::get<1>(ln_pt_surf_pairs2[i_ln]).size(); i_surf++) {
            surf_gln_pairs[ std::get<1>(ln_pt_surf_pairs2[i_ln])[i_surf] ].push_back( {1, ln_gID} );
        }
    }
    cout << endl << "    Created " << gln_surf_pairs.get_max_ind() - gln_surf_pairs.get_min_ind() << " unique gmsh lines." << endl;
    cout << "    Complete!" << endl;
    

    


    // Synchronize the lines that were created
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Synchronizing..." << endl;
    gmsh::model::occ::synchronize();
    cout << "    Complete!" << endl;




    
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Creating surfaces..." << endl;
    // Create a vector where given the gmsh surface tag, it gives the AR volumes that the surface belongs to. This should be updated everytime a gmsh surface is defined
    DisplacedVector<std::vector<int>> gsurf_vol_pairs;
    DisplacedVector<std::vector<double>> gsurf_normal_pairs;

    std::vector<std::vector<int>> vol_gsurf_pairs; for (int i = 0; i < N_AR; i++) { vol_gsurf_pairs.push_back( std::vector<int>() ); }

    for (int i_surf = 0; i_surf < surf_gln_pairs.size(); i_surf++) {
        cout << "\r    Considering surface  i_surf = " << i_surf << "/" << surf_gln_pairs.size() - 1 << "." << std::flush;

        // Remove duplicate lines and make sure there are at least 3 lines left to create a surface with
        removeDuplicates(surf_gln_pairs[i_surf], [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.second < b.second; });
        assert (surf_gln_pairs[i_surf].size() >= 3);
        
        // Record the first line tag in curveLoop. Then, figure out the orientation of the other lines
        std::vector<int> curveLoop;
        curveLoop.push_back( surf_gln_pairs[i_surf][0].second );
        
        // Get the boundaries of the lines to figure out their orientations
        std::vector<std::vector<std::pair<int, int>>> bndry_pts; for (int i = 0; i < surf_gln_pairs[i_surf].size(); i++) { bndry_pts.push_back( std::vector<std::pair<int, int>>() ); }
        for (int i_ln = 0; i_ln < surf_gln_pairs[i_surf].size(); i_ln++) {
            gmsh::model::getBoundary({surf_gln_pairs[i_surf][i_ln]}, bndry_pts[i_ln], false, false, false);
        }
        

        // Figure out the orientation of the other surface lines
        int target_pt = bndry_pts[0][1].second; // Get the pt tag of the [1] point of the [0] line
        std::vector<int> used_inds = {0};
        std::vector<int> all_pt_tags = {bndry_pts[0][1].second}; // For assigning the surface an AR volume index
        for (int i_ln = 1; i_ln < bndry_pts.size(); i_ln++) {
            for (int i_ln2 = 1; i_ln2 < bndry_pts.size(); i_ln2++) {
                
                // If the line was already added to curveLoop, skip it. Each line can only be used once
                bool wasUsed = false;
                for (int i_u = 0; i_u < used_inds.size(); i_u++) { if (i_ln2 == used_inds[i_u]){ wasUsed = true; break; } }
                if (wasUsed) { continue; }
                
                if (bndry_pts[i_ln2][0].second == target_pt) {
                    curveLoop.push_back( surf_gln_pairs[i_surf][i_ln2].second );
                    target_pt = bndry_pts[i_ln2][1].second;
                    all_pt_tags.push_back( target_pt );
                    used_inds.push_back( i_ln2 );
                    break;
                }
                if (bndry_pts[i_ln2][1].second == target_pt) {
                    curveLoop.push_back( -surf_gln_pairs[i_surf][i_ln2].second );
                    target_pt = bndry_pts[i_ln2][0].second;
                    all_pt_tags.push_back( target_pt );
                    used_inds.push_back( i_ln2 );
                    break;
                }
            }
        }
        // Ensure that all lines were used, and that there are as many points recorded in all_pt_tags as there are lines in bndry_pts
        assert (used_inds.size() == bndry_pts.size()); assert (all_pt_tags.size() == bndry_pts.size());
        
        
        // Create the curve loop and corresponding plane surface
        int curveLoop_tag = gmsh::model::occ::addCurveLoop(curveLoop);
        int surf_gID = gmsh::model::occ::addPlaneSurface({curveLoop_tag});

                    


        // Add the gmsh surface tag to gln_gsurf_pairs
        for (int i_gln = 0; i_gln < surf_gln_pairs[i_surf].size(); i_gln++) { gln_gsurf_pairs[surf_gln_pairs[i_surf][i_gln].second].push_back( surf_gID ); }


        // Figure out which AR volume the surface is in by checking the AR volume of the points/point tags (all_pt_tags) with gpt_vol_pairs[gmsh_pt_tag].
        //cout << "Surface " << i_surf << endl; // For debugging
        if (i_surf == 0) { gsurf_vol_pairs.SetMinInd(surf_gID); gsurf_normal_pairs.SetMinInd(surf_gID); }
        gsurf_normal_pairs.push_back( surf_normals[i_surf] );

        
        std::vector<int> vol_inds;
        for (int i_pt = 0; i_pt < all_pt_tags.size(); i_pt++) {
            // Add the point's volume indices to vol_inds
            vol_inds.insert(vol_inds.end(), gpt_vol_pairs[all_pt_tags[i_pt]].begin(), gpt_vol_pairs[all_pt_tags[i_pt]].end());
            if (vol_inds.size() == 1) { break; }
            //cout << "vol_inds after adding vols: "; for (int i = 0; i < vol_inds.size(); i++) {cout << vol_inds[i] << " ";} cout << endl; // For debugging
            if (i_pt == 0) { continue; }
            // Remove the numbers that only show up once (as well as one of the numbers that have a duplicate in vol_inds)
            std::vector<int> vol_inds2;
            std::sort(vol_inds.begin(), vol_inds.end());
            //cout << "vol_inds after sort: "; for (int i = 0; i < vol_inds.size(); i++) {cout << vol_inds[i] << " ";} cout << endl; // For debugging
            for (int i = 0; i < vol_inds.size(); i++) { if (i == 0) { continue; } if (vol_inds[i] == vol_inds[i - 1]) { vol_inds2.push_back( vol_inds[i] ); } }
            vol_inds = vol_inds2;
            //cout << "vol_inds after deletion: "; for (int i = 0; i < vol_inds.size(); i++) {cout << vol_inds[i] << " ";} cout << endl; // For debugging
        }
        assert (vol_inds.size() > 0);
        // Add the AR volume indices to the corresponding surface in gsurf_vol_pairs
        gsurf_vol_pairs.push_back( vol_inds );
        // Add the gmsh surface tag to the corresponding volumes in vol_gsurf_pairs
        for (int i = 0; i < vol_inds.size(); i++) { vol_gsurf_pairs[vol_inds[i]].push_back( surf_gID ); }
    }
    cout << endl << "    Created " << gsurf_vol_pairs.get_max_ind() - gsurf_vol_pairs.get_min_ind() << " gmsh surfaces." << endl;
    cout << "    Complete!" << endl;





    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Creating AR corner points inside the domain..." << endl;
    AR_corner_pts.updateIsInside();
    std::vector<std::vector<int>> vol_gpt_AR_corners; for (int i = 0; i < N_AR; i++) { vol_gpt_AR_corners.push_back( std::vector<int>() ); }
    for (int i_pt = 0; i_pt < AR_corner_pts.isInside.size(); i_pt++)
    {
        // If the AR corner point is inside the computational domain, and it has not yet been created, create it and add the point to the various variables
        if (AR_corner_pts.isInside[i_pt] == 1 && std::get<1>(AR_corner_pts.AR_corner_pts[i_pt]) == false) {
            std::vector<double> &ARcoords = std::get<0>(AR_corner_pts.AR_corner_pts[i_pt]);
            pt_gID = gmsh::model::occ::addPoint(ARcoords[0], ARcoords[1], ARcoords[2], 1.0); // Create the point
            gpt_vol_pairs.push_back( assignARIndex(ARcoords, AR_planes, tol) );
            for (int i = 0; i < gpt_vol_pairs[pt_gID].size(); i++) {
                vol_gpt_AR_corners[ gpt_vol_pairs[pt_gID][i] ].push_back( pt_gID );
            }
        }
    }
    cout << "    Complete!" << endl;

    



    // Synchronize the curve loops and surfaces that were created
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Synchronizing..." << endl;
    gmsh::model::occ::synchronize();
    cout << "    Complete!" << endl;
    




    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Finding additional points, lines, and surfaces to close the volumes in each AR..." << endl;
    // Go through the AR boxes. Get the unique lines of the surfaces in each volume, and particularly, the unique lines on the 3 surfaces of the planes defined with smaller coordinate values
    // These lines, along with lines created by parts of the AR planes, will create surfaces that close the volumes in each AR
    
    // Initiate a variable to hold the gmsh line tags for unique lines on each bounding surface.
    std::vector<std::vector<std::vector<int>>> bounding_surfs_glns;
    std::vector<std::tuple< std::pair<int, int>, int, int>> bounding_surfs_new_lns_gtags; // in the tuple is get<0> = pt tags, get<1> = i_AR, get<2> = i_ARbs
    for (int i_AR = 0; i_AR < AR_bounding_surfs.size(); i_AR++) {
        bounding_surfs_glns.push_back( std::vector<std::vector<int>>() );
        for (int i_ARbs = 0; i_ARbs < AR_bounding_surfs[i_AR].size(); i_ARbs++) {
            bounding_surfs_glns[i_AR].push_back( std::vector<int>() );
        }
    }

    
    std::vector<std::pair<int, std::vector<double>>> new_lns_coords;
    std::vector<std::pair<int, int>> new_lns_gtags;
    for (int i_AR = 0; i_AR < AR_bounding_surfs.size(); i_AR++) {
        cout << "\r    Considering AR  i_AR = " << i_AR << "/" << AR_bounding_surfs.size() - 1 << std::flush;
        
        // Get the lines that bound the group of contiguous STL surfaces in the AR volume
        // These lines will be the unique lines in the set of all STL triangle lines within the AR volume
        // For STL's of closed geometry, such lines should lie on the AR bounding surfaces 
        std::vector<std::pair<int, int>> bndry_lns;
        {
            std::vector<std::pair<int, int>> surfs_gmsh_format, bndry_lns_dummy;
            std::vector<int> &gsurf_tags = vol_gsurf_pairs[i_AR];
            for (int i = 0; i < gsurf_tags.size(); i++) { surfs_gmsh_format.push_back({2, gsurf_tags[i]}); }            
            gmsh::model::getBoundary(surfs_gmsh_format, bndry_lns_dummy, false, false, false);

            std::sort(bndry_lns_dummy.begin(), bndry_lns_dummy.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.second < b.second; });
            bool skip = false;
            for (int i = 0; i < bndry_lns_dummy.size(); i++) {
                if (skip) { if (i < bndry_lns_dummy.size() - 1) { assert (bndry_lns_dummy[i].second != bndry_lns_dummy[i + 1].second); } skip = false; continue; }
                if (i == bndry_lns_dummy.size() - 1) { bndry_lns.push_back( bndry_lns_dummy[i] ); continue; }
                if (bndry_lns_dummy[i].second == bndry_lns_dummy[i + 1].second) { skip = true; continue; }
                bndry_lns.push_back( bndry_lns_dummy[i] );
            }
        }
        
        // If there are unique lines, get the coordinates of the line end points
        std::vector<std::vector<std::pair<int, std::vector<double>>>> pt_coords;
        for (int i_ln = 0; i_ln < bndry_lns.size(); i_ln++) {
            pt_coords.push_back( std::vector<std::pair<int, std::vector<double>>>() );
            std::vector<std::pair<int, int>> bndry_pts;
            gmsh::model::getBoundary({bndry_lns[i_ln]}, bndry_pts, false, false, false);
            for (int i_pt = 0; i_pt < bndry_pts.size(); i_pt++) {
                pt_coords[i_ln].push_back( std::pair<int, std::vector<double>>() );
                pt_coords[i_ln][i_pt].first = bndry_pts[i_pt].second;
                std::vector<double> pt_para_coords;
                gmsh::model::getValue(0, bndry_pts[i_pt].second, pt_para_coords, pt_coords[i_ln][i_pt].second); // I do not know what the parametric coords are. They might not be relevant to points; only lines, surfaces, or volumes?
            }
        }
        
        // Add the AR corner points (the ones that have not yet been defined) to pt_coords. The AR might contain no unique lines, but can still contain corner points
        for (int i_pt = 0; i_pt < vol_gpt_AR_corners[i_AR].size(); i_pt++) {
            std::vector<std::pair<int, std::vector<double>>> new_pt; new_pt.push_back( std::pair<int, std::vector<double>>() );
            std::vector<double> pt_para_coords;
            
            new_pt[0].first = vol_gpt_AR_corners[i_AR][i_pt];
            gmsh::model::getValue(0, new_pt[0].first, pt_para_coords, new_pt[0].second);
            pt_coords.push_back( new_pt );
        }



        //cout << "Points in the AR (Number added by AR corners = " << vol_gpt_AR_corners[i_AR].size() << "):" << endl;
        //for (int i_ln = 0; i_ln < pt_coords.size(); i_ln++) {
        //    for (int i_pt = 0; i_pt < pt_coords[i_ln].size(); i_pt++) {
        //        cout << pt_coords[i_ln][i_pt].second[0] << " " << pt_coords[i_ln][i_pt].second[1] << " " << pt_coords[i_ln][i_pt].second[2] << endl;
        //    }
        //}



        // If there are no unique line points or AR corner points, then 1.) there are no STL surfaces in the AR or 2.) the STL surfaces in the AR form a closed volume---where the
        // computational space is inside the volume---that does not touch the AR bounding surfaces by more than a line (i.e., the AR bounding surfaces do not cut the STL geometry;
        // the geometry only touchs the AR bounding surfaces by a line or point)
        if (pt_coords.size() == 0) { continue; }

        

        // Now, we go through the AR bounding surfaces and see if the points of each unique line lie on each surface---and in particular, if they lie on the edges of the AR bounding surfaces
        // If points do lie on the AR bounding surface, we will need to make a surface that coincides with the bounding surface and closes the geometry inside the AR
        // In the case where points lie on the edges of the AR bounding surface, we will need to make lines and points to make the surface that closes the geometry inside the AR
        for (int i_ARbs = 0; i_ARbs < AR_bounding_surfs[i_AR].size(); i_ARbs++) {
            cout << "\r    Considering AR  i_AR = " << i_AR << "/" << AR_bounding_surfs.size() - 1 << "  and bounding surface  i_ARbs = " << i_ARbs << "/" << AR_bounding_surfs[i_AR].size() - 1 << "." << std::flush;

            int &dim = AR_bounding_surfs[i_AR][i_ARbs].dim;
            double &val = AR_bounding_surfs[i_AR][i_ARbs].coord;
            std::vector<int> &alt_dim = AR_bounding_surfs[i_AR][i_ARbs].alt_dim;
            std::vector<std::vector<double>> &alt_val = AR_bounding_surfs[i_AR][i_ARbs].alt_dim_min_max;


            // See if any of the unique line points lie on surface i_ARbs in AR_bounding_surfs 
            std::vector<std::pair<int, int>> ln_pt_of_pt_coords_on_bs;
            for (int i_ln = 0; i_ln < pt_coords.size(); i_ln++) {
                for (int i_pt = 0; i_pt < pt_coords[i_ln].size(); i_pt++) {
                    if (compareValues(pt_coords[i_ln][i_pt].second[dim], val, tol)) { ln_pt_of_pt_coords_on_bs.push_back( {i_ln, i_pt} ); }
                }
            }


            // If no unique line points lie on the current bounding surface (AR_bounding_surfs[i_AR][i_ARbs]), go to analyze the next bounding surface
            if (ln_pt_of_pt_coords_on_bs.size() == 0) { continue; }



            for (int i_ln = 0; i_ln < bndry_lns.size(); i_ln++) {
                bounding_surfs_glns[i_AR][i_ARbs].push_back( bndry_lns[i_ln].second );
            }



            
            // If there are unique line points on the current bounding surface, divide those points by which bounding surface edges they lie on
            // If the points are not on the edges, forget about them for now. We will return to them when we create a surface with them/their lines
            //cout << "ln_pt_of_pt_coords_on_bs.size(): " << ln_pt_of_pt_coords_on_bs.size() << endl;


            // Separate the ln-pt pairs by which edge of the bounding surface they lie on
            std::vector<std::vector<std::pair<int, int>>> ln_pt_of_pt_coords_on_bl; for (int i = 0; i < 4; i++) { ln_pt_of_pt_coords_on_bl.push_back( std::vector<std::pair<int, int>>() ); }
            std::vector<std::pair<int, int>> ln_pt_of_pt_coords_off_bl;
            for (int i_pair = 0; i_pair < ln_pt_of_pt_coords_on_bs.size(); i_pair++) {
                
                int &i_ln = ln_pt_of_pt_coords_on_bs[i_pair].first;
                int &i_pt = ln_pt_of_pt_coords_on_bs[i_pair].second;
                
                std::vector<double> &pnt = pt_coords[i_ln][i_pt].second;
                bool notChosen = true;
                if      ( (compareValues(pnt[alt_dim[0]], alt_val[0][0], tol) || pnt[alt_dim[0]] > alt_val[0][0]) &&
                          (compareValues(pnt[alt_dim[0]], alt_val[0][1], tol) || pnt[alt_dim[0]] < alt_val[0][1]) &&
                           compareValues(pnt[alt_dim[1]], alt_val[1][0], tol) ) { ln_pt_of_pt_coords_on_bl[0].push_back( {i_ln, i_pt} ); notChosen = false; }
                else if ( (compareValues(pnt[alt_dim[0]], alt_val[0][0], tol) || pnt[alt_dim[0]] > alt_val[0][0]) &&
                          (compareValues(pnt[alt_dim[0]], alt_val[0][1], tol) || pnt[alt_dim[0]] < alt_val[0][1]) &&
                           compareValues(pnt[alt_dim[1]], alt_val[1][1], tol) ) { ln_pt_of_pt_coords_on_bl[1].push_back( {i_ln, i_pt} ); notChosen = false; }
                if      ( (compareValues(pnt[alt_dim[1]], alt_val[1][0], tol) || pnt[alt_dim[1]] > alt_val[1][0]) &&
                          (compareValues(pnt[alt_dim[1]], alt_val[1][1], tol) || pnt[alt_dim[1]] < alt_val[1][1]) &&
                           compareValues(pnt[alt_dim[0]], alt_val[0][0], tol) ) { ln_pt_of_pt_coords_on_bl[2].push_back( {i_ln, i_pt} ); notChosen = false; }
                else if ( (compareValues(pnt[alt_dim[1]], alt_val[1][0], tol) || pnt[alt_dim[1]] > alt_val[1][0]) &&
                          (compareValues(pnt[alt_dim[1]], alt_val[1][1], tol) || pnt[alt_dim[1]] < alt_val[1][1]) &&
                           compareValues(pnt[alt_dim[0]], alt_val[0][1], tol) ) { ln_pt_of_pt_coords_on_bl[3].push_back( {i_ln, i_pt} ); notChosen = false; }
                if (notChosen) { ln_pt_of_pt_coords_off_bl.push_back( {i_ln, i_pt} ); };
            }

            // Remove both end points of lines that lie completely on a bounding surface edge. Keep track of these points and then remove the points (regardless of the line they are on) entirely from consideration
            for (int i_onln = 0; i_onln < ln_pt_of_pt_coords_on_bl.size(); i_onln++) {
                //cout << "ln_pt_of_pt_coords_on_bl[" << i_onln << "].size() before removing: " << ln_pt_of_pt_coords_on_bl[i_onln].size() << endl;
                //cout << "Points are:" << endl;
                //for (int i = 0; i < ln_pt_of_pt_coords_on_bl[i_onln].size(); i++) {
                //    std::vector<double> &dummy = pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i].second ].second;
                //    cout << dummy[0] << " " << dummy[1] << " " << dummy[2] << endl;
                //}
                
                
                int i_para, i_perp, para_dim, perp_dim, perp_dim_val;
                if      (i_onln == 0) { i_para = 0; i_perp = 1; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][0]; }
                else if (i_onln == 1) { i_para = 0; i_perp = 1; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][1]; }
                else if (i_onln == 2) { i_para = 1; i_perp = 0; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][0]; }
                else if (i_onln == 3) { i_para = 1; i_perp = 0; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][1]; }
                

                // Removing the end points of lines that lie completely on a bounding surface edge, while taking note of the point coords that were deleted (these will be used to delete more points after)
                std::vector<std::vector<double>> remove_coords;
                {
                    std::vector<std::pair<int, int>> new_vec;
                    std::sort(ln_pt_of_pt_coords_on_bl[i_onln].begin(), ln_pt_of_pt_coords_on_bl[i_onln].end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.first < b.first; });
                    bool skip = false;
                    for (int i = 0; i < ln_pt_of_pt_coords_on_bl[i_onln].size(); i++) {
                        if (skip) {
                            if (i != ln_pt_of_pt_coords_on_bl[i_onln].size() - 1) { assert (ln_pt_of_pt_coords_on_bl[i_onln][i].first != ln_pt_of_pt_coords_on_bl[i_onln][i + 1].first); }
                            skip = false;
                            continue;
                        }
                        if (i == ln_pt_of_pt_coords_on_bl[i_onln].size() - 1) { new_vec.push_back( ln_pt_of_pt_coords_on_bl[i_onln][i] ); continue; }
                        if (ln_pt_of_pt_coords_on_bl[i_onln][i].first == ln_pt_of_pt_coords_on_bl[i_onln][i + 1].first) {
                            skip = true;
                            remove_coords.push_back( pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i].second ].second );
                            remove_coords.push_back( pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i + 1].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i + 1].second ].second );
                            //for (int j = 0; j < end_coords.size(); j++) {
                            //    if (compareVectors(end_coords[j], pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i].second ].second), tol) { remove_end_coords[j] = true; }
                            //    if (compareVectors(end_coords[j], pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i + 1].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i + 1].second ].second), tol) { remove_end_coords[j] = true; } }
                            continue;
                        }
                        new_vec.push_back( ln_pt_of_pt_coords_on_bl[i_onln][i] );
                    }
                    ln_pt_of_pt_coords_on_bl[i_onln] = new_vec;
                }

                // If a point was removed because it was part of a line that was completely on the bounding surface edge, remove that point entirely
                {
                    std::vector<std::pair<int, int>> new_vec;
                    for (int i = 0; i < ln_pt_of_pt_coords_on_bl[i_onln].size(); i++) {
                        bool keepIt = true;
                        std::vector<double> &coord_vec = pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i].second ].second;
                        for (int i_r = 0; i_r < remove_coords.size(); i_r++) {
                            if (compareVectors(remove_coords[i_r], coord_vec, tol)) {
                                keepIt = false;
                                break;
                            }
                        }
                        if (keepIt) { new_vec.push_back( ln_pt_of_pt_coords_on_bl[i_onln][i] ); }
                    }
                    ln_pt_of_pt_coords_on_bl[i_onln] = new_vec;
                    
                }
                


                // Maybe get the normal here? and then if the point is an end point that appears twice with 0 normal, we disregard it because it means there is a non-unique
                // line on the bounding edge that contains the point (i.e., STL geometry forms a bracket---whose surfaces have perpendicular normals to the edge direction---on the edge).
                // If the point is an end point that has 0 normal, but only appears once, it is one of the AR corner points that was added, and it should be used to build the AR

                
                // Get the STL surface normal vectors for each point
                std::vector<std::vector<double>> norm_vecs; for (int i = 0; i < ln_pt_of_pt_coords_on_bl[i_onln].size(); i++) { norm_vecs.push_back( std::vector<double>() ); }
                for (int i_pair = 0; i_pair < ln_pt_of_pt_coords_on_bl[i_onln].size(); i_pair++)
                {
                    std::vector<std::vector<double>> end_coords;
                    for (int i = 0; i < 2; i++) {
                        end_coords.push_back( std::vector<double>(3) );
                        end_coords[i][dim] = val; end_coords[i][para_dim] = alt_val[i_para][i]; end_coords[i][perp_dim] = perp_dim_val;
                    }
                
                    
                    // If the pair is one of the added AR corner points, make the normal vector
                    if (ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first >= bndry_lns.size()) {
                        std::vector<double> &pnt_coord = pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].second ].second;
                        std::vector<double> new_norm = {0.0, 0.0, 0.0};
                        if (compareVectors(end_coords[0], pnt_coord, tol)) {
                            new_norm[para_dim] = -1.0;
                            norm_vecs[i_pair] = new_norm; }
                        else if (compareVectors(end_coords[1], pnt_coord, tol)) {
                            new_norm[para_dim] = 1.0;
                            norm_vecs[i_pair] = new_norm; }
                        else {
                            cerr << "gmsh_utils.cpp: " << FUNCTION_NAME << ": CRITICAL ERROR: A point seems to have been added to pt_coords through the 'AR corner point' method, but the point is not actually an AR corner point on the current AR bounding surface edge. ";
                            cerr << "The point is {" << pnt_coord[0] << ", " << pnt_coord[1] << ", " << pnt_coord[2] << "}. i_AR = " << i_AR << ". i_onln = " << i_onln << "." << endl;
                            exit(1); }
                        continue;
                    }

                    int &gln = bndry_lns[ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first].second;
                    // Get the possible gmsh surfaces for the line
                    std::vector<int> &gsurfs = gln_gsurf_pairs[gln]; assert (gsurfs.size() == 2); // There should be 2 surfaces for every line in the closed STL geometry
                    // Each surface belongs to some AR volumes. Go through the surfaces and see which lies in the current AR volume. Since bndry_lns only considers lines appearing once among the surfaces in the current volume, there should only be 1 surface that belongs to the current volume
                    bool foundNormalVec = false;
                    int N_norms_found = 0;
                    for (int i_gs = 0; i_gs < gsurfs.size(); i_gs++) {
                        // Get the possible volumes for the surface
                        std::vector<int> &vols = gsurf_vol_pairs[gsurfs[i_gs]];
                        // See if the surface's volume is the AR volume being analyzed. If so, save the surface normal
                        for (int i_vol = 0; i_vol < vols.size(); i_vol++) { if (vols[i_vol] == i_AR) { norm_vecs[i_pair] = gsurf_normal_pairs[gsurfs[i_gs]]; N_norms_found += 1; } }
                    }
                    assert (N_norms_found == 1);
                }
                
                std::vector<std::tuple< std::vector<double>, std::pair<int, int>, std::vector<double> >> pt_normal_pairs;
                for (int i_pair = 0; i_pair < ln_pt_of_pt_coords_on_bl[i_onln].size(); i_pair++) {
                    pt_normal_pairs.push_back( std::tuple<std::vector<double>, std::pair<int, int>, std::vector<double>>() );
                    std::get<0>(pt_normal_pairs[pt_normal_pairs.size() - 1]) = pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].second ].second;
                    std::get<1>(pt_normal_pairs[pt_normal_pairs.size() - 1]) = ln_pt_of_pt_coords_on_bl[i_onln][i_pair];
                    std::get<2>(pt_normal_pairs[pt_normal_pairs.size() - 1]) = norm_vecs[i_pair];
                }
                std::sort(pt_normal_pairs.begin(), pt_normal_pairs.end(),
                    [](const std::tuple<std::vector<double>, std::pair<int, int>, std::vector<double>> &a,
                       const std::tuple<std::vector<double>, std::pair<int, int>, std::vector<double>> &b) {
                        if (std::get<0>(a)[0] == std::get<0>(b)[0]) {
                            if (std::get<0>(a)[1] == std::get<0>(b)[1]) {
                                return std::get<0>(a)[2] < std::get<0>(b)[2];
                            }
                            else { return std::get<0>(a)[1] < std::get<0>(b)[1]; }
                        }
                        else { return std::get<0>(a)[0] < std::get<0>(b)[0]; }
                    });
                std::vector<std::pair<int, int>> new_vec;
                bool skip = false;
                for (int i_pair = 0; i_pair < pt_normal_pairs.size(); i_pair++) {
                    std::vector<double> &current_coord = std::get<0>(pt_normal_pairs[i_pair]); std::pair<int, int> &current_pair = std::get<1>(pt_normal_pairs[i_pair]);
                    //cout << "Coords and their normals: {" << current_coord[0] << " " << current_coord[1] << " " << current_coord[2] << "}, " << std::get<2>(pt_normal_pairs[i_pair])[para_dim] << endl;
                    if (skip) {
                        if (i_pair < pt_normal_pairs.size() - 1) {
                            std::vector<double> &next_coord = std::get<0>(pt_normal_pairs[i_pair + 1]);
                            assert (!compareVectors(current_coord, next_coord, tol));
                        }
                        skip = false;
                        continue;
                    }
                    if (i_pair == pt_normal_pairs.size() - 1) { new_vec.push_back( current_pair ); continue; }
                    std::vector<double> &next_coord = std::get<0>(pt_normal_pairs[i_pair + 1]);
                    double &current_norm_para = std::get<2>(pt_normal_pairs[i_pair])[para_dim]; double &next_norm_para = std::get<2>(pt_normal_pairs[i_pair + 1])[para_dim];
                    // If two coordinates are the same and have a normal of 0 in the parallel dimension, skip; do not add either
                    if (compareVectors(current_coord, next_coord, tol) && compareValues(current_norm_para, 0.0, tol) && compareValues(next_norm_para, 0.0, tol)) { skip = true; continue; }
                    // If two coordinates are the same and have the same normal in the parallel dimension, add one and skip the other. No need for doubles
                    //if (compareVectors(current_coord, next_coord, tol) && compareValues(current_norm_para, next_norm_para, tol)) { new_vec.push_back( current_pair ); skip = true; continue; }
                    // If two coordinates are the same and the normals in the parallel dimension are pointed the same way, add one and skip the other. No need for doubles
                    if (compareVectors(current_coord, next_coord, tol) && current_norm_para * next_norm_para > 0.0) { new_vec.push_back( current_pair ); skip = true; continue; }
                    new_vec.push_back( current_pair );
                }
                ln_pt_of_pt_coords_on_bl[i_onln] = new_vec;




                //cout << "ln_pt_of_pt_coords_on_bl[" << i_onln << "].size(): " << ln_pt_of_pt_coords_on_bl[i_onln].size() << endl;

                
            }
            




            // For the points on the bounding surface edges, 
            for (int i_onln = 0; i_onln < ln_pt_of_pt_coords_on_bl.size(); i_onln++) {
                // See the normal of the line's STL. If the normal is not 0 in the dim direction, then march to the next point in that direction (either another point or the AR corner). Make the line, and the point in the case it is the corner.
                
                int i_para, i_perp, para_dim, perp_dim, perp_dim_val;
                if      (i_onln == 0) { i_para = 0; i_perp = 1; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][0]; }
                else if (i_onln == 1) { i_para = 0; i_perp = 1; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][1]; }
                else if (i_onln == 2) { i_para = 1; i_perp = 0; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][0]; }
                else if (i_onln == 3) { i_para = 1; i_perp = 0; para_dim = alt_dim[i_para]; perp_dim = alt_dim[i_perp]; perp_dim_val = alt_val[i_perp][1]; }
                    
                std::vector<std::vector<double>> end_coords;
                for (int i = 0; i < 2; i++) {
                    end_coords.push_back( std::vector<double>(3) );
                    end_coords[i][dim] = val; end_coords[i][para_dim] = alt_val[i_para][i]; end_coords[i][perp_dim] = perp_dim_val;
                }
                


                //cout << "i_onln: " << i_onln << endl;
                //cout << "Points:" << endl;
                //for (int i_pair = 0; i_pair < ln_pt_of_pt_coords_on_bl[i_onln].size(); i_pair++) {
                //    std::vector<double> &pnt_coord = pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].second ].second;
                //    cout << pnt_coord[0] << " " << pnt_coord[1] << " " << pnt_coord[2] << endl;
                //}
                


                // Get the STL surface normal vectors for each point
                std::vector<std::vector<double>> norm_vecs; for (int i = 0; i < ln_pt_of_pt_coords_on_bl[i_onln].size(); i++) { norm_vecs.push_back( std::vector<double>() ); }
                for (int i_pair = 0; i_pair < ln_pt_of_pt_coords_on_bl[i_onln].size(); i_pair++) {
                    // If the pair is one of the added AR corner points, make the normal vector 0
                    //if (ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first >= bndry_lns.size()) { norm_vecs[i_pair] = {0.0, 0.0, 0.0}; continue; }
                    if (ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first >= bndry_lns.size()) {
                        std::vector<double> &pnt_coord = pt_coords[ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first ][ ln_pt_of_pt_coords_on_bl[i_onln][i_pair].second ].second;
                        std::vector<double> new_norm = {0.0, 0.0, 0.0};
                        if (compareVectors(end_coords[0], pnt_coord, tol)) {
                            new_norm[para_dim] = -1.0;
                            norm_vecs[i_pair] = new_norm; }
                        else if (compareVectors(end_coords[1], pnt_coord, tol)) {
                            new_norm[para_dim] = 1.0;
                            norm_vecs[i_pair] = new_norm; }
                        else {
                            cerr << "gmsh_utils.cpp: " << FUNCTION_NAME << ": CRITICAL ERROR: A point seems to have been added to pt_coords through the 'AR corner point' method, but the point is not actually an AR corner point on the current AR bounding surface edge. ";
                            cerr << "The point is {" << pnt_coord[0] << ", " << pnt_coord[1] << ", " << pnt_coord[2] << "}. i_AR = " << i_AR << ". i_onln = " << i_onln << "." << endl;
                            exit(1); }
                        continue;
                    }

                    int &gln = bndry_lns[ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first].second;
                    // Get the possible gmsh surfaces for the line
                    std::vector<int> &gsurfs = gln_gsurf_pairs[gln]; assert (gsurfs.size() == 2); // There should be 2 surfaces for every line in the closed STL geometry
                    // Each surface belongs to some AR volumes. Go through the surfaces and see which lies in the current AR volume. Since bndry_lns only considers lines appearing once among the surfaces in the current volume, there should only be 1 surface that belongs to the current volume
                    bool foundNormalVec = false;
                    int N_norms_found = 0;
                    for (int i_gs = 0; i_gs < gsurfs.size(); i_gs++) {
                        // Get the possible volumes for the surface
                        std::vector<int> &vols = gsurf_vol_pairs[gsurfs[i_gs]];
                        // See if the surface's volume is the AR volume being analyzed. If so, save the surface normal
                        for (int i_vol = 0; i_vol < vols.size(); i_vol++) { if (vols[i_vol] == i_AR) { norm_vecs[i_pair] = gsurf_normal_pairs[gsurfs[i_gs]]; N_norms_found += 1; } }
                    }
                    assert (N_norms_found == 1);
                }
                

                // Now consider each pair in ln_pt_of_pt_coords_on_bl and use each pair
                std::vector<int> used_i_pair;
                for (int i_pair = 0; i_pair < ln_pt_of_pt_coords_on_bl[i_onln].size(); i_pair++) {
                    // If the i_pair was already considered as the end point to a new ln, go to the next i_pair
                    bool wasUsed = false; for (int i = 0; i < used_i_pair.size(); i++) { if (i_pair == used_i_pair[i]) { wasUsed = true; break; } } if (wasUsed) { continue; }
                    
                    int &ln_ind = ln_pt_of_pt_coords_on_bl[i_onln][i_pair].first;
                    int &pt_ind = ln_pt_of_pt_coords_on_bl[i_onln][i_pair].second;
                    
                    
                    double &start_pt_coord = pt_coords[ln_ind][pt_ind].second[para_dim];
                    int &start_gpt = pt_coords[ln_ind][pt_ind].first;
                    

                    
                    //cout << "i_onln: " << i_onln << ", i_pair: " << i_pair << " of " << ln_pt_of_pt_coords_on_bl[i_onln].size() - 1 << ", " << pt_coords[ln_ind][pt_ind].second[0] << " " << pt_coords[ln_ind][pt_ind].second[1] << " " << pt_coords[ln_ind][pt_ind].second[2] << ", " << norm_vecs[i_pair][para_dim] << endl;



                    if (compareValues(norm_vecs[i_pair][para_dim], 0.0, tol)) { used_i_pair.push_back( i_pair ); continue; }
                    else if (norm_vecs[i_pair][para_dim] < 0) {  // NOTE: Normals are assumed to point outward of the pore-space; so we want to move in the opposite direction
                        // Make the appropriate corner point of the AR the first potential end point 
                        std::vector<double> end_pt_coords(3); end_pt_coords[dim] = val; end_pt_coords[para_dim] = alt_val[i_para][1]; end_pt_coords[perp_dim] = perp_dim_val;

                        // If the point (start_pt) is the end point for this direction of normal, then consider the point used and skip
                        if (compareValues(start_pt_coord, end_pt_coords[para_dim], tol)) { continue; }

                        int end_pair = -1, end_gpt = -1;
                        // See if there is a closer end point in the pairs
                        for (int i_pair2 = 0; i_pair2 < ln_pt_of_pt_coords_on_bl[i_onln].size(); i_pair2++) {
                            // If the i_pair2 was already considered as an end point, go to the next i_pair2
                            bool wasUsed2 = false; for (int i = 0; i < used_i_pair.size(); i++) { if (i_pair2 == used_i_pair[i] || i_pair == i_pair2) { wasUsed2 = true; break; } } if (wasUsed2) { continue; }

                            int &ln_ind2 = ln_pt_of_pt_coords_on_bl[i_onln][i_pair2].first;
                            int &pt_ind2 = ln_pt_of_pt_coords_on_bl[i_onln][i_pair2].second;
                            // If the potential end point is closer to the start point than the end point, make it the new end point 
                            std::vector<double> &potential_end_pt_coords = pt_coords[ln_ind2][pt_ind2].second;
                            int potential_end_gpt = pt_coords[ln_ind2][pt_ind2].first;
                            if (potential_end_pt_coords[para_dim] > start_pt_coord && (compareValues(potential_end_pt_coords[para_dim],  end_pt_coords[para_dim], tol) || potential_end_pt_coords[para_dim] < end_pt_coords[para_dim])) {
                                end_pt_coords = potential_end_pt_coords; end_pair = i_pair2; end_gpt = potential_end_gpt; }
                        }
                        if (end_pair > -1) {
                            used_i_pair.push_back( end_pair );
                            new_lns_gtags.push_back( {start_gpt, end_gpt} );
                            std::tuple<std::pair<int, int>, int, int> dummy_entry;
                            std::get<0>(dummy_entry) = {start_gpt, end_gpt}; std::get<1>(dummy_entry) = i_AR; std::get<2>(dummy_entry) = i_ARbs;
                            bounding_surfs_new_lns_gtags.push_back( dummy_entry );
                        }
                        else { exit(1); new_lns_coords.push_back( {start_gpt, end_pt_coords} ); }
                        used_i_pair.push_back( i_pair );
                    }
                    else if (norm_vecs[i_pair][para_dim] > 0) {
                        // Make the first potential end point the appropriate corner point of the AR
                        std::vector<double> end_pt_coords(3); end_pt_coords[dim] = val; end_pt_coords[para_dim] = alt_val[i_para][0]; end_pt_coords[perp_dim] = perp_dim_val;
                        
                        // If the point (start_pt) is the end point for this direction of normal, then consider the point used and skip
                        if (compareValues(start_pt_coord, end_pt_coords[para_dim], tol)) { continue; }

                        int end_pair = -1, end_gpt = -1;
                        // See if there is a closer end point in the pairs
                        for (int i_pair2 = 0; i_pair2 < ln_pt_of_pt_coords_on_bl[i_onln].size(); i_pair2++) {
                            // If the i_pair2 was already considered as an end point, go to the next i_pair2
                            bool wasUsed2 = false; for (int i = 0; i < used_i_pair.size(); i++) { if (i_pair2 == used_i_pair[i] || i_pair == i_pair2) { wasUsed2 = true; break; } } if (wasUsed2) { continue; }

                            int &ln_ind2 = ln_pt_of_pt_coords_on_bl[i_onln][i_pair2].first;
                            int &pt_ind2 = ln_pt_of_pt_coords_on_bl[i_onln][i_pair2].second;
                            // If the potential end point is closer to the start point than the end point, make it the new end point 
                            std::vector<double> &potential_end_pt_coords = pt_coords[ln_ind2][pt_ind2].second;
                            int potential_end_gpt = pt_coords[ln_ind2][pt_ind2].first;
                            if (potential_end_pt_coords[para_dim] < start_pt_coord && (compareValues(potential_end_pt_coords[para_dim],  end_pt_coords[para_dim], tol) || potential_end_pt_coords[para_dim] > end_pt_coords[para_dim])) {
                                end_pt_coords = potential_end_pt_coords; end_pair = i_pair2; end_gpt = potential_end_gpt; }
                        }
                        if (end_pair > -1) {
                            used_i_pair.push_back( end_pair );
                            new_lns_gtags.push_back( {start_gpt, end_gpt} );
                            std::tuple<std::pair<int, int>, int, int> dummy_entry;
                            std::get<0>(dummy_entry) = {start_gpt, end_gpt}; std::get<1>(dummy_entry) = i_AR; std::get<2>(dummy_entry) = i_ARbs;
                            bounding_surfs_new_lns_gtags.push_back( dummy_entry );
                        }
                        else { exit(1); new_lns_coords.push_back( {start_gpt, end_pt_coords} ); }
                        used_i_pair.push_back( i_pair );
                    }
                }
                //cout << used_i_pair.size() << " " << ln_pt_of_pt_coords_on_bl[i_onln].size() << endl;
                assert (used_i_pair.size() == ln_pt_of_pt_coords_on_bl[i_onln].size());
            }



        }
    }
    cout << endl << "    Complete!" << endl;





    // No new coordinates should need to be created since we are now creating the AR corner points before the preivous for loop
    assert (new_lns_coords.size() == 0);





    // Synchronize the curve loops and surfaces that were created
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Synchronizing..." << endl;
    gmsh::model::occ::synchronize();
    cout << "    Complete!" << endl;




    
    cout << "gmsh_utils.cpp: " << FUNCTION_NAME << ": Creating additional lines to close the volumes in each AR..." << endl;
    // Sort the gmsh point tags for each line
    for (int i_ln = 0; i_ln < bounding_surfs_new_lns_gtags.size(); i_ln++) {
        std::pair<int, int> &pt_pair = std::get<0>(bounding_surfs_new_lns_gtags[i_ln]);
        if (pt_pair.first > pt_pair.second) {
            int gtag_temp = pt_pair.second;
            pt_pair.second = pt_pair.first;
            pt_pair.first = gtag_temp;
        }
    }
    // Sort the lines in new_lns_gtags by the tags.
    std::sort(bounding_surfs_new_lns_gtags.begin(), bounding_surfs_new_lns_gtags.end(),
        [](const std::tuple<std::pair<int, int>, int, int> &a, const std::tuple<std::pair<int, int>, int, int> &b) {
            if (std::get<0>(a).first == std::get<0>(b).first) { return std::get<0>(a).second < std::get<0>(b).second; }
            else { return std::get<0>(a).first < std::get<0>(b).first; }
        });

    // Create gmsh lines from new_lns_gtags and collect the tags (associated with each surface) in ???
    int N_new_lns_gtags = 0;
    for (int i_ln = 0; i_ln < bounding_surfs_new_lns_gtags.size(); i_ln++) {
        
        std::pair<int, int> &current_pt_pair = std::get<0>(bounding_surfs_new_lns_gtags[i_ln]);
        
        if (i_ln == 0) {
            ln_gID = gmsh::model::occ::addLine(current_pt_pair.first, current_pt_pair.second);
            int &i_AR = std::get<1>(bounding_surfs_new_lns_gtags[i_ln]);
            int &i_ARbs = std::get<2>(bounding_surfs_new_lns_gtags[i_ln]);
            bounding_surfs_glns[i_AR][i_ARbs].push_back( ln_gID );
            N_new_lns_gtags += 1;
            continue;
        }
        
        std::pair<int, int> &prev_pt_pair = std::get<0>(bounding_surfs_new_lns_gtags[i_ln - 1]);
        if (current_pt_pair.first == prev_pt_pair.first && current_pt_pair.second == prev_pt_pair.second) { }
        else { ln_gID = gmsh::model::occ::addLine(current_pt_pair.first, current_pt_pair.second); }
        
        int &i_AR = std::get<1>(bounding_surfs_new_lns_gtags[i_ln]);
        int &i_ARbs = std::get<2>(bounding_surfs_new_lns_gtags[i_ln]);        
        bounding_surfs_glns[i_AR][i_ARbs].push_back( ln_gID );
        N_new_lns_gtags += 1;
    }
    cout << "    Created " << N_new_lns_gtags << " additional lines." << endl;
    cout << "    Complete!" << endl;




    
    // Synchronize the curve loops and surfaces that were created
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Synchronizing..." << endl;
    gmsh::model::occ::synchronize();
    cout << "    Complete!" << endl;





    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Creating additional surfaces to close the volumes in each AR..." << endl;
    // Now that the lines have been created, go through bounding_surfs_glns and create the surfaces to close the geometries in each AR
    int surf_gID;
    for (int i_AR = 0; i_AR < bounding_surfs_glns.size(); i_AR++) {
        for (int i_ARbs = 0; i_ARbs < bounding_surfs_glns[i_AR].size(); i_ARbs++) {
            cout << "\r    Considering AR  i_AR = " << i_AR << "/" << bounding_surfs_glns.size() - 1 << "  and surface bounding  i_ARbs = " << i_ARbs << "/" << bounding_surfs_glns[i_AR].size() - 1 << "." << std::flush;
            
            // If there are no gmsh lines on the AR bounding surface, go to the next bounding surface
            if (bounding_surfs_glns[i_AR][i_ARbs].size() == 0) { continue; }
            
            // Get the dim and val of the bounding surface
            int &dim = AR_bounding_surfs[i_AR][i_ARbs].dim;
            double &val = AR_bounding_surfs[i_AR][i_ARbs].coord;
            std::vector<int> &alt_dim = AR_bounding_surfs[i_AR][i_ARbs].alt_dim;
            std::vector<std::vector<double>> &alt_val = AR_bounding_surfs[i_AR][i_ARbs].alt_dim_min_max;
            //cout << "dim: " << dim << ", val: " << val << endl;

            // Get the gmsh tags of the end points of each line. These will be used to
            // figure out how the lines should be oriented. Neglect the lines with points
            // that do not have a coordinate equal to val at dimension dim
            std::vector<std::vector<std::pair<int, int>>> bndry_pts;
            std::vector<int> bndry_glns, bndry_vol;
            bool first_line = true;
            for (int i_ln = 0; i_ln < bounding_surfs_glns[i_AR][i_ARbs].size(); i_ln++) {
                // Get the gmsh tags of the end points
                std::vector<std::pair<int, int>> bndry_pts_dummy;
                gmsh::model::getBoundary({{1, bounding_surfs_glns[i_AR][i_ARbs][i_ln]}}, bndry_pts_dummy, false, false, false);
                
                // Get the coordinates of the end points
                std::vector<double> pt_para_coords; std::vector<std::vector<double>> end_coords;
                for (int i = 0; i < bndry_pts_dummy.size(); i++) { end_coords.push_back( std::vector<double>() );
                    gmsh::model::getValue(0, bndry_pts_dummy[i].second, pt_para_coords, end_coords[i]); }
                
                // Only keep the end points if the dim coordinate of both points is val
                if (!compareValues(end_coords[0][dim], val, tol) || !compareValues(end_coords[1][dim], val, tol)) { continue; }

                bndry_pts.push_back( bndry_pts_dummy );
                bndry_glns.push_back( bounding_surfs_glns[i_AR][i_ARbs][i_ln] );
                //cout << "i_ln: " << i_ln << endl;
                //cout << "{" << end_coords[0][0] << ", " << end_coords[0][1] << ", " << end_coords[0][2] << "}" << endl;
                //cout << "{" << end_coords[1][0] << ", " << end_coords[1][1] << ", " << end_coords[1][2] << "}" << endl;
                
                // Get the volume indices for the points; only keep them if the volume is shared between all points
                for (int i_pt = 0; i_pt < bndry_pts_dummy.size(); i_pt++) {
                    std::vector<int> &vol_inds = gpt_vol_pairs[bndry_pts_dummy[i_pt].second];
                    bndry_vol.insert(bndry_vol.end(), vol_inds.begin(), vol_inds.end());
                }
                assert (bndry_vol.size() != 0);

                //cout << "bndry_vol: " << endl;
                //for (int i = 0; i < bndry_vol.size(); i++) {
                //    cout << bndry_vol[i] << endl;
                //}
                
                
                std::vector<int> bndry_vol2;
                std::sort(bndry_vol.begin(), bndry_vol.end());
                //cout << "bndry_vol sorted: i_ln " << i_ln << endl;
                //for (int i = 0; i < bndry_vol.size(); i++) {
                //    cout << bndry_vol[i] << endl;
                //}

                if (first_line) {
                    for (int i = 0; i < bndry_vol.size(); i++) {
                        if (i == 0) { continue; }
                        if (bndry_vol[i] == bndry_vol[i - 1]) { bndry_vol2.push_back( bndry_vol[i] ); }
                    }
                    first_line = false;
                } else {
                    for (int i = 0; i < bndry_vol.size(); i++) {
                        if (i == 0 || i == 1) { continue; }
                        if (bndry_vol[i] == bndry_vol[i - 1] && bndry_vol[i] == bndry_vol[i - 2]) { bndry_vol2.push_back( bndry_vol[i] ); }
                    }
                }
                
                bndry_vol = bndry_vol2;
                //cout << "bndry_vol after delete duplicates: " << endl;
                //for (int i = 0; i < bndry_vol.size(); i++) {
                //    cout << bndry_vol[i] << endl;
                //}
                
            }
            assert (bndry_vol.size() > 0);

            
            // Find the curve loops in the set of glns on the AR bounding surface
            std::vector<std::vector<int>> curveLoops = findCurveLoops(bndry_glns, tol, false);
            //cout << "Number of curve loops found: " << curveLoops.size() << endl;
            if (curveLoops.size() == 0) { continue; }



            if (surf_gID == 69766) {
                cout << endl;
                cout << "bndry_glns.size(): " << bndry_glns.size() << endl;
                for (int i = 0; i < bndry_glns.size(); i++) {
                    cout << " " << bndry_glns[i];
                }
                cout << endl;
                

                cout << "curveLoops.size(): " << curveLoops.size() << endl;
                for (int i = 0; i < curveLoops.size(); i++) {
                    cout << "curveLoop " << i << ":";
                    for (int jj = 0; jj < curveLoops[i].size(); jj++) {
                        cout << " " << curveLoops[i][jj];
                    }
                    cout << endl;
                }
                exit(1);
            }

            // Do raycasting if there is more than 1 curveloop to figure out which curve loops will create a surface together
            std::vector<std::vector<std::vector<int>>> curveLoops_arranged;
            if (curveLoops.size() > 1) { curveLoops_arranged = arrangeCurveLoopsForSurfaces(curveLoops, dim, tol, false); }
            else { curveLoops_arranged = {curveLoops}; }



            // Create the gmsh curve loops and corresponding plane surfaces.
            // If all points of the curve loop lie on the boundary of the bounding surface, the bounding surface could already be covered with defined STL surfaces. In
            // this case, look at the coordinates of one of the surfaces of the lines. If all points lie on the bounding surface, then do not define a surface with the curve loop
            
            for (int i_clg = 0; i_clg < curveLoops_arranged.size(); i_clg++) {
                
                std::vector<int> curveLoop_group_tags;
                
                for (int i_cl = 0; i_cl < curveLoops_arranged[i_clg].size(); i_cl++)
                {
                    // Get the unique gmsh point tags of the lines in the curveloop
                    std::vector<int> glns;
                    for (int i_ln = 0; i_ln < curveLoops_arranged[i_clg][i_cl].size(); i_ln++) {
                        int gln = curveLoops_arranged[i_clg][i_cl][i_ln]; if (gln < 0) { gln *= -1.0; } glns.push_back( gln ); }
                    removeDuplicates(glns, [](const int &a, const int &b){ return a < b; });
                    
                    // Get the coordinates of the points in the curve loop
                    std::vector<std::vector<double>> curveLoop_pt_coords;
                    for (int i_ln = 0; i_ln < glns.size(); i_ln++) {
                        // Get the gmsh point tags of the gmsh line end points
                        std::vector<std::pair<int, int>> gpts_dummy; gmsh::model::getBoundary({{1, glns[i_ln]}}, gpts_dummy, false, false, false);
                        // Get the coordinates of the end points
                        std::vector<double> pt_para_coords;
                        for (int i = 0; i < gpts_dummy.size(); i++) { curveLoop_pt_coords.push_back( std::vector<double>() );
                            gmsh::model::getValue(0, gpts_dummy[i].second, pt_para_coords, curveLoop_pt_coords[curveLoop_pt_coords.size() - 1]); }
                    }
                
                    // See if any points are not on the edge of the bounding surface
                    bool allPtsOnEdgeOfBoundingSurface = true;
                    for (int i_pt = 0; i_pt < curveLoop_pt_coords.size(); i_pt++) {
                        if (compareValues(curveLoop_pt_coords[i_pt][alt_dim[0]], alt_val[0][0], tol) ||
                            compareValues(curveLoop_pt_coords[i_pt][alt_dim[0]], alt_val[0][1], tol) ||
                            compareValues(curveLoop_pt_coords[i_pt][alt_dim[1]], alt_val[1][0], tol) ||
                            compareValues(curveLoop_pt_coords[i_pt][alt_dim[1]], alt_val[1][1], tol)) { continue; }
                        allPtsOnEdgeOfBoundingSurface = false;
                        break;
                    }   
                    
                    // If not all points are on the edges of the bounding surface, define the curve loop/surface and move onto the next
                    if (!allPtsOnEdgeOfBoundingSurface) {
                        curveLoop_group_tags.push_back( gmsh::model::occ::addCurveLoop(curveLoops_arranged[i_clg][i_cl]) );
                        continue;
                    }
                    
                    // If all points are on the edges of the bounding surface, look at the normal of the surface of one of the lines. If all points of
                    // that surface lie on the bounding surface (i.e., the dim coordinates are val), do not create the curve 
                    // loop/surface. Keep in mind that some of the lines do not have a surface yet (i.e., the ones created between
                    // AR corners and points for closing the AR volumes). If all lines do not have a surface, define the curve loop/surface
                    std::vector<double> normal_vec;
                    for (int i_ln = 0; i_ln < glns.size(); i_ln++) {
                        // The gmsh line tag could be larger than what is avaiable in gln_gsurf_pairs if the line is one that was made to connect an AR corner point. In this case,
                        // the surface of the line has not been created yet. As such, we cannot use this line to evaluate whether its surface lies on the bounding surface
                        if (glns[i_ln] > gln_gsurf_pairs.get_max_ind()) { continue; }

                        // If the line is defined in gln_gsurf_pairs, get the possible gmsh surfaces for the line
                        std::vector<int> &gsurfs = gln_gsurf_pairs[glns[i_ln]]; assert (gsurfs.size() == 2); // There should be 2 surfaces for every line in the closed STL geometry
                        // Each surface belongs entirely to some AR volume(s) (at most, 2). Go through the surfaces and see which lies in the current AR volume
                        int N_norms_found = 0;
                        for (int i_gs = 0; i_gs < gsurfs.size(); i_gs++) {
                            // Get the possible volumes for the surface
                            std::vector<int> &vols = gsurf_vol_pairs[gsurfs[i_gs]];
                            // See if the surface's volume is the AR volume being analyzed. If so, save the surface normal
                            for (int i_vol = 0; i_vol < vols.size(); i_vol++) { if (vols[i_vol] == i_AR) { normal_vec = gsurf_normal_pairs[gsurfs[i_gs]]; N_norms_found += 1; } }
                        }
                        assert (N_norms_found == 1);
                        break;
                    }
                    
                    // If the normal vector was not defined in the previous section, it is probably because none of the gmsh lines belong to a surface; they were created to
                    // attach the AR corner points. If this is the case, create the curve loop and surface
                    if (normal_vec.size() == 0) {
                        curveLoop_group_tags.push_back( gmsh::model::occ::addCurveLoop(curveLoops_arranged[i_clg][i_cl]) );
                        continue;
                    }

                    // If the normal of the surface points directly in the dim direction, then the surface lies on the bounding surface, and the curve loop should not be created
                    if (compareValues(normal_vec[dim], 1.0, tol) || compareValues(normal_vec[dim], -1.0, tol)) { continue; }

                    // If the normal of the surface points does not point directly in the dim direction, then the surface does not lie on the bounding surface, and the curve loop should be created
                    curveLoop_group_tags.push_back( gmsh::model::occ::addCurveLoop(curveLoops_arranged[i_clg][i_cl]) );
                }

                if (curveLoop_group_tags.size() > 0) {
                    // Create the surface with the curve loop group tags
                    
                    if (surf_gID == 69744) {
                        curveLoop_group_tags[1] *= -1;
                    }

                    surf_gID = gmsh::model::occ::addPlaneSurface(curveLoop_group_tags);
                    cout << "surf_gID: " << surf_gID << endl;

                    // Add the AR volume indices to the corresponding surface in gsurf_vol_pairs
                    gsurf_vol_pairs.push_back( bndry_vol );
                    //gsurf_normal_pairs.push_back( surf_normals[i_surf] );
                    for (int i_v = 0; i_v < bndry_vol.size(); i_v++) { vol_gsurf_pairs[bndry_vol[i_v]].push_back( surf_gID ); }
                }
            }


        }
    }
    cout << endl << "    Complete!" << endl;







    
    // Synchronize the curve loops and surfaces that were created
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Synchronizing..." << endl;
    gmsh::model::occ::synchronize();
    cout << "    Complete!" << endl;





    //gmsh::write("my_model.geo_unrolled");
    //exit(1);







    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Creating gmsh volumes in each AR..." << endl;
    // Go through the AR and try to create a gmsh volume with the surfaces
    std::vector<std::pair<int, int>> gvol_gformat;
    int gvol, gvol_first = -1;
    for (int i_AR = 0; i_AR < vol_gsurf_pairs.size(); i_AR++) {
        cout << "\r    Considering AR  i_AR = " << i_AR << "/" << vol_gsurf_pairs.size() - 1 << "." << std::flush;
        
        // If there are no surfaces in i_AR, go to the next i_AR
        if (vol_gsurf_pairs[i_AR].size() == 0) { continue; }

        // Find the surface loops within the gmsh surface tags of each AR
        std::vector<std::vector<int>> surfaceLoops = findSurfaceLoops(vol_gsurf_pairs[i_AR]);
        

        // Create the gmsh surface loops and cut volumes from the surface loops found
        for (int i_sl = 0; i_sl < surfaceLoops.size(); i_sl++) {
            int gsurfLoop = gmsh::model::occ::addSurfaceLoop(surfaceLoops[i_sl]);
            gvol = gmsh::model::occ::addVolume({gsurfLoop});
            
            // Change gvol_first if it is the first gvol created in order to count the gmsh volumes created
            if (gvol_first == -1) { gvol_first = gvol; }

            // Add the gmsh volume tag to the output, which will be in gmsh format (similar to the gmsh::model::occ::cut() function)
            gvol_gformat.push_back( {3, gvol} );
        }
    }
    cout << endl << "    Created " << gvol_gformat[gvol_gformat.size() - 1].second - (gvol_first - 1) << " gmsh volumes." << endl;
    cout << "    Complete!" << endl;
    



    //std::vector<std::pair<int, int>> all_entities;
    //gmsh::model::getEntities(all_entities, 2);
    //cout << "all_entities.size(): " << all_entities.size() << endl;
    



    // Synchronize the curve loops and surfaces that were created
    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Synchronizing..." << endl;
    gmsh::model::occ::synchronize();
    cout << "    Complete!" << endl;

    
    

    



    //gmsh::write("my_model.geo_unrolled");
    //exit(1);

    cout << FILE_NAME << ": " << FUNCTION_NAME << ": Successfully imported and cut STL geometry! Exiting " << FUNCTION_NAME << "." << endl;
    
    return gvol_gformat;
}




*/


















/*
void getTriPointsFromSTL(const string &STL_file_path, const double &geo_scale,
    std::vector<std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>> &coord_ln_surf_pairs,
    std::vector<std::vector<double>> &surf_normals, int &tri_count)
{
    // Create the file stream
    ifstream file(STL_file_path);
    if (!file) { cerr << "gmsh_utils.cpp: getTriPointsFromSTL(): CRITICAL ERROR: Error opening STL file for reading. Check that the file path is correct." << endl; exit(1); }

    // Read the file, obtain the points, and give initial line/surface IDs
    cout << "gmsh_utils.cpp: getTriPointsFromSTL(): Reading STL file... " << endl;
    string line;
    tri_count = -1;
    while (getline(file, line))
    {
        if (line.substr(0, 8) == "endsolid") { break; }
        //if (line == "outer loop") { tri_count += 1; }
        if (line.substr(0, 12) == "facet normal") {
            tri_count += 1;
            surf_normals.push_back( std::vector<double>() );
            getCoordFromSTLLine(line, surf_normals[tri_count], 1.0, 12);
        }
        if (line == "endloop") {
            std::get<1>(coord_ln_surf_pairs[coord_ln_surf_pairs.size() - 3]) = {(int)coord_ln_surf_pairs.size() - 3, (int)coord_ln_surf_pairs.size() - 1};
            std::get<1>(coord_ln_surf_pairs[coord_ln_surf_pairs.size() - 2]) = {(int)coord_ln_surf_pairs.size() - 3, (int)coord_ln_surf_pairs.size() - 2};
            std::get<1>(coord_ln_surf_pairs[coord_ln_surf_pairs.size() - 1]) = {(int)coord_ln_surf_pairs.size() - 2, (int)coord_ln_surf_pairs.size() - 1};
            
            std::get<2>(coord_ln_surf_pairs[coord_ln_surf_pairs.size() - 3]) = {tri_count};
            std::get<2>(coord_ln_surf_pairs[coord_ln_surf_pairs.size() - 2]) = {tri_count};
            std::get<2>(coord_ln_surf_pairs[coord_ln_surf_pairs.size() - 1]) = {tri_count};
        }
        if (line.substr(0, 6) == "vertex") {
            coord_ln_surf_pairs.push_back( std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>() );
            getCoordFromSTLLine(line, std::get<0>(coord_ln_surf_pairs[coord_ln_surf_pairs.size() - 1]), geo_scale);
        }
    }
    tri_count += 1; // Add one more, since we started at -1 to make the ID start at 0
    cout << "gmsh_utils.cpp: getTriPointsFromSTL(): Complete!" << endl;
}

void getCoordFromSTLLine(const string &line, std::vector<double> &coord, const double &scale, const int pos)
{
    size_t first_space = line.find(' ', pos);
    size_t first_letter_pos = first_space + 1;
    size_t second_space = line.find(' ', first_letter_pos);
    size_t second_letter_pos = second_space + 1;
    size_t third_space = line.find(' ', second_letter_pos);
    size_t third_letter_pos = third_space + 1;
    size_t end_line = line.find('\n', third_letter_pos);
    double x_coord = stod( line.substr(first_letter_pos, second_space - first_letter_pos) );
    double y_coord = stod( line.substr(second_letter_pos, third_space - second_letter_pos) );
    double z_coord = stod( line.substr(third_letter_pos, end_line - third_letter_pos) );
    coord = {x_coord * scale, y_coord * scale, z_coord * scale};    
}
*/




// ========================================
// Functions for making/editting averaging region arrangments
// ========================================  
void MakeARMesh_Uniform2DRectangular(std::vector<int> N_AR, std::vector<double> ell, std::vector<double> L,
    std::vector<std::pair<int, int>> &ARs_uc, std::vector<std::vector<std::vector<std::vector<double>>>> &ARs_uc_geo,
    std::vector<std::vector<int>> &AR_neighbors)
{
    int last_ind = 0;
    for (int j = 0; j < N_AR[1]; ++j)
    {
        for (int i = 0; i < N_AR[0]; ++i)
        {
            // Define the bottom left corner of the AR
            double x = i * ell[0] - L[0]/2;
            double y = j * ell[1] - L[1]/2;

            // Create the AR
            int AR = gmsh::model::occ::addRectangle(x, y, 0.0, ell[0], ell[1]);
            
            // Synchronize the model
            gmsh::model::occ::synchronize();
            
            // Store the AR so that the geometry may be cut into it later
            ARs_uc.push_back({2, AR}); // (dim = 2 for surface)

            // Store the coordinates of the AR's points so that we can distinguish the AR interfaces from 
            // geometry surfaces after cutting the AR. ARs_uc_geo is:
            //      2D: ARs_uc_geo[AR, lines, points, coords]
            ARs_uc_geo.push_back(getSurfaceCoords(AR));


            // Define AR_neighbors- (not really needed... is it?)
            AR_neighbors.push_back( std::vector<int>() );
            if (N_AR[1] > 1) {
                if (j == 0) {
                    AR_neighbors[last_ind].push_back( i + (j + 1)*N_AR[0] ); }
                else if (j == N_AR[1] - 1) {
                    AR_neighbors[last_ind].push_back( i + (j - 1)*N_AR[0] ); }
                else {
                    AR_neighbors[last_ind].push_back( i + (j - 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i + (j + 1)*N_AR[0] ); } }
            
            if (N_AR[0] > 1) {
                if (i == 0) {
                    AR_neighbors[last_ind].push_back( (i + 1) + j*N_AR[0] ); }
                else if (i == N_AR[0] - 1) {
                    AR_neighbors[last_ind].push_back( (i - 1) + j*N_AR[0] ); }
                else {
                    AR_neighbors[last_ind].push_back( (i - 1) + j*N_AR[0] );
                    AR_neighbors[last_ind].push_back( (i + 1) + j*N_AR[0] ); } }

            // Consider Diagonals
            if (N_AR[1] > 1 && N_AR[0] > 1) {
                if (i == 0 && j == 0) {
                    AR_neighbors[last_ind].push_back( i + 1 + (j + 1)*N_AR[0] ); }
                else if (i == N_AR[0] - 1 && j == 0) {
                    AR_neighbors[last_ind].push_back( i - 1 + (j + 1)*N_AR[0] ); }
                else if (i == 0 && j == N_AR[1] - 1) {
                    AR_neighbors[last_ind].push_back( i + 1 + (j - 1)*N_AR[0] ); }
                else if (i == N_AR[0] - 1 && j == N_AR[1] - 1) {
                    AR_neighbors[last_ind].push_back( i - 1 + (j - 1)*N_AR[0] ); }
                else if (j == 0) {
                    AR_neighbors[last_ind].push_back( i - 1 + (j + 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i + 1 + (j + 1)*N_AR[0] ); }
                else if (j == N_AR[1] - 1) {
                    AR_neighbors[last_ind].push_back( i - 1 + (j - 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i + 1 + (j - 1)*N_AR[0] ); }
                else if (i == 0) {
                    AR_neighbors[last_ind].push_back( i + 1 + (j - 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i + 1 + (j + 1)*N_AR[0] ); }
                else if (i == N_AR[0] - 1) {
                    AR_neighbors[last_ind].push_back( i - 1 + (j - 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i - 1 + (j + 1)*N_AR[0] ); }
                else {
                    AR_neighbors[last_ind].push_back( i - 1 + (j - 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i + 1 + (j - 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i - 1 + (j + 1)*N_AR[0] );
                    AR_neighbors[last_ind].push_back( i + 1 + (j + 1)*N_AR[0] ); }
            }
            last_ind += 1;
        }
    }
}

void MakeARMesh_Uniform3DRectangular(const std::vector<int> &N_AR, const std::vector<double> &ell,
    const std::vector<double> &L, std::vector<std::pair<int, int>> &ARs_uc,
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &ARs_uc_geo)
{
    for (int k = 0; k < N_AR[2]; ++k)
    {
        for (int j = 0; j < N_AR[1]; ++j)
        {
            for (int i = 0; i < N_AR[0]; ++i)
            {
                // Define the bottom left corner of the AR
                double x = i * ell[0] - L[0]/2;
                double y = j * ell[1] - L[1]/2;
                double z = k * ell[2] - L[2]/2;

                // Create the AR
                int AR = gmsh::model::occ::addBox(x, y, z, ell[0], ell[1], ell[2]);
                
                // Store the AR so that the geometry may be cut into it later
                ARs_uc.push_back({3, AR});

                // Store the coordinates of the AR's points so that we can distinguish the AR interfaces from 
                // geometry surfaces after cutting the AR. ARs_uc_geo is:
                //      3D: ARs_uc_geo[AR, surface, lines, points, coords]
                //ARs_uc_geo.push_back(getVolumeCoords(AR));
            }
        }
    }
    
    // Synchronize the model
    gmsh::model::occ::synchronize();
}

/*
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
*/
// NOW USING THIS ONE
bool createMergedARGroups(const std::vector<std::pair<int, int>> &ARs,
    double min_vol_threshold, double max_AR_length, double max_AR_length_ratio,
    std::vector<std::vector<int>> &AR_tags_toBeMerged, std::vector<double> &AR_pore_space_vols)
{
    // In AR_tags_toBeMerged, the inner list is {merge host tag (i.e., the big AR), merging tag (i.e., the small AR)}
    bool anything_Merged = false;

    for (int i_ar = 0; i_ar < ARs.size(); i_ar++)
    {
        // Get AR information
        int dim = ARs[i_ar].first;
        int tag = ARs[i_ar].second;

        // Compute AR volume/area
        double AR_pore_space_vol;
        gmsh::model::occ::getMass(dim, tag, AR_pore_space_vol);

        // Decide if the AR must be merged (i.e., it is smaller than the threshold)
        if (AR_pore_space_vol < min_vol_threshold)
        {
            // Get the max/min x/y/z coordinates of the AR
            std::vector<double> min_max_x_y_z = getEntityMinMaxCoords(ARs[i_ar]); // min x, max x, min y, max y, min z, max z ---> min_max_x_y[0], min_max_x_y[1], min_max_x_y[2], min_max_x_y[3], min_max_x_y[4], min_max_x_y[5] 
            
            // See if a "host" AR can be found to merge the small AR to
            bool found_host = false;
            for (int j_ar = 0; j_ar < ARs.size(); j_ar++)
            {
                if (i_ar == j_ar) { continue; } // AR cannot host itself
                
                // Get the tag of the potential host AR
                int tag2 = ARs[j_ar].second;

                // If host AR, AR[j_ar], is not adjacent to small AR, AR[i_ar], try the next AR
                if (!areAdjacentEntities(dim, tag, tag2)) { continue; }
                
                // For now, we require the host AR to satisfy the min_vol_threshold itself; otherwise, try the next AR
                double ARj_vol;
                gmsh::model::occ::getMass(dim, tag2, ARj_vol);
                if (ARj_vol < min_vol_threshold) { continue; }


                // TODO: check AR_tags_toBeMerged to see if host AR has already been merged. Then, consider merged host AR
                //       as a potential host, not just AR[j_ar]


                // If potential host AR is adjacent and large enough, get the max/min x/y/z coordinates of the host AR
                std::vector<double> new_min_max_x_y_z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                std::vector<double> min_max_x_y_z2 = getEntityMinMaxCoords(ARs[j_ar]);
                new_min_max_x_y_z[0] = min_max_x_y_z[0] < min_max_x_y_z2[0] ? min_max_x_y_z[0] : min_max_x_y_z2[0];
                new_min_max_x_y_z[1] = min_max_x_y_z[1] > min_max_x_y_z2[1] ? min_max_x_y_z[1] : min_max_x_y_z2[1];
                new_min_max_x_y_z[2] = min_max_x_y_z[2] < min_max_x_y_z2[2] ? min_max_x_y_z[2] : min_max_x_y_z2[2];
                new_min_max_x_y_z[3] = min_max_x_y_z[3] > min_max_x_y_z2[3] ? min_max_x_y_z[3] : min_max_x_y_z2[3];
                new_min_max_x_y_z[4] = min_max_x_y_z[4] < min_max_x_y_z2[4] ? min_max_x_y_z[4] : min_max_x_y_z2[4];
                new_min_max_x_y_z[5] = min_max_x_y_z[5] > min_max_x_y_z2[5] ? min_max_x_y_z[5] : min_max_x_y_z2[5];
                

                // TODO: again, host AR might be a merged AR. Get the min/max x/y/z for the merged AR (i.e., AR_tags_toBeMerged)


                // If merging the host and small ARs would make the new AR too big, try the next AR
                if ( new_min_max_x_y_z[1] - new_min_max_x_y_z[0] > max_AR_length ||
                     new_min_max_x_y_z[3] - new_min_max_x_y_z[2] > max_AR_length ||
                     new_min_max_x_y_z[5] - new_min_max_x_y_z[4] > max_AR_length ) { continue; }
                
                // If merging the host and small ARs would create an AR with too large a length scale ratio, try the next AR
                //if (new_min_max_x_y_z[4] == 0.0 && new_min_max_x_y_z[5] == 0.0) { // For 2D
                if (dim == 2) { // For 2D
                    if ( (new_min_max_x_y_z[1] - new_min_max_x_y_z[0])/(new_min_max_x_y_z[3] - new_min_max_x_y_z[2]) >= max_AR_length_ratio ||
                         (new_min_max_x_y_z[3] - new_min_max_x_y_z[2])/(new_min_max_x_y_z[1] - new_min_max_x_y_z[0]) >= max_AR_length_ratio ) { continue; } }
                else { // For 3D
                    if ( (new_min_max_x_y_z[1] - new_min_max_x_y_z[0])/(new_min_max_x_y_z[3] - new_min_max_x_y_z[2]) >= max_AR_length_ratio ||
                         (new_min_max_x_y_z[3] - new_min_max_x_y_z[2])/(new_min_max_x_y_z[1] - new_min_max_x_y_z[0]) >= max_AR_length_ratio ||
                         (new_min_max_x_y_z[1] - new_min_max_x_y_z[0])/(new_min_max_x_y_z[5] - new_min_max_x_y_z[4]) >= max_AR_length_ratio ||
                         (new_min_max_x_y_z[5] - new_min_max_x_y_z[4])/(new_min_max_x_y_z[1] - new_min_max_x_y_z[0]) >= max_AR_length_ratio ||
                         (new_min_max_x_y_z[3] - new_min_max_x_y_z[2])/(new_min_max_x_y_z[5] - new_min_max_x_y_z[4]) >= max_AR_length_ratio ||
                         (new_min_max_x_y_z[5] - new_min_max_x_y_z[4])/(new_min_max_x_y_z[3] - new_min_max_x_y_z[2]) >= max_AR_length_ratio ) { continue; } }
                
                
                // If everything is satisfied, merge the host and small ARs
                bool new_entry = true;
                for (int i = 0; i < AR_tags_toBeMerged.size(); i++) {
                    if (AR_tags_toBeMerged[i][0] == tag2) {
                        AR_tags_toBeMerged[i].push_back( tag );
                        AR_pore_space_vols[i] += AR_pore_space_vol;
                        new_entry = false;
                        break; } }
                if (new_entry) { AR_tags_toBeMerged.push_back( {tag2, tag} ); AR_pore_space_vols.push_back( ARj_vol + AR_pore_space_vol ); }
                
                found_host = true;
                anything_Merged = true;
                break;
            }
            
            if (!found_host) { AR_tags_toBeMerged.push_back( {tag} ); AR_pore_space_vols.push_back( AR_pore_space_vol ); }

        }
        else
        {
            bool new_entry = true;
            for (int i = 0; i < AR_tags_toBeMerged.size(); i++) {
                if (AR_tags_toBeMerged[i][0] == tag) { new_entry = false; break; } }
            if (new_entry) { AR_tags_toBeMerged.push_back( {tag} ); AR_pore_space_vols.push_back( AR_pore_space_vol ); }
        }
    }

    return anything_Merged;
}







// ========================================
// Functions for separating averaging region lines from cut geometry lines 
// ========================================  
/*
// This funtion is no longer used. It's faster/more reliable to use getNonARInterfaceEntityTags.
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
*/

/*
// This funtion is no longer used. It's faster/more reliable to use getNonARInterfaceEntityTags.
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
*/

// Function that removes AR interface lines/surfaces from the resulting cut geometry
void getNonARInterfaceEntityTags(const std::vector<std::pair<int, int>> &ARs,
    std::vector<std::vector<int>> &entity_tags_in_AR, int AR_dim)
{
    // First, create a vector of pairs {AR, entity tag}. These will be used to find line/surface entities that
    // touch 2 ARs (i.e., this indicates the entity is an AR interface, not a geo cut entity)
    
    // -------For OMP implementation----------
    //int N_threads = omp_get_num_threads();
    //std::vector<std::vector<std::pair<int, int>>> AR_ln_pairs_OMP;
    //for (int i = 0; i < omp_get_num_threads(); i++) { AR_ln_pairs_OMP.push_back( std::vector<std::pair<int, int>>() ); }
    // Add stuff to AR_ln_pairs_OMP[thread_ID]. After "#pragma omp parallel", combine thread vectors into AR_ln_pairs.
    //#pragma omp parallel
    //{
    //    // Get the current thread's ID
    //    int thread_ID = omp_get_thread_num();
    // -------For OMP implementation----------
    
    // Go through the ARs
    std::vector<std::pair<int, int>> sorted_AR_entity_pairs, AR_entities;
    for (int i_ar = 0; i_ar < ARs.size(); i_ar++)
    {
        // Get the line/surface tags in each AR, and store them in AR_entities
        assert (ARs[i_ar].first == AR_dim);
        gmsh::model::getBoundary({ARs[i_ar]}, AR_entities, false, false, false);
        
        // Go through the lines/surfaces, and add them with their AR number to sorted_AR_entity_pairs
        for (int i_ar_ln = 0; i_ar_ln < AR_entities.size(); i_ar_ln++) { sorted_AR_entity_pairs.push_back( {i_ar, AR_entities[i_ar_ln].second} ); }
    }

    // Sort the {AR, entity tag} pairs based on the entity tag
    std::sort(sorted_AR_entity_pairs.begin(), sorted_AR_entity_pairs.end(), 
              [](const std::pair<int, int>& a, const std::pair<int, int>& b) { return a.second < b.second; });
    
    // Initialize entity_tags_in_AR
    entity_tags_in_AR.clear();
    for (int i = 0; i < ARs.size(); i++) { entity_tags_in_AR.push_back( std::vector<int>() ); }

    // Go through the sorted list of pairs
    bool skip = false;
    //std::vector<std::pair<int, int>> AR_interface_entities;
    for (int i = 0; i < sorted_AR_entity_pairs.size(); i++)
    {
        if (skip) { skip = false; continue; }
        
        // If you did not skip on the very last entity tag, then it is unique and not an AR interface
        if (i == sorted_AR_entity_pairs.size() - 1) { entity_tags_in_AR[sorted_AR_entity_pairs[i].first].push_back( sorted_AR_entity_pairs[i].second ); continue; }
        // Check if the entity tag matches the previous tag. If it does, something is wrong. Either an entity is duplicated 3 times, or "skip" is not working
        if (i != 0) { if (sorted_AR_entity_pairs[i].second == sorted_AR_entity_pairs[i-1].second) { cerr << "gmsh_utils.cpp: getNonARInterfaceEntityTags(): CRITICAL ERROR: Found an entity (line for 2D, surface for 3D) that touches 3 AR. Should not be possible." << endl; exit(1); } }
        
        // If the entity tag matches the next, then skip; it is an AR interface
        if (sorted_AR_entity_pairs[i].second == sorted_AR_entity_pairs[i+1].second) { skip = true; continue; } //AR_interface_entities.push_back( sorted_AR_entity_pairs[i] );
        // If the entity tag does not match the next (i.e., it is unique), then it is not an AR interface
        else { entity_tags_in_AR[sorted_AR_entity_pairs[i].first].push_back( sorted_AR_entity_pairs[i].second ); }
    }
}

void getPostCutTagsForLine(const std::vector<std::vector<double>> &target_ln_input,
    const std::vector<std::vector<std::vector<double>>> &cut_geo_ln_coords,
    const std::vector<int> &entity_tags_in_AR, const std::vector<double> &domain_boundaries,
    std::vector<int> &ln_tags, std::vector<int> &skip_ind)
{
    // Check to see if at least 1 of the target line points is within the domain (on the boundary is not good enough). If not, return.
    std::vector<std::vector<double>> target_ln;
    if (!(target_ln_input[0][0] > domain_boundaries[0] && target_ln_input[0][0] < domain_boundaries[1] &&
        target_ln_input[0][1] > domain_boundaries[2] && target_ln_input[0][1] < domain_boundaries[3]))
    {
        if (!(target_ln_input[1][0] > domain_boundaries[0] && target_ln_input[1][0] < domain_boundaries[1] &&
            target_ln_input[1][1] > domain_boundaries[2] && target_ln_input[1][1] < domain_boundaries[3])) { return; }
        target_ln = {target_ln_input[1], target_ln_input[0]};
    }
    else { target_ln = target_ln_input; }


    // Go through the lines post-cut lines to see if their endpoints match the target line's
    int second_pt_ind = -1;
    bool skip = false;
    for (int i_ln = 0; i_ln < cut_geo_ln_coords.size(); i_ln++) {
        
        // Check if the current line index should be skipped due to its previous analysis (i.e., it has already been identified as part of the target line)
        for (int i_skip = 0; i_skip < skip_ind.size(); i_skip++) { if (i_ln == skip_ind[i_skip]) { skip = true; break; } }
        if (skip) { continue; }
        
        // See if the current line cut_geo_ln_coords[i_ln] has an end point that matches the first target line end point. If not, go to the next i_ln
        if (compareVectors(target_ln[0], cut_geo_ln_coords[i_ln][0])) { second_pt_ind = 1; }
        else if (compareVectors(target_ln[0], cut_geo_ln_coords[i_ln][1])) { second_pt_ind = 0; }
        else { continue; }
        
        // See if the second point of cut_geo_ln_coords[i_ln] is the same as the second end point of the target line.
        // If so, cut_geo_ln_coords[i_ln] is the line; store the tag and return from the function
        if (compareVectors(target_ln[1], cut_geo_ln_coords[i_ln][second_pt_ind])) { ln_tags.push_back( entity_tags_in_AR[i_ln] ); return; }

        // If not, first check if the line in cut_geo_ln_coords is parallel with the target line. If not, it is not the correct line; move onto the next i_ln
        if (!areParallel(target_ln, cut_geo_ln_coords[i_ln])) { continue; }

        // If it is parallel, save the tag, as cut_geo_ln_coords[i_ln] is (probably) part of the target line. Add i_ln to the skip lines. Then, create a line
        // from the second point in cut_geo_ln_coords[i_ln] and the second target line point, and recall this function. (i.e., this will go through the lines
        // again and find a line that matches the second point in cut_geo_ln_coords[i_ln], and hopefully the second target line point, or at least find a line
        // that matches the first point and is parallel, in which case the function can be recursively recalled again until all line tags that make up the
        // target line are found)
        ln_tags.push_back( entity_tags_in_AR[i_ln] );
        skip_ind.push_back( i_ln );
        getPostCutTagsForLine( {cut_geo_ln_coords[i_ln][second_pt_ind], target_ln[1]}, cut_geo_ln_coords, entity_tags_in_AR, domain_boundaries, ln_tags, skip_ind);
        return;
    }
    
    // If the function gets to here, then the cut lines have run out, and the target line can not be found. This should not happen.
    cerr << "gmsh_utils.cpp: getNewTagsForCut(): CRITICAL ERROR: Entity geometry (target_ln) not found in cut geometry (cut_geo_ln_coords). This should not happen; a line that is not in the cut geometry is being looked for." << endl;
    exit(1);
}

void getToolGeoTags(const std::vector<std::vector<int>> &cut_inds,
    const std::vector<std::vector<int>> &entity_tags_in_AR,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo,
    std::vector<double> domain_boundaries,
    std::vector<std::vector<int>> &saved_tags)
{
    // Define necessary variables
    std::vector<int> entity_tags_in_AR_flat;
    std::vector<std::vector<std::vector<double>>> cut_geo_ln_coords, domain_boundary_lns;

    // Get the line coordinates of all post-cut lines (except for interface lines, since they are not in entity_tags_in_AR)
    for (int i_ar = 0; i_ar < entity_tags_in_AR.size(); i_ar++) {
        for (int i_tag = 0; i_tag < entity_tags_in_AR[i_ar].size(); i_tag++) {
            cut_geo_ln_coords.push_back( getLineCoords(entity_tags_in_AR[i_ar][i_tag]) );
            entity_tags_in_AR_flat.push_back( entity_tags_in_AR[i_ar][i_tag] ); } }
    
    // Get the domain boundary lines from the domain boundaries
    std::vector<double> pt_minmin = {domain_boundaries[0], domain_boundaries[2], 0.0};
    std::vector<double> pt_maxmin = {domain_boundaries[1], domain_boundaries[2], 0.0};
    std::vector<double> pt_minmax = {domain_boundaries[0], domain_boundaries[3], 0.0};
    std::vector<double> pt_maxmax = {domain_boundaries[1], domain_boundaries[3], 0.0};
    domain_boundary_lns.push_back( {pt_minmin, pt_maxmin} );
    domain_boundary_lns.push_back( {pt_maxmin, pt_maxmax} );
    domain_boundary_lns.push_back( {pt_maxmax, pt_minmax} );
    domain_boundary_lns.push_back( {pt_minmax, pt_minmin} );

    // Define necessary variables
    std::vector<std::vector<std::vector<double>>> lns_to_find;

    // Go through the groups of cut surfaces to be assigned to different physical groups
    for (int i_pg = 0; i_pg < cut_inds.size(); i_pg++) {// pg is physical group
        saved_tags.push_back( std::vector<int>() ); // We will save the line tags corresponding to each physical group in a different vector
        for (int i_ind = 0; i_ind < cut_inds[i_pg].size(); i_ind++) {
            lns_to_find = tool_geo[ cut_inds[i_pg][i_ind] ]; // These are the lines corresponding to a single pre-cut cut surf. Find the corresponding post-cut line tags for each one. There can be multiple pre-cut cut surfs assigned to make up a physical group

            for (int i_ln = 0; i_ln < lns_to_find.size(); i_ln++) {
                std::vector<int> skip_ind;

                // The following can turn into
                //     getPostCutTagsForLine(lns_to_find[i_ln], cut_geo_ln_coords, entity_tags_in_AR_flat, domain_boundaries, saved_tags[saved_tags.size() - 1], skip_ind); 
                // once the intersecting line/domain problem is figured out
                
                // Make sure at least 1 point of lns_to_find[i_ln] is within the area (i.e., domain) covered by ARs. Make sure that point is provided as the first point to getPostCutTagsForLine
                
                // Check the first point; if it is within the domain
                if (lns_to_find[i_ln][0][0] > domain_boundaries[0] && lns_to_find[i_ln][0][0] < domain_boundaries[1] &&
                    lns_to_find[i_ln][0][1] > domain_boundaries[2] && lns_to_find[i_ln][0][1] < domain_boundaries[3])
                {
                    // If the first point is within the domain, call getPostCutTagsForLine to find the line tags associated with lns_to_find[i_ln] (i.e., those that are within the domain)
                    getPostCutTagsForLine(lns_to_find[i_ln], cut_geo_ln_coords, entity_tags_in_AR_flat, domain_boundaries, saved_tags[saved_tags.size() - 1], skip_ind);
                }
                else
                {
                    // If the first point is not within the domain, check the second point
                    if (lns_to_find[i_ln][1][0] > domain_boundaries[0] && lns_to_find[i_ln][1][0] < domain_boundaries[1] &&
                        lns_to_find[i_ln][1][1] > domain_boundaries[2] && lns_to_find[i_ln][1][1] < domain_boundaries[3])
                    {
                        // If the second point is within the domain, call getPostCutTagsForLine with the points flipped to find the line tags associated with lns_to_find[i_ln] (i.e., those that are within the domain)
                        getPostCutTagsForLine({lns_to_find[i_ln][1], lns_to_find[i_ln][0]}, cut_geo_ln_coords, entity_tags_in_AR_flat, domain_boundaries, saved_tags[saved_tags.size() - 1], skip_ind);
                    }
                    else
                    {
                        // If the second point is also not within the domain, skip the line for now. Need to write more code for this
                        // If lns_to_find[i_ln] intersect the domain, then we should stop the code.
                        if (linesIntersect(lns_to_find[i_ln], domain_boundary_lns[0]) || linesIntersect(lns_to_find[i_ln], domain_boundary_lns[1]) ||
                            linesIntersect(lns_to_find[i_ln], domain_boundary_lns[2]) || linesIntersect(lns_to_find[i_ln], domain_boundary_lns[3]))
                        {
                            cerr << "gmsh_utils.cpp: getToolGeoTags(): CRITICAL ERROR: The end points of tool cut line (lns_to_find[i_ln]) exist outside the domain, but intersect the domain. This might cause error in the labeling of boundaries. Further code should be written for this case." << endl;
                            exit(1);
                        }
                        continue;
                    }
                    
                }
            }

        }
    }
}

void getDomainBoundaryLineTags(const std::vector<std::vector<int>> &ln_tags_in_AR,
    std::vector<std::vector<int>> &domain_boundary_tags, const std::vector<std::pair<int, int>> &ARs,
    std::vector<std::vector<int>> &AR_tags_neighbors, double L_x, double L_y, double tol)
{
    bool verbose = false;
    int tag;
    std::vector<int> top_ln_tags, right_ln_tags, bottom_ln_tags, left_ln_tags, unused_ln_tags;
    std::vector<int> top_AR_tags, right_AR_tags, bottom_AR_tags, left_AR_tags;
    
    for (int i_ar = 0; i_ar < ln_tags_in_AR.size(); i_ar++) {
        for (int i_ln = 0; i_ln < ln_tags_in_AR[i_ar].size(); i_ln++) {
            
            tag = ln_tags_in_AR[i_ar][i_ln];
        
            // Check if the line lies on the inlet (i.e., the left-most side of the domain).
            if (isCoincides(tag, 0, -L_x/2, tol)) { left_ln_tags.push_back( tag ); left_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added curve " << tag << " to the left lines vector." << std::endl; }
            // Check if the line lies on the outlet (i.e., the right-most side of the domain).
            else if (isCoincides(tag, 0, L_x/2, tol)) { right_ln_tags.push_back( tag ); right_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added curve " << tag << " to the right lines vector." << std::endl; }
            // Check if the line lies on the top boundary
            else if (isCoincides(tag, 1, L_y/2, tol)) { top_ln_tags.push_back( tag ); top_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added curve " << tag << " to the top lines vector." << std::endl; }
            // Check if the line lies on the bottom boundary
            else if (isCoincides(tag, 1, -L_y/2, tol)) { bottom_ln_tags.push_back( tag ); bottom_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added curve " << tag << " to the bottom lines vector." << std::endl; }
            // If none of the previous, put it in the unused tags. These are the cut geometry line tags.
            else { unused_ln_tags.push_back( tag ); if (verbose) std::cout << "Added curve " << tag << " to the unused lines vector." << std::endl; }
        }
    }

    domain_boundary_tags.push_back( top_ln_tags );
    domain_boundary_tags.push_back( right_ln_tags );
    domain_boundary_tags.push_back( bottom_ln_tags );
    domain_boundary_tags.push_back( left_ln_tags );
    domain_boundary_tags.push_back( unused_ln_tags );

    AR_tags_neighbors.push_back( top_AR_tags );
    AR_tags_neighbors.push_back( right_AR_tags );
    AR_tags_neighbors.push_back( bottom_AR_tags );
    AR_tags_neighbors.push_back( left_AR_tags );
}

void getDomainBoundarySurfaceTags(const std::vector<std::vector<int>> &surf_tags_in_AR,
    std::vector<std::vector<int>> &domain_boundary_tags, const std::vector<double> &L)
{
    bool verbose = false;
    int tag;
    std::vector<int> top_ln_tags, bottom_ln_tags, left_ln_tags, right_ln_tags, front_ln_tags, back_ln_tags, unused_ln_tags;
    for (int i_ar = 0; i_ar < surf_tags_in_AR.size(); i_ar++)
    {
        // Create a physical group from the AR
        //AR_tags.push_back( gmsh::model::addPhysicalGroup(dim, {tag}) );
        for (int i_srfc = 0; i_srfc < surf_tags_in_AR[i_ar].size(); i_srfc++)
        {
            tag = surf_tags_in_AR[i_ar][i_srfc];

            // Check if the surface lies on the inlet (i.e., the left-most side of the domain).
            if (isCoincidesSurface(tag, 0, -L[0]/2)) { left_ln_tags.push_back( tag ); if (verbose) std::cout << "Added surface " << tag << " to the left surfaces vector." << std::endl; }
            // Check if the surface lies on the outlet (i.e., the right-most side of the domain).
            else if (isCoincidesSurface(tag, 0, L[0]/2)) { right_ln_tags.push_back( tag ); if (verbose) std::cout << "Added surface " << tag << " to the right surfaces vector." << std::endl; }
            // Check if the surface lies on the top boundary
            else if (isCoincidesSurface(tag, 1, L[1]/2)) { top_ln_tags.push_back( tag ); if (verbose) std::cout << "Added surface " << tag << " to the top surfaces vector." << std::endl; }
            // Check if the surface lies on the bottom boundary
            else if (isCoincidesSurface(tag, 1, -L[1]/2)) { bottom_ln_tags.push_back( tag ); if (verbose) std::cout << "Added surface " << tag << " to the bottom surfaces vector." << std::endl; }
            // Check if the surface lies on the front boundary
            else if (isCoincidesSurface(tag, 2, L[2]/2)) { front_ln_tags.push_back( tag ); if (verbose) std::cout << "Added surface " << tag << " to the front surfaces vector." << std::endl; }
            // Check if the surface lies on the back boundary
            else if (isCoincidesSurface(tag, 2, -L[2]/2)) { back_ln_tags.push_back( tag ); if (verbose) std::cout << "Added surface " << tag << " to the back surfaces vector." << std::endl; }
            // If none of the previous, put it in the unused tags
            else { unused_ln_tags.push_back( tag ); if (verbose) std::cout << "Added surface " << tag << " to the unused surfaces vector." << std::endl; }
        }
    }

    domain_boundary_tags.push_back( top_ln_tags );
    domain_boundary_tags.push_back( bottom_ln_tags );
    domain_boundary_tags.push_back( left_ln_tags );
    domain_boundary_tags.push_back( right_ln_tags );
    domain_boundary_tags.push_back( front_ln_tags );
    domain_boundary_tags.push_back( back_ln_tags );
    domain_boundary_tags.push_back( unused_ln_tags );
}
