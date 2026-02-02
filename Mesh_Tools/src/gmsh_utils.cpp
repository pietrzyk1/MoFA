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

// Define a function that takes in a set of entity tags (as {dim, tag} pairs) and returns the minimum/maximum x/y/z coordinates of the set's boundary points.
std::vector<double> getEntityMinMaxCoords(const std::vector<std::pair<int, int>> &dimtag)
{
    std::vector<double> minmaxcoords;
    int dim;
    for (int i = 0; i < dimtag.size(); i++) {
        // If i == 0, get the dimension of the provided entities. If i != 0, assert all entities have the same dimension
        if (i == 0) { dim = dimtag[i].first; } else { assert (dimtag[i].first == dim); }
        
        // Get the min/max coords for the entity
        std::vector<double> minmaxcoords_ent = getEntityMinMaxCoords(dimtag[i]); // output: {min x, max x, min y, max y, min z, max z}
        
        // Compare the min/max coords of the entity to the saved coords
        if (i == 0) { minmaxcoords = minmaxcoords_ent; }
        else {
            for (int j = 0; j < minmaxcoords_ent.size(); j++) {
                if (j % 2 == 0) { if (minmaxcoords_ent[j] < minmaxcoords[j]) { minmaxcoords[j] = minmaxcoords_ent[j]; } }
                else { if (minmaxcoords_ent[j] > minmaxcoords[j]) { minmaxcoords[j] = minmaxcoords_ent[j]; } }
            }
        }
    }

    return minmaxcoords;
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

// NOW USING THIS ONE. NOTE: this depends on "gmsh::model::occ::getMass", which needs properly oriented geometry (i.e., normals point out of volumes)
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
        AR_pore_space_vol = std::abs(AR_pore_space_vol);


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
                ARj_vol = std::abs(ARj_vol);
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






bool createMergedARGroups2(const std::vector<std::pair<int, int>> &ARs,
    double min_vol_threshold, double max_AR_length, double max_AR_length_ratio,
    std::vector<std::vector<int>> &AR_tags_toBeMerged, std::vector<double> &AR_pore_space_vols)
{
    // In AR_tags_toBeMerged, the inner list is {merge host tag (i.e., the big AR), merging tag (i.e., the small AR)}
    bool anything_Merged = false;
    int dim;


    // First, compute the area/volume of each AR and reorder them so that the AR with the smallest area/volume is considered first
    std::vector<std::pair<int, double>> ARs_reordered_mass(ARs.size());
    for (int i_AR = 0; i_AR < ARs.size(); i_AR++) {
        // Get AR information
        if (i_AR == 0) { dim = ARs[i_AR].first; } else { assert (ARs[i_AR].first == dim); };
        int tag = ARs[i_AR].second;

        // Compute AR area/volume
        double AR_pore_space_vol;
        gmsh::model::occ::getMass(dim, tag, AR_pore_space_vol);
        assert (AR_pore_space_vol > 0.0); // If the mass is < 0, the normals of the AR region bounding lines/surfaces are not all pointed outwards
        
        // Assign the tag-mass pair
        ARs_reordered_mass[i_AR] = {ARs[i_AR].second, AR_pore_space_vol};
    }

    // Sort the AR tags by their mass (smallest to largest)
    std::sort(ARs_reordered_mass.begin(), ARs_reordered_mass.end(),
        [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
        return a.second < b.second; });
    
    
    // Initiate the merge group indices. Given an index corresponding to ARs_reordered_mass, provide the index of AR_tags_toBeMerged where the AR is associated
    std::vector<int> MG_inds(ARs.size()); for (int i = 0; i < MG_inds.size(); i++) { MG_inds[i] = -1; }


    // Go through the AR and merge them to create a set of AR with a larger minimum AR area/volume than what was obtained with arbitrarily defining a grid of AR (i.e., no tiny ARs)
    int last_ind = -1;
    for (int i_ar = 0; i_ar < ARs_reordered_mass.size(); i_ar++)
    {
        // Get AR information
        int tag = ARs_reordered_mass[i_ar].first;
        double smallAR_poreVol = ARs_reordered_mass[i_ar].second;
        
        // Decide if the AR must be merged (i.e., it is smaller than the threshold)
        if (smallAR_poreVol < min_vol_threshold && MG_inds[i_ar] == -1)
        {
            // Get the max/min x/y/z coordinates of the AR
            std::vector<double> min_max_x_y_z = getEntityMinMaxCoords({dim, tag}); // min x, max x, min y, max y, min z, max z ---> min_max_x_y[0], min_max_x_y[1], min_max_x_y[2], min_max_x_y[3], min_max_x_y[4], min_max_x_y[5] 
            
            // See if a "host" AR can be found to merge the small AR to
            bool found_host = false;
            for (int j_ar = 0; j_ar < ARs_reordered_mass.size(); j_ar++)
            {
                // AR cannot host itself
                if (i_ar == j_ar) { continue; }
                
                // Get the tag of the potential host AR
                int tag2 = ARs_reordered_mass[j_ar].first;

                // If host AR, AR[j_ar], is not adjacent to small AR, AR[i_ar], try the next AR
                if (!areAdjacentEntities(dim, tag, tag2)) { continue; }
                
                
                // If the host AR is adjacent to the small AR, get the min/max coords for the host AR/host AR's merge group
                std::vector<double> min_max_x_y_z2;
                if (MG_inds[j_ar] != -1) {
                    int N_AR_MG = AR_tags_toBeMerged[MG_inds[j_ar]].size();
                    std::vector<std::pair<int, int>> mergeGroupDimTag(N_AR_MG);
                    for (int i_m = 0; i_m < N_AR_MG; i_m++) { mergeGroupDimTag[i_m] = {dim, AR_tags_toBeMerged[MG_inds[j_ar]][i_m]}; }
                    min_max_x_y_z2 = getEntityMinMaxCoords(mergeGroupDimTag);
                }
                else { min_max_x_y_z2 = getEntityMinMaxCoords({dim, tag2}); }

                // If the small AR were to be merged with the potential host, get the max/min x/y/z coordinates that would exist with the merge
                std::vector<double> new_min_max_x_y_z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                new_min_max_x_y_z[0] = min_max_x_y_z[0] < min_max_x_y_z2[0] ? min_max_x_y_z[0] : min_max_x_y_z2[0];
                new_min_max_x_y_z[1] = min_max_x_y_z[1] > min_max_x_y_z2[1] ? min_max_x_y_z[1] : min_max_x_y_z2[1];
                new_min_max_x_y_z[2] = min_max_x_y_z[2] < min_max_x_y_z2[2] ? min_max_x_y_z[2] : min_max_x_y_z2[2];
                new_min_max_x_y_z[3] = min_max_x_y_z[3] > min_max_x_y_z2[3] ? min_max_x_y_z[3] : min_max_x_y_z2[3];
                new_min_max_x_y_z[4] = min_max_x_y_z[4] < min_max_x_y_z2[4] ? min_max_x_y_z[4] : min_max_x_y_z2[4];
                new_min_max_x_y_z[5] = min_max_x_y_z[5] > min_max_x_y_z2[5] ? min_max_x_y_z[5] : min_max_x_y_z2[5];
                
                // If the merge would make the length(s) of the merge group (i.e., the new AR) too big, move on to checking the next potential host AR
                if ( new_min_max_x_y_z[1] - new_min_max_x_y_z[0] > max_AR_length ||
                     new_min_max_x_y_z[3] - new_min_max_x_y_z[2] > max_AR_length ||
                     new_min_max_x_y_z[5] - new_min_max_x_y_z[4] > max_AR_length ) { continue; }
                
                // If the merge would create a merge group (i.e., a new AR) with too large a length scale ratio, move on to checking the next potential host AR
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
                

                // If everything is satisfied, merge the small AR with the host AR/merge group
                bool new_entry = true;
                if (MG_inds[j_ar] != -1) {
                    AR_tags_toBeMerged[MG_inds[j_ar]].push_back( tag );
                    AR_pore_space_vols[MG_inds[j_ar]] += smallAR_poreVol;
                    MG_inds[i_ar] = MG_inds[j_ar];
                }
                else {
                    AR_tags_toBeMerged.push_back( {tag2, tag} );
                    AR_pore_space_vols.push_back( smallAR_poreVol + ARs_reordered_mass[j_ar].second );
                    MG_inds[i_ar] = AR_tags_toBeMerged.size() - 1;
                    MG_inds[j_ar] = AR_tags_toBeMerged.size() - 1;
                }

                found_host = true;
                anything_Merged = true;
                break;
            }
            
            if (!found_host) { AR_tags_toBeMerged.push_back( {tag} ); AR_pore_space_vols.push_back( smallAR_poreVol ); }
        }
        else
        {
            if (smallAR_poreVol > min_vol_threshold) { last_ind = i_ar; break; }
        }
    }

    // Add the remaining AR (i.e., the AR large enough to not need to merge), only if they were not used in a merge
    if (last_ind != -1) {
        for (int i_NM = last_ind; i_NM < ARs_reordered_mass.size(); i_NM++) {
            if (MG_inds[i_NM] == -1) {
                AR_tags_toBeMerged.push_back( {ARs_reordered_mass[i_NM].first} );
                AR_pore_space_vols.push_back( ARs_reordered_mass[i_NM].second );
                MG_inds[i_NM] = AR_tags_toBeMerged.size() - 1;
            }
        }
    }

    return anything_Merged;
}
*/

bool createMergedARGroups(const std::vector<std::pair<int, int>> &ARs,
    double min_vol_threshold, double max_AR_length, double max_AR_length_ratio,
    std::vector<std::vector<int>> &AR_tags_toBeMerged, std::vector<double> &AR_pore_space_vols)
{
    // Initialize variables
    int max_iterations = 20; // Maximum number of iterations for merging the AR/AR groups
    int dim; // Dimension of the ARs (i.e., 2 or 3 for 2D or 3D)
    int N_AR_init = ARs.size();

    
    // Fill the initial AR_mergeGroups_mass_pairs vector with the indices of the ARs (as if they were their
    // own merge groups, e.g. {{1}, {2}, ...}) and the areas/volumes of each AR
    std::vector<std::pair<std::vector<int>, double>> AR_mergeGroups_mass_pairs(N_AR_init);
    for (int i_AR = 0; i_AR < N_AR_init; i_AR++) {
        // Get AR information
        if (i_AR == 0) { dim = ARs[i_AR].first; } else { assert (ARs[i_AR].first == dim); };
        int tag = ARs[i_AR].second;

        // Compute AR area/volume
        double AR_pore_space_vol;
        gmsh::model::occ::getMass(dim, tag, AR_pore_space_vol);
        assert (AR_pore_space_vol > 0.0); // If the mass is < 0, the normals of the AR region bounding lines/surfaces are not all pointed outwards
        
        // Assign the tag-mass pair
        AR_mergeGroups_mass_pairs[i_AR] = {{ARs[i_AR].second}, AR_pore_space_vol};
    }

    // Merge ARs into groups according to the provided min_vol_threshold, max_AR_length, and max_AR_length_ratio.
    // Iterations will then merge the AR merge groups until the maximum number of iterations is reached, or if
    // no more merges can be made
    bool ARsWereMerged = false;
    int total_iters = 0;
    for (int i_iter = 0; i_iter < max_iterations; i_iter++) {
        total_iters = i_iter + 1;
        cout << "\r    createMergedARGroups(): Merge iteration " << total_iters << " (out of " << max_iterations << " maximum iterations)." << std::flush;
        bool anything_Merged_iter = createMergedARGroups_core(dim, AR_mergeGroups_mass_pairs, min_vol_threshold, max_AR_length, max_AR_length_ratio);
        if (!anything_Merged_iter) { break; }
        ARsWereMerged = true;
    }
    cout << "\r    createMergedARGroups(): Total merge iterations: " << total_iters << " (out of " << max_iterations << " maximum iterations)." << std::flush;
    cout << endl;
    
    // Extract the final AR merge groups and volumes
    AR_tags_toBeMerged.clear(); AR_tags_toBeMerged.resize(AR_mergeGroups_mass_pairs.size());
    AR_pore_space_vols.clear(); AR_pore_space_vols.resize(AR_mergeGroups_mass_pairs.size());    
    for (int i = 0; i < AR_mergeGroups_mass_pairs.size(); i++) {
        AR_tags_toBeMerged[i] = AR_mergeGroups_mass_pairs[i].first;
        AR_pore_space_vols[i] = AR_mergeGroups_mass_pairs[i].second;
    }
    
    return ARsWereMerged;
}

bool createMergedARGroups_core(int dim, std::vector<std::pair<std::vector<int>, double>> &AR_mergeGroups_mass_pairs,
    double min_vol_threshold, double max_AR_length, double max_AR_length_ratio)
{
    // Define a bool of whether anything was merged or not
    bool anything_Merged = false;

    // Define the initial number of AR (merge) groups
    int N_MG = AR_mergeGroups_mass_pairs.size();
    
    // Declare vectors for the input AR merge groups and AR volumes
    std::vector<std::vector<int>> AR_unmerged(N_MG);
    std::vector<double> AR_unmerged_vol(N_MG);
    
    // Initiate a vector that holds each AR's/initial merge group's new merge group index. In other words, given
    // an index corresponding to AR_unmerged, provide the index of AR_tags_toBeMerged where the AR is associated
    std::vector<int> MG_inds(N_MG); for (int i = 0; i < N_MG; i++) { MG_inds[i] = -1; }

    // Declare vectors for the output AR merge groups and AR volumes
    std::vector<std::vector<int>> AR_tags_toBeMerged;
    std::vector<double> AR_pore_space_vols;


    // Sort the input AR/merge group tags by their masses (smallest to largest)
    std::sort(AR_mergeGroups_mass_pairs.begin(), AR_mergeGroups_mass_pairs.end(),
        [](const std::pair<std::vector<int>, double> &a, const std::pair<std::vector<int>, double> &b) {
        return a.second < b.second; });
        
    // Extract the ordered AR merge IDs and volumes into AR_unmerged and AR_unmerged_vol. Then, clear AR_mergeGroups_mass_pairs
    // so that the solution can be stored in it later
    for (int i = 0; i < N_MG; i++) {
        AR_unmerged[i] = AR_mergeGroups_mass_pairs[i].first;
        AR_unmerged_vol[i] = AR_mergeGroups_mass_pairs[i].second;
    }
    AR_mergeGroups_mass_pairs.clear();


    // Go through the AR/merge groups and merge them to create a set of new merge groups with a larger minimum AR area/volume than what was provided
    int last_ind = -1;
    for (int i_ar = 0; i_ar < N_MG; i_ar++)
    {
        // Decide if the AR/merge group must be merged to another AR/merge group (i.e., its volume is smaller than the threshold)
        if (AR_unmerged_vol[i_ar] < min_vol_threshold && MG_inds[i_ar] == -1)
        {
            // Get the max/min x/y/z coordinates of the small AR/ merge group
            std::vector<double> min_max_x_y_z; // min x, max x, min y, max y, min z, max z ---> min_max_x_y[0], min_max_x_y[1], min_max_x_y[2], min_max_x_y[3], min_max_x_y[4], min_max_x_y[5] 
            {
                std::vector<std::pair<int, int>> mergeGroupDimTag(AR_unmerged[i_ar].size());
                for (int i_m = 0; i_m < AR_unmerged[i_ar].size(); i_m++) { mergeGroupDimTag[i_m] = {dim, AR_unmerged[i_ar][i_m]}; }
                min_max_x_y_z = getEntityMinMaxCoords(mergeGroupDimTag);
            }
            
            // See if a "host" AR/merge group can be found to merge the small AR/merge group to
            bool found_host = false;
            for (int j_ar = 0; j_ar < N_MG; j_ar++)
            {
                // The small AR/merge group cannot host itself
                if (i_ar == j_ar) { continue; }
                
                // If host AR/merge group, AR_unmerged[j_ar], is not adjacent to the small AR/merge group, AR_unmerged[i_ar], try the next AR/merge group
                bool areAdjacent = false;
                for (int i_adj = 0; i_adj < AR_unmerged[i_ar].size(); i_adj++) {
                    for (int j_adj = 0; j_adj < AR_unmerged[j_ar].size(); j_adj++) {
                        if (areAdjacentEntities(dim, AR_unmerged[i_ar][i_adj], AR_unmerged[j_ar][j_adj])) { areAdjacent = true; break; } }
                    if (areAdjacent) { break; } }
                if (!areAdjacent) { continue; }
                
                // If the host AR/merge group is adjacent to the small AR/merge group, get the min/max coords for the host AR/merge group
                std::vector<double> min_max_x_y_z2;
                {
                    std::vector<std::pair<int, int>> mergeGroupDimTag(AR_unmerged[j_ar].size());
                    for (int i_m = 0; i_m < AR_unmerged[j_ar].size(); i_m++) { mergeGroupDimTag[i_m] = {dim, AR_unmerged[j_ar][i_m]}; }
                    min_max_x_y_z2 = getEntityMinMaxCoords(mergeGroupDimTag);
                }

                // If the small AR/merge group were to be merged with the potential host, get the max/min x/y/z coordinates that would exist with the merge
                std::vector<double> new_min_max_x_y_z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                new_min_max_x_y_z[0] = min_max_x_y_z[0] < min_max_x_y_z2[0] ? min_max_x_y_z[0] : min_max_x_y_z2[0];
                new_min_max_x_y_z[1] = min_max_x_y_z[1] > min_max_x_y_z2[1] ? min_max_x_y_z[1] : min_max_x_y_z2[1];
                new_min_max_x_y_z[2] = min_max_x_y_z[2] < min_max_x_y_z2[2] ? min_max_x_y_z[2] : min_max_x_y_z2[2];
                new_min_max_x_y_z[3] = min_max_x_y_z[3] > min_max_x_y_z2[3] ? min_max_x_y_z[3] : min_max_x_y_z2[3];
                new_min_max_x_y_z[4] = min_max_x_y_z[4] < min_max_x_y_z2[4] ? min_max_x_y_z[4] : min_max_x_y_z2[4];
                new_min_max_x_y_z[5] = min_max_x_y_z[5] > min_max_x_y_z2[5] ? min_max_x_y_z[5] : min_max_x_y_z2[5];

                // If the merge would make the length(s) of the new merge group too big, move on to checking the next potential host AR/merge group
                if ( new_min_max_x_y_z[1] - new_min_max_x_y_z[0] > max_AR_length ||
                     new_min_max_x_y_z[3] - new_min_max_x_y_z[2] > max_AR_length ||
                     new_min_max_x_y_z[5] - new_min_max_x_y_z[4] > max_AR_length ) { continue; }
                
                // If the merge would create a new merge group with too large a length scale ratio, move on to checking the next potential host AR/merge group
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
                

                // If everything is satisfied, merge the small AR/merge group with the host AR/merge group
                if (MG_inds[j_ar] != -1) {
                    std::vector<int> &host = AR_tags_toBeMerged[MG_inds[j_ar]];
                    host.insert(host.end(), AR_unmerged[i_ar].begin(), AR_unmerged[i_ar].end());
                    AR_pore_space_vols[MG_inds[j_ar]] += AR_unmerged_vol[i_ar];
                    MG_inds[i_ar] = MG_inds[j_ar];
                }
                else {
                    std::vector<int> new_merge = AR_unmerged[j_ar]; new_merge.insert(new_merge.end(), AR_unmerged[i_ar].begin(), AR_unmerged[i_ar].end());
                    AR_tags_toBeMerged.push_back( new_merge );
                    AR_pore_space_vols.push_back( AR_unmerged_vol[i_ar] + AR_unmerged_vol[j_ar] );
                    MG_inds[i_ar] = AR_tags_toBeMerged.size() - 1;
                    MG_inds[j_ar] = AR_tags_toBeMerged.size() - 1;
                }

                found_host = true;
                anything_Merged = true;
                break;
            }
            
            if (!found_host) { AR_tags_toBeMerged.push_back( AR_unmerged[i_ar] ); AR_pore_space_vols.push_back( AR_unmerged_vol[i_ar] ); }
        }
        else
        {
            if (AR_unmerged_vol[i_ar] > min_vol_threshold) { last_ind = i_ar; break; }
        }
    }

    // Add the remaining ARs/merge groups (i.e., the ARs/merge groups large enough to not need to merge), only if they were not used in a merge
    if (last_ind != -1) {
        for (int i_NM = last_ind; i_NM < N_MG; i_NM++) {
            if (MG_inds[i_NM] == -1) {
                AR_tags_toBeMerged.push_back( AR_unmerged[i_NM] );
                AR_pore_space_vols.push_back( AR_unmerged_vol[i_NM] );
                MG_inds[i_NM] = AR_tags_toBeMerged.size() - 1;
            }
        }
    }

    // Finally, fill AR_mergeGroups_mass_pairs with the new merge groups and volumes
    assert (AR_tags_toBeMerged.size() == AR_pore_space_vols.size());
    AR_mergeGroups_mass_pairs.resize(AR_tags_toBeMerged.size());
    for (int i = 0; i < AR_tags_toBeMerged.size(); i++) { AR_mergeGroups_mass_pairs[i] = {AR_tags_toBeMerged[i], AR_pore_space_vols[i]}; }

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
    std::vector<std::vector<int>> &domain_boundary_tags, const std::vector<std::pair<int, int>> &ARs,
    std::vector<std::vector<int>> &AR_tags_neighbors, const std::vector<double> &L)
{
    bool verbose = false;
    int tag;
    std::vector<int> top_ln_tags, bottom_ln_tags, left_ln_tags, right_ln_tags, front_ln_tags, back_ln_tags, unused_ln_tags;
    std::vector<int> top_AR_tags, right_AR_tags, bottom_AR_tags, left_AR_tags, front_AR_tags, back_AR_tags;

    for (int i_ar = 0; i_ar < surf_tags_in_AR.size(); i_ar++)
    {
        // Create a physical group from the AR
        //AR_tags.push_back( gmsh::model::addPhysicalGroup(dim, {tag}) );
        for (int i_srfc = 0; i_srfc < surf_tags_in_AR[i_ar].size(); i_srfc++)
        {
            tag = surf_tags_in_AR[i_ar][i_srfc];

            // Check if the surface lies on the inlet (i.e., the left-most side of the domain).
            if (isCoincidesSurface(tag, 0, -L[0]/2)) { left_ln_tags.push_back( tag ); left_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added surface " << tag << " to the left surfaces vector." << std::endl; }
            // Check if the surface lies on the outlet (i.e., the right-most side of the domain).
            else if (isCoincidesSurface(tag, 0, L[0]/2)) { right_ln_tags.push_back( tag ); right_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added surface " << tag << " to the right surfaces vector." << std::endl; }
            // Check if the surface lies on the top boundary
            else if (isCoincidesSurface(tag, 1, L[1]/2)) { top_ln_tags.push_back( tag ); top_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added surface " << tag << " to the top surfaces vector." << std::endl; }
            // Check if the surface lies on the bottom boundary
            else if (isCoincidesSurface(tag, 1, -L[1]/2)) { bottom_ln_tags.push_back( tag ); bottom_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added surface " << tag << " to the bottom surfaces vector." << std::endl; }
            // Check if the surface lies on the front boundary
            else if (isCoincidesSurface(tag, 2, L[2]/2)) { front_ln_tags.push_back( tag ); front_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added surface " << tag << " to the front surfaces vector." << std::endl; }
            // Check if the surface lies on the back boundary
            else if (isCoincidesSurface(tag, 2, -L[2]/2)) { back_ln_tags.push_back( tag ); back_AR_tags.push_back( ARs[i_ar].second ); if (verbose) std::cout << "Added surface " << tag << " to the back surfaces vector." << std::endl; }
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
    
    AR_tags_neighbors.push_back( top_AR_tags );
    AR_tags_neighbors.push_back( right_AR_tags );
    AR_tags_neighbors.push_back( bottom_AR_tags );
    AR_tags_neighbors.push_back( left_AR_tags );
    AR_tags_neighbors.push_back( front_AR_tags );
    AR_tags_neighbors.push_back( back_AR_tags );
}
