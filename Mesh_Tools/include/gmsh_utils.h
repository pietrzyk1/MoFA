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
#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include <gmsh.h>
//#include "mfem.hpp"
#include "JSON_IO.h"





// ========================================
// Declare global variables
// ========================================  

// Declare the global tolerance for the gmsh functions
constexpr double tol_global = 1e-6;
// Define Pi
const double pi = 4.0 * std::atan(1.0);





// ========================================
// Declare functions
// ========================================

std::vector<double> getPointCoords(int pt_tag);
// Here:
//      pt_tag: a point tag

std::vector<std::vector<double>> getLineCoords(int ln_tag);
// Here:
//      ln_tag: a line tag

std::vector<std::vector<std::vector<double>>> getSurfaceCoords(int srfc_tag);
// Here:
//      srfc_tag: a surface tag

std::vector<std::vector<std::vector<std::vector<double>>>> getVolumeCoords(int vol_tag);
// Here:
//      vol_tag: a volume tag

//std::vector<double> getMinMaxLineCoords(int tag); // Can be incorporated in getEntityMinMaxCoords (but also is not used to my knowledge...)
// Here:
//      tag: a line tag

std::vector<double> getMinMaxSurfaceCoords(int tag); // probably no longer needed due to getEntityMinMaxCoords
// Here:
//      tag: a surface tag

std::vector<double> getEntityMinMaxCoords(const std::pair<int, int> &dimtag);
// Here:
//      dimtag: an entity's dim-tag pair






std::vector<double> cross(const std::vector<double>& a, const std::vector<double> &b);
// Here:
//      a: a typical vector in space with 3 components.
//      b: a typical vector in space with 3 components.

double norm(const std::vector<double>& a);
// Here:
//      a: a typical vector in space with 3 components.

std::vector<double> computeDirectionVector(const std::vector<double> &a, const std::vector<double> &b);
// Here:
//      a: a vector containing 3 doubles: the x, y, and z coordinates of a point.
//      b: a vector containing 3 doubles: the x, y, and z coordinates of a point.

std::vector<double> computeSurfaceNormal_3Point(const std::vector<std::vector<double>> &a);
// Here:
//      a: a vector of 3 points (i.e., a size-3 vector containing vectors of 3 doubles: the x, y, and z coordinates of the surfaces's corner points)





bool compareValues(const double a, const double b, const double tol = 1e-6);
// Here:
//      a: a double to compare with b
//      b: a double to compare with a
//      tol: the tolerance for double comparisons

bool compareVectors(const std::vector<double> &a, const std::vector<double> &b, const double tol = tol_global);
// Here:
//      a: a vector of doubles to compare with vector b
//      b: a vector of doubles to compare with vector a
//      tol: the tolerance for double comparisons

bool compareLines(const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b, const double tol = tol_global);
// Here:
//      a: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points.
//      b: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points.
//      tol: the tolerance for double comparisons

std::vector<std::vector<double>> getThreeSurfacePoints(const std::vector<std::vector<std::vector<double>>> &a);
// Here:
//      a: a vector of multiple vectors of 2 vectors containing 3 doubles: the x, y, and z coordinates of the line's end points.

bool areParallel(const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b, const double tol = tol_global);
// Here:
//      a: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points.
//      b: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points.
//      tol: the tolerance for double comparisons

bool areParallel(const std::vector<std::vector<std::vector<double>>> a, const std::vector<std::vector<std::vector<double>>> b, const double tol = tol_global);
// Here:
//      a: a vector of multiple vectors of 2 vectors containing 3 doubles: the x, y, and z coordinates of the line's end points.
//      b: a vector of multiple vectors of 2 vectors containing 3 doubles: the x, y, and z coordinates of the line's end points.
//      tol: the tolerance for double comparisons

bool areParallel(const std::vector<std::vector<std::vector<double>>> a, const std::vector<double> b, const double tol);
// Here:
//      a: a vector of multiple vectors of 2 vectors containing 3 doubles: the x, y, and z coordinates of the line's end points.
//      b: a typical vector in space with 3 components.
//      tol: the tolerance for double comparisons

bool isCoincides(const int ln_tag, const int dim, const double test_val, const double tol = tol_global);
// Here:
//      ln_tag: the line tag
//      dim: 0, 1, or 2 to check the x, y, or z dimension, respectively
//      test_val: the coordinate value to check for alignment along
//      tol: the tolerance for double comparisons

bool isCoincidesSurface(const int srfc_tag, const int dim, const double test_val, const double tol = tol_global);
// Here:
//      srfc_tag: the surface tag
//      dim: 0, 1, or 2 to check the x, y, or z dimension, respectively
//      test_val: the coordinate value to check for alignment along
//      tol: the tolerance for double comparisons

bool compareLineEndPoints(const std::vector<std::vector<double>> a, const std::vector<std::vector<double>> b, const double tol = tol_global);
// Here:
//      a: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points
//      b: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points
//      tol: the tolerance for double comparisons

bool compareSurfaceEndPoints(const std::vector<std::vector<std::vector<double>>> a, const std::vector<std::vector<std::vector<double>>> b, const double tol = tol_global);
// Here:
//      a: a vector of multiple vectors of 2 vectors containing 3 doubles: the x, y, and z coordinates of the line's end points
//      b: a vector of multiple vectors of 2 vectors containing 3 doubles: the x, y, and z coordinates of the line's end points
//      tol: the tolerance for double comparisons

bool areAdjacentEntities(const int entity_dim, const int A_tag, const int B_tag);
// Here:
//      entity_dim: the dimension of the geometric entities (i.e., 1 for lines, 2 for surfaces, 3 for volumes)
//      A_tag: the tag of entity A
//      B_tag: the tag of entity B
bool areAdjacentEntities(const std::vector<std::pair<int, int>> &A_in_dimtag, const std::vector<std::pair<int, int>> &B_in_dimtag,
    std::vector<std::pair<int, int>> &A_out_dimtag, std::vector<std::pair<int, int>> &B_out_dimtag);
// Here:
//      A_in_dimtag: a vector of dimension-tag pairs of A entities to get the boundaries for (i.e., this is the first entry in gmsh::model::getBoundary())
//      B_in_dimtag: a vector of dimension-tag pairs of B entities to get the boundaries for (i.e., this is the first entry in gmsh::model::getBoundary())
//      A_out_dimtag: a vector of dimension-tag boundary pairs of A (i.e., this is the second entry in gmsh::model::getBoundary())
//      B_out_dimtag: a vector of dimension-tag boundary pairs of B (i.e., this is the second entry in gmsh::model::getBoundary())

bool shareAnyBoundaryEntities(const int entity_dim, const int A_tag, const int B_tag);
// Here:
//      entity_dim: the dimension of the geometric entities (i.e., 1 for lines, 2 for surfaces, 3 for volumes)
//      A_tag: the tag of entity A
//      B_tag: the tag of entity B

bool linesIntersect(const std::vector<std::vector<double>>& ln1, const std::vector<std::vector<double>>& ln2, double tol = 1e-10);
// Here:
//      ln1: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points
//      ln2: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points
//      tol: the tolerance for double comparisons

bool isCollinear(const std::vector<std::vector<double>> &ln_pts, const std::vector<double> &pt, double tol = 1e-10);
// Here:
//      ln_pts: a vector of 2 vectors that contain 3 doubles: the x, y, and z coordinates of the line's end points
//      pt: a vector of the x, y, and z coordinates of a point that is being checked for collinearity with the line
//      tol: the tolerance for double comparisons





int getPhysicalGroup(int &physGroupTag, const string side, const int dim, std::vector<int> &top_tags, std::vector<int> &right_tags, std::vector<int> &bottom_tags, std::vector<int> &left_tags);
// Here:
//      physGroupTag: an integer describing the physical group tag
//      side: a string that tells which side of the domain will be in the physical group
//      dim: the dimension of the geometric entity (i.e., 1 for lines, 2 for surfaces)
//      top_tags: the group of line tags describing the top of the domain
//      right_tags: the group of line tags describing the right of the domain
//      bottom_tags: the group of line tags describing the bottom of the domain
//      left_tags: the group of line tags describing the left of the domain





int createCircle(const double x_center, const double y_center, const double z_center, const double r, const double dx_input);
// Here:
//      x_center: the x coordinate of the circle's centerpoint
//      y_center: the y coordinate of the circle's centerpoint
//      z_center: the z coordinate of the circle's centerpoint
//      r: the radius of the circle
//      dx_input: the edge length between points on the circumferance of the circle





std::vector<int> getCutSurfaceFromFile(const string &geometry_file_path, const double scale = 1.0);
// Here:
//      geometry_file_path: the file path to the text file with the points of the surfaces

std::vector<int> makeLines(const std::vector<double> &x, const std::vector<double> &y);
// Here:
//      x: the x coordinates
//      y: the y coordinates
//      lines: an empty array that will be filled with the line numbers

void load2DToolGeo(std::vector<std::pair<int, int>> &tool_sfs, std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo, const double geo_scale, const string geometry_file_path = "None");
void load2DToolGeo(std::vector<std::pair<int, int>> &tool_sfs, std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo, const double geo_scale, const double x = 0.0, const double y = 0.0, const double z = 0.0);
// Here:
//      tool_sfs: a vector of pairs, where the first number in the pair is 2 (for surface; create load3DToolGeo for 3, volumes) and the second is the surface tag
//      tool_geo: a nested vector of the coordinates of the relevant surface points. It has indices [surface, lines for those surfaces, points in those lines, coords of those points]
//      geo_scale: the scale with which to scale the loaded/imported geometry
//      x: the x-coordinate of the center of the default circle that is loaded by this function when no geometry file path is provided
//      y: the y-coordinate of the center of the default circle that is loaded by this function when no geometry file path is provided
//      z: the z-coordinate of the center of the default circle that is loaded by this function when no geometry file path is provided
//      geometry_file_path: the file path to a tool geometry file





void MakeARMesh_Uniform2DRectangular(std::vector<int> N_AR, std::vector<double> ell, std::vector<double> L,
    std::vector<std::pair<int, int>> &ARs_uc, std::vector<std::vector<std::vector<std::vector<double>>>> &ARs_uc_geo,
    std::vector<std::vector<int>> &AR_neighbors);
//void MakeARMesh_Uniform2DRectangular(const int N_AR_x, const double &ell_x, const double &L_x,
//    const int N_AR_y, const double &ell_y, const double &L_y, std::vector<std::pair<int, int>> &ARs_uc,
//    std::vector<std::vector<std::vector<std::vector<double>>>> &ARs_uc_geo);
// Here:
//      N_AR = {N_AR_x, N_AR_y}: the number of ARs in each direction
//      ell = {ell_x, ell_y}: the length of an AR in each direction
//      L = {L_x, L_y}: the total length of the domain in each direction
//      ARs_uc: an array that is filled by this function
//      ARs_uc_geo: an array of the AR points that is filled by this function
//      AR_neighbors: given the AR number, this variable outputs a list of AR numbers that neighbor the given AR (including diagonals)

void MakeARMesh_Uniform3DRectangular(const std::vector<int> &N_AR, const std::vector<double> &ell,
    const std::vector<double> &L, std::vector<std::pair<int, int>> &ARs_uc,
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &ARs_uc_geo);
// Here:
//      N_AR = {N_AR_x, N_AR_y, N_AR_z}: the number of ARs in each direction
//      ell = {ell_x, ell_y, ell_z}: the length of an AR in each direction
//      L = {L_x, L_y, L_z}: the total length of the domain in each direction
//      ARs_uc: an array that is filled by this function
//      ARs_uc_geo: an array of the AR points that is filled by this function

// USING BELOW NOW
/*
bool createMergedARGroups(const std::vector<std::pair<int, int>> &ARs,
    double min_area_threshold, double max_AR_length, double max_AR_length_ratio,
    std::vector<std::vector<int>> &AR_tags_toBeMerged, std::vector<double> &AR_pore_space_areas);
// Here:
//      ARs: a vector of surfaces that exist in the ARs after the cut. They have the form {{dim, tag}, ...}
//      min_area_threshold: the minimum allowable AR area
//      max_AR_length: the maximum allowable edge length of an AR
//      max_AR_length_ratio: the maximum allowable ratio between AR edge lengths (i.e., how long vs. how wide an AR is)
//      AR_tags_toBeMerged: a vector of the vectors of AR tags to be merged (i.e., {{1}, {2,3}, {4,6}, {5}, ...})
//      AR_pore_space_areas: a vector of the merged AR areas
*/
bool createMergedARGroups(const std::vector<std::pair<int, int>> &ARs,
    double min_vol_threshold, double max_AR_length, double max_AR_length_ratio,
    std::vector<std::vector<int>> &AR_tags_toBeMerged, std::vector<double> &AR_pore_space_vols);
// Here:
//      ARs: a vector of surfaces that exist in the ARs after the cut. They have the form {{dim, tag}, ...}
//      min_vol_threshold: the minimum allowable AR area/volume
//      max_AR_length: the maximum allowable edge length of an AR
//      max_AR_length_ratio: the maximum allowable ratio between AR edge lengths (i.e., how long vs. how wide vs. how tall an AR is)
//      AR_tags_toBeMerged: a vector of the vectors of AR tags to be merged (i.e., {{1}, {2,3}, {4,6}, {5}, ...})
//      AR_pore_space_vols: a vector of the merged AR areas/volumes





/*
// This funtion is no longer used. It's faster/more reliable to use getNonARInterfaceEntityTags.
void separateARandGeometrySurfaceLines(const std::vector<std::pair<int, int>> &ARs,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo,
    std::vector<std::vector<std::pair<int, int>>> &ARs_lns,
    std::vector<std::vector<std::pair<int, int>>> &tool_sfs_lns);
// Here:
//      ARs: a vector of surfaces that exist in the ARs after the cut. They have the form {{dim, tag}, ...}
//      tool_geo: a nested vector of the coordinates of the cut surface points. It has indices [surface, lines for those surfaces, points in those lines, coords of those points]
//      ARs_lns: a nested vector of: first level = AR, second level = AR boundary lines, and then the pair is {dim tag}... I think...
//      tool_sfs_lns: a nested vector of: first level = tool surfaces, second level = surface lines, and then the pair is {dim tag}... I think...
*/
/*
// This funtion is no longer used. It's faster/more reliable to use getNonARInterfaceEntityTags.
void separateARandGeometryVolumeSurfaces(const std::vector<std::pair<int, int>> &ARs,
    const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &tool_geo_3D,
    std::vector<std::vector<std::pair<int, int>>> &ARs_srfcs,
    std::vector<std::vector<std::pair<int, int>>> &tool_vol_srfcs);
// Here:
//      ARs: a vector of surfaces that exist in the ARs after the cut. They have the form {{dim, tag}, ...}
//      tool_geo_3D: a nested vector of the coordinates of the cut volume points. It has indices [volume, surfaces for those volumes, lines for those surfaces, points in those lines, coords of those points]
//      ARs_srfcs: a nested vector of: first level = AR, second level = AR boundary surfaces, and then the pair is {dim, tag}... I think...
//      tool_vol_srfcs: a nested vector of: first level = tool volumes, second level = volume surfaces, and then the pair is {dim, tag}... I think...
*/

// Function that removes AR interface lines/surfaces from the resulting cut geometry
void getNonARInterfaceEntityTags(const std::vector<std::pair<int, int>> &ARs,
    std::vector<std::vector<int>> &entity_tags_in_AR, int AR_dim);
// Here:
//      ARs: a vector of surfaces that exist in the ARs after the cut. They have the form {{dim, tag}, ...}
//      entity_tags_in_AR: a nested vector of: first level = AR, second level = AR interface entity tags
//      AR_dim: the dimension of the AR (i.e., 2 for surfaces, 3 for volumes)

void getPostCutTagsForLine(const std::vector<std::vector<double>> &target_ln,
    const std::vector<std::vector<std::vector<double>>> &cut_geo_ln_coords,
    const std::vector<int> &entity_tags_in_AR, const std::vector<double> &domain_boundaries,
    std::vector<int> &ln_tags, std::vector<int> &skip_ind);
// Here:
//      target_ln: a vector of the end point coordinates of the line that we want the post-cut line tags for
//      cut_geo_ln_coords: a nested vector that holds the end point coordinates of each post-cut line. Indices are [line, end points, x-y-z coord]
//      entity_tags_in_AR: a vector that contains the post-cut line tags. The order of this variable corresponds to the order of lines in cut_geo_ln_coords
//      domain_boundaries: a vector that contains the domain boundaries: {min x-coord, max x-coord, min y-coord, max y-coord}
//      ln_tags: (output) a vector of the post-cut line tags of lines that make up the target line
//      skip_ind: a vector of the post-cut line indices to skip during analysis (because such lines have already been determined as portions of the target line)

void getToolGeoTags(const std::vector<std::vector<int>> &cut_inds,
    const std::vector<std::vector<int>> &entity_tags_in_AR,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo,
    std::vector<double> domain_boundaries, std::vector<std::vector<int>> &saved_tags);
// Here:
//      cut_inds: Sets of pre-cut surface indices (i.e., indices of tool_geo) that are desired to be defined with a unique physical group. For example, {{2}, {3,5}} says define one
//          physical group for surface tool_geo[2], and one physical group that contains surfaces tool_geo[3] and tool_geo[5]
//      entity_tags_in_AR: a nested vector that contains the post-cut line tags for each AR
//      tool_geo: a nested vector of the coordinates of the relevant surface points. It has indices [surface, lines for those surfaces, points in those lines, coords of those points]
//      domain_boundaries: a vector of left, right, bottom, and top boundaries of the rectangular domain (assumes a rectangular domain). The notation is: {left-most x coord, right-most x coord, bottom-most y coord, top-most y coord}
//      saved_tags: (output) a nested vector of post-cut line tags that make up the surfaces specified in cut_inds for each physical group

void getDomainBoundaryLineTags(const std::vector<std::vector<int>> &ln_tags_in_AR,
    std::vector<std::vector<int>> &domain_boundary_tags, const std::vector<std::pair<int, int>> &ARs,
    std::vector<std::vector<int>> &AR_tags_neighbors, double L_x, double L_y, double tol = tol_global);
// Here:
//      ln_tags_in_AR: a nested vector of: first level = AR, second level = AR interface line tags
//      domain_boundary_tags: a nested vector of: first level = top, bottom, left, right, and unused line tags, second level = line tags
//      ARs: a vector gmsh formatted AR surface tags. They have the form {{dim, tag}, ...}
//      AR_tags_neighbors: a nested vector of AR inds that touch the top, bottom, left, and right domain boundaries. First level = top, bottom, left, and right; second level = AR inds
//      L_x: the length of the domain
//      L_y: the height of the domain

void getDomainBoundarySurfaceTags(const std::vector<std::vector<int>> &surf_tags_in_AR,
    std::vector<std::vector<int>> &domain_boundary_tags, const std::vector<double> &L);
// Here:
//      surf_tags_in_AR: a nested vector of: first level = AR, second level = AR interface surface tags
//      domain_boundary_tags: a nested vector of: first level = top, bottom, left, right, front, back, and unused surface tags, second level = surface tags
//      L = {L_x, L_y, L_z}: the length of the domain in each direction
