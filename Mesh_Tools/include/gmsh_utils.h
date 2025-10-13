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

std::vector<double> getMinMaxLineCoords(int tag);
// Here:
//      tag: a line tag

std::vector<double> getMinMaxSurfaceCoords(int tag);
// Here:
//      tag: a surface tag





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

bool compareVectors(const std::vector<double> a, const std::vector<double> b, const double tol = tol_global);
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





void MakeARMesh_Uniform2DRectangular(const int N_AR_x, const double &ell_x, const double &L_x,
    const int N_AR_y, const double &ell_y, const double &L_y, std::vector<std::pair<int, int>> &ARs_uc,
    std::vector<std::vector<std::vector<std::vector<double>>>> &ARs_uc_geo);
// Here:
//      N_AR_x: the number of ARs in the x-direction
//      ell_x: the length of an AR in the x-direction
//      L_x: the total length of the domain in the x-direction
//      N_AR_y: the number of ARs in the y-direction
//      ell_y: the length of an AR in the y-direction
//      L_y: the total length of the domain in the y-direction
//      ARs_uc: an array that is filled by this function
//      ARs_uc_geo: an array of the AR points that is filled by this function

void MakeARMesh_Uniform3DRectangular(const int N_AR_x, const double &ell_x, const double &L_x,
    const int N_AR_y, const double &ell_y, const double &L_y,
    const int N_AR_z, const double &ell_z, const double &L_z,
    std::vector<std::pair<int, int>> &ARs_uc,
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &ARs_uc_geo);
// Here:
//      N_AR_x: the number of ARs in the x-direction
//      ell_x: the length of an AR in the x-direction
//      L_x: the total length of the domain in the x-direction
//      N_AR_y: the number of ARs in the y-direction
//      ell_y: the length of an AR in the y-direction
//      L_y: the total length of the domain in the y-direction
//      N_AR_z: the number of ARs in the z-direction
//      ell_z: the length of an AR in the z-direction
//      L_z: the total length of the domain in the z-direction
//      ARs_uc: an array that is filled by this function
//      ARs_uc_geo: an array of the AR points that is filled by this function

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





void separateARandGeometrySurfaceLines(const std::vector<std::pair<int, int>> &ARs,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &tool_geo,
    std::vector<std::vector<std::pair<int, int>>> &ARs_lns,
    std::vector<std::vector<std::pair<int, int>>> &tool_sfs_lns);
// Here:
//      ARs: a vector of surfaces that exist in the ARs after the cut. They have the form {{dim, tag}, ...}
//      tool_geo: a nested vector of the coordinates of the cut surface points. It has indices [surface, lines for those surfaces, points in those lines, coords of those points]
//      ARs_lns: a nested vector of: first level = AR, second level = AR boundary lines, and then the pair is {dim tag}... I think...
//      tool_sfs_lns: a nested vector of: first level = tool surfaces, second level = surface lines, and then the pair is {dim tag}... I think...

void separateARandGeometryVolumeSurfaces(const std::vector<std::pair<int, int>> &ARs,
    const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &tool_geo_3D,
    std::vector<std::vector<std::pair<int, int>>> &ARs_srfcs,
    std::vector<std::vector<std::pair<int, int>>> &tool_vol_srfcs);
// Here:
//      ARs: a vector of surfaces that exist in the ARs after the cut. They have the form {{dim, tag}, ...}
//      tool_geo_3D: a nested vector of the coordinates of the cut volume points. It has indices [volume, surfaces for those volumes, lines for those surfaces, points in those lines, coords of those points]
//      ARs_srfcs: a nested vector of: first level = AR, second level = AR boundary surfaces, and then the pair is {dim, tag}... I think...
//      tool_vol_srfcs: a nested vector of: first level = tool volumes, second level = volume surfaces, and then the pair is {dim, tag}... I think...
