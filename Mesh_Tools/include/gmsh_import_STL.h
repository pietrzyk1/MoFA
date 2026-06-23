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
#include <vector>
#include <string>
#include <gmsh.h>





// ========================================
// Support functions for importing STL geometry into gmsh
// ========================================

std::vector<int> assignARIndex(const std::vector<double> &coord, const std::vector<std::vector<double>> &AR_planes, const double tol);

template <typename T, typename lamfunc>
void removeDuplicates(std::vector<T> &vec, lamfunc func);

template <typename T, typename lamfunc>
void getUniqueEntries(std::vector<T> &vec, lamfunc func);





// ========================================
// Functions for finding curve loops in a vector of gmsh line tags
// ========================================

std::vector<std::vector<int>> findCurveLoops(std::vector<int> &glns, double tol, bool verbose = false);

void findCurveLoopsRecursive(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns,
    const std::vector<std::vector<std::pair<int, int>>> &bndry_pts, std::vector<std::vector<int>> &curveLoops, const double tol, bool verbose = false);

void obtainFirstLineForCurveLoopFinder(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns,
    const std::vector<std::vector<std::pair<int, int>>> &bndry_pts, std::vector<int> &curveLoop, bool verbose = false);

void followCurveLoop(int &target_pt, std::vector<int> &unused_inds, const std::vector<int> &glns, const std::vector<std::vector<std::pair<int, int>>> &bndry_pts,
    std::vector<int> &curveLoop, const string &add_mode, bool verbose = false);

int getCurveLoopEndPointCoords(const std::vector<int> &curveLoop, std::vector<std::vector<double>> &end_coords);





// ========================================
// Functions for finding surface loops in a vector of gmsh surface tags
// ========================================

std::vector<std::vector<int>> findSurfaceLoops(const std::vector<int> &gsurfs, bool verbose = false);

void findSurfaceLoopsRecursive(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<std::vector<int>> &surfaceLoops, bool verbose = false);

void obtainFirstSurfaceForSurfaceLoopFinder(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<int> &surfaceLoop, bool verbose = false);

void followSurfaceLoop(std::vector<int> &target_lns, std::vector<int> &unused_inds, const std::vector<int> &gsurfs,
    const std::vector<std::vector<int>> &bndry_glns, std::vector<int> &surfaceLoop, bool verbose = false);





// ========================================
// Functions for finding surface loops in a vector of gmsh surface tags
// ========================================

std::vector<std::vector<std::vector<int>>> arrangeCurveLoopsForSurfaces(const std::vector<std::vector<int>> &curveLoops, int dim, const double tol, bool verbose = false);





// ========================================
// Function for getting the triangle points and normals from an STL file
// ========================================

void getTriPointsFromSTL(const string &STL_file_path, const double &geo_scale,
    std::vector<std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>> &coord_ln_surf_pairs,
    std::vector<std::vector<double>> &surf_normals, int &tri_count);
// Here:
//      STL_file_path: the file path to an STL file
//      geo_scale: the scale with which to scale the loaded/imported geometry
//      coord_ln_surf_pairs: a tuple of information for each point: coordinates (vector<double>), indices of the lines that the point is a part of (vector<int>), and indices of the surfaces that the points belong to (the second vector<int>)
//      surf_normals: a vector of the normal vectors for each STL surface
//      tri_count: the number of triangles/surfaces defined in the STL

void getTrianglesFromSTL(const string &STL_file_path, const double &geo_scale,
    std::vector<std::vector<std::vector<double>>> &tri_coords,
    std::vector<std::vector<double>> &tri_normals, int &tri_count, std::vector<double> &min_max_x_y_z);
// Here:
//      STL_file_path: the file path to an STL file
//      geo_scale: the scale with which to scale the loaded/imported geometry
//      tri_coords: a vector of triangles described by another vector of their point coordinates
//      tri_normals: a vector of the normal vectors for each STL triangle surface
//      tri_count: the number of triangles/surfaces defined in the STL
//      min_max_x_y_z: the minimum and maximum x y and z coordinate values

void getCoordFromSTLLine(const string &line, std::vector<double> &coord, const double &scale, const int pos = 0);
// Here:
//      line: a line (string) from an STL file
//      coord: an empty vector of doubles to put the coordinate in
//      scale: the scale with which to scale the loaded/imported geometry
//      pos: an initial position from which to read the line





// ========================================
// Function for importing STL cut geometry into gmsh
// ========================================

int createCutVolumeFromSTL(const string &STL_file_path, const double &geo_scale);
// Here:
//      STL_file_path: the file path to an STL file
//      geo_scale: the scale with which to scale the loaded/imported geometry





// ========================================
// Function for importing STL geometry in gmsh (primarily for homogenization closure problem simulations)
// ========================================

std::vector<std::pair<int, int>> splitRectangularSTL(const string &STL_file_dir, const string &STL_file_name, const double &geo_scale);
// Here:
//      STL_file_dir: the file directory to the STL file for mesh creation
//      STL_file_name: the STL file name for mesh creation
//      geo_scale: the scale with which to scale the loaded/imported geometry

void writeSTLFile(const string &STL_file_path, const std::vector<std::vector<std::vector<double>>> tri_coords,
    std::vector<std::vector<double>> tri_normals);
// Here:
//      STL_file_path: the file path to the STL file to be written
//      tri_coords: a nested vector with the structure: { triangles { points {x,y,z coords} } }
//      tri_normals: a nested vector with the normal vectors of each triangle


void changeElementTags(const string msh_file_path, const string new_msh_file_path,
    const std::vector<std::vector<size_t>> &elemTags_split, const std::vector<int> &new_pgs,
    const std::vector<int> &new_elem_entities);
// Here:
//      msh_file_path: the file path to the STL mesh file to be editted
//      new_msh_file_path: the file path to the new/editted STL mesh file
//      elemTags_split: a nested vector of the element tags grouped by those belonging to the same physical group
//      new_pgs: the new physical group tags
//      new_elem_entities: the new element entity tags

std::vector<std::pair<int, int>> importSTLVolume(const string &STL_file_path, const double &geo_scale);
// Here:
//      STL_file_path: the file path to an STL file
//      geo_scale: the scale with which to scale the loaded/imported geometry





// ========================================
// Function for importing and cutting STL geometry in gmsh
// ========================================

std::vector<std::pair<int, int>> importAndCutSTLVolume(const string &STL_file_path, const double &geo_scale, const std::vector<std::vector<double>> &AR_planes);
// Here:
//      STL_file_path: the file path to an STL file
//      geo_scale: the scale with which to scale the loaded/imported geometry
//      AR_planes: lists of x, y, and z coordinates that describe planes that cut the geometry into AR regions (e.g., {{1, 2}, {3}, {1}} cuts the geometry with planes x = 1, x = 2, y = 3, and z = 1)
