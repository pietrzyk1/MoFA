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
#include "JSON_IO.h"



using namespace std;



istringstream createErrorDict(const string &project_dir, string sim_type)
{
    JSONDict errorDict;

    JSONDict errorDict_output;
    errorDict_output["directory"] = project_dir + "output/absolute_error/"; // The absolute error output directory
    errorDict_output["file name"] = "abs_error.txt"; // The absolute error output file name
    errorDict["output path"] = &errorDict_output;

    JSONDict errorDict_avgpore;
    if (sim_type == "transport") {
        errorDict_avgpore["directory"] = project_dir + "output/avg_porescale_solution/"; // The directory where the average porescale solution is stored
    }
    else {
        errorDict_avgpore["directory"] = project_dir + "output/avg_unsteady_stokes_solution/"; // The directory where the average porescale solution is stored
    }
    errorDict_avgpore["file name"] = "avg_sol.txt"; // The average porescale solution file name
    errorDict["avg porescale path"] = &errorDict_avgpore;

    JSONDict errorDict_upscaled;
    if (sim_type == "transport") {
        errorDict_upscaled["directory"] = project_dir + "output/upscaled_solution/"; // The directory where the upscaled solution is stored
    }
    else {
        errorDict_upscaled["directory"] = project_dir + "output/upscaled_unsteady_stokes_solution/"; // The directory where the upscaled solution is stored
    }
    errorDict_upscaled["file name"] = "upscaled_sol.txt"; // The upscaled solution file name
    errorDict["upscaled path"] = &errorDict_upscaled;

    return errorDict.saveToStream();
}

