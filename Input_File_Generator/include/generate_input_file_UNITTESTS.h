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
#include "generate_input_file_TOOLS.h"



using namespace std;



// ===============================================================
//   Define a function for generating a default/tutorial config
//   file for the MoFA model formulation process.
// ===============================================================

void createUnitTestUpscaledConfig(string &project_dir, string &config_file_name, const int unit_test_number)
{
    // ===============================================================
    //   Define the config file dictionary.
    // ===============================================================
    JSONDict configFile;
    
    
    // ===============================================================
    //   Define the mesh dictionary.
    // ===============================================================
    JSONDict meshDict;
    
    
    if (unit_test_number == 0)
    {
        // Unit-test specs: 2D-to-1D, heterogeneous media (tutorial porous media geometry), 10x1 AR, traditional MoFA, Pe = \epsilon^{-1}
        int is3D = 0;
        vector<int> periodicity = {1, 0, 0};
        istringstream meshDict_stream = createMeshDict(project_dir, "upscaled_mesh", "upscaled",
            1, is3D, periodicity, {0.02, 0.02},
            0.1, {10, 1}, {10.0, 1.0}, {1.0, 1.0},
            0, 0.2, 1.75, 2.0,
            project_dir + "data/", "Tutorial_Geometry.txt", 1.0,
            {}, {});
        meshDict.loadFromStream( meshDict_stream );
        configFile["mesh"] = &meshDict;
    }
    else if (unit_test_number == 1)
    {
        // Unit-test specs: 2D-to-2D, heterogeneous media (.tif geometry), 10x10 AR with merge, localized mesh, generalized residual, Pe = 
        int is3D = 0;
        vector<int> periodicity = {0, 0, 0};
        istringstream meshDict_stream = createMeshDict(project_dir, "upscaled_mesh", "upscaled",
            1, is3D, periodicity, {0.05, 0.05},
            0.1, {10, 10, 0}, {10.0, 10.0, 0.0}, {1.0, 1.0, 0.0},
            1, 0.2, 1.75, 2.0,
            project_dir + "data/", "Sintered_Geometry.txt", 10.0,
            {}, {});
        meshDict.loadFromStream( meshDict_stream );
        configFile["mesh"] = &meshDict;
    }
    else
    {
        cerr << "generate_input_file_UNITTESTS.h: createUnitTestUpscaledConfig(): CRITICAL ERROR: Provided unit test ID unrecognized. Unit test ID: " << unit_test_number << "." << endl;
        exit(1);
    }



    // ===============================================================
    //   Define the stokes solver dictionary.
    // ===============================================================
    JSONDict stokesDict;
    

    if (unit_test_number == 0)
    {
        // Unit-test specs: 2D-to-1D, heterogeneous media (tutorial porous media geometry), 10x1 AR, traditional MoFA, Pe = \epsilon^{-1}
        vector<int> periodicity = {1, 0, 0};
        istringstream stokesDict_stream = createStokesSolverDict(project_dir,
            "upscaled_u_sol.gf", "upscaled_p_sol.gf", "upscaled_mesh.mesh",
            2, 1, periodicity, 1.0,
            0, {337.84, 0.0, 0.0}, 0.0,
            10000, 1.0e-10, 0.0);
        stokesDict.loadFromStream( stokesDict_stream );
        configFile["stokes"] = &stokesDict;
    }
    else if (unit_test_number == 1)
    {
        // Unit-test specs: 2D-to-2D, heterogeneous media (.tif geometry), 10x10 AR with merge, localized mesh, generalized residual, Pe = 
        vector<int> periodicity = {0, 0, 0};
        istringstream stokesDict_stream = createStokesSolverDict(project_dir,
            "upscaled_u_sol.gf", "upscaled_p_sol.gf", "upscaled_mesh.mesh",
            2, 1, periodicity, 1.0,
            1, {0.0, 0.0, 0.0}, 0.2320,
            10000, 1.0e-10, 0.0);
        stokesDict.loadFromStream( stokesDict_stream );
        configFile["stokes"] = &stokesDict;
    }
    else
    {
        cerr << "generate_input_file_UNITTESTS.h: createUnitTestUpscaledConfig(): CRITICAL ERROR: Provided unit test ID unrecognized. Unit test ID: " << unit_test_number << "." << endl;
        exit(1);
    }
    


    // ===============================================================
    //   Define the closure problem solver dictionary.
    // ===============================================================
    JSONDict closureDict;


    if (unit_test_number == 0)
    {
        // Unit-test specs: 2D-to-1D, heterogeneous media (tutorial porous media geometry), 10x1 AR, traditional MoFA, Pe = \epsilon^{-1}
        vector<int> periodicity = {0, 0, 0};
        istringstream transportClosureDict_stream = createTransportClosureDict(project_dir,
            2, periodicity, 1, 1, 0,
            1.0, 1.0, {0.0},
            0, 1, 0,
            0.0, 1.0, 0.0,
            10000, 1.0e-13, 1.0e-8);
        closureDict.loadFromStream( transportClosureDict_stream );
        configFile["scalar closure"] = &closureDict;
    }
    else if (unit_test_number == 1)
    {
        // Unit-test specs: 2D-to-2D, heterogeneous media (.tif geometry), 10x10 AR with merge, localized mesh, generalized residual, Pe = 
        vector<int> periodicity = {0, 0, 0};
        istringstream transportClosureDict_stream = createTransportClosureDict(project_dir,
            2, periodicity, 1, 1, 0,
            1.0, 1.0, {0.0},
            1, 1, 0,
            0.0, 1.0, 0.15,
            10000, 1.0e-13, 1.0e-8);
        closureDict.loadFromStream( transportClosureDict_stream );
        configFile["scalar closure"] = &closureDict;
    }
    else
    {
        cerr << "generate_input_file_UNITTESTS.h: createUnitTestUpscaledConfig(): CRITICAL ERROR: Provided unit test ID unrecognized. Unit test ID: " << unit_test_number << "." << endl;
        exit(1);
    }



    // ===============================================================
    //   Define the upscaled solver dictionary.
    // ===============================================================
    JSONDict upscaledDict;

    
    if (unit_test_number == 0)
    {
        // Unit-test specs: 2D-to-1D, heterogeneous media (tutorial porous media geometry), 10x1 AR, traditional MoFA, Pe = \epsilon^{-1}
        istringstream transportUpscaledDict_stream = createTransportUpscaledDict(project_dir,
            1.0,
            1000, 0.0001, 100);
        upscaledDict.loadFromStream( transportUpscaledDict_stream );
        configFile["upscaled"] = &upscaledDict;
    }
    else if (unit_test_number == 1)
    {
        // Unit-test specs: 2D-to-2D, heterogeneous media (.tif geometry), 10x10 AR with merge, localized mesh, generalized residual, Pe = 
        istringstream transportUpscaledDict_stream = createTransportUpscaledDict(project_dir,
            1.0,
            10000, 0.00001, 100);
        upscaledDict.loadFromStream( transportUpscaledDict_stream );
        configFile["upscaled"] = &upscaledDict;
    }
    else
    {
        cerr << "generate_input_file_UNITTESTS.h: createUnitTestUpscaledConfig(): CRITICAL ERROR: Provided unit test ID unrecognized. Unit test ID: " << unit_test_number << "." << endl;
        exit(1);
    }



    // ===============================================================
    //   Define the error calculator dictionary.
    // ===============================================================
    JSONDict errorDict;
    istringstream dict = createErrorDict(project_dir, "transport");
    errorDict.loadFromStream( dict );
    configFile["error calc"] = &errorDict;



    // ===============================================================
    //   Save the config file.
    // ===============================================================
    configFile.saveToFile(project_dir + "config/" + config_file_name);
}










// ===============================================================
//   Define a function for generating a default/tutorial config
//   file for the porescale model formulation process.
// ===============================================================

void createUnitTestPorescaleConfig(string &project_dir, string &config_file_name, const int unit_test_number)
{
    
    // ===============================================================
    //   Define the config file dictionary.
    // ===============================================================
    JSONDict configFile;
    

    
    // ===============================================================
    //   Define the mesh dictionary.
    // ===============================================================
    JSONDict meshDict;


    if (unit_test_number == 0)
    {
        // For running the tutorial unit-test (typical 2D, heterogeneous media, 10x1 AR, traditional MoFA, Pe = \epsilon^{-1})
        int is3D = 0;
        vector<int> periodicity = {1, 0, 0};
        istringstream meshDict_stream = createMeshDict(project_dir, "porescale_mesh", "porescale",
            1, is3D, periodicity, {0.002, 0.002},
            0.1, {10, 1}, {1.0, 0.1}, {0.1, 0.1},
            0, 0.002, 0.175, 2.0,
            project_dir + "data/", "Tutorial_Geometry.txt", 0.1,
            {}, {});
        meshDict.loadFromStream( meshDict_stream );
        configFile["mesh"] = &meshDict;
    }
    else if (unit_test_number == 1)
    {
        // Unit-test specs: 2D-to-2D, heterogeneous media (.tif geometry), 10x10 AR with merge, localized mesh, generalized residual, Pe = 
        int is3D = 0;
        vector<int> periodicity = {0, 0, 0};
        istringstream meshDict_stream = createMeshDict(project_dir, "porescale_mesh", "porescale",
            1, is3D, periodicity, {0.005, 0.005},
            0.1, {10, 10, 0}, {1.0, 1.0, 0.0}, {0.1, 0.1, 0.0},
            1, 0.002, 0.175, 2.0,
            project_dir + "data/", "Sintered_Geometry.txt", 1.0,
            {}, {});
        meshDict.loadFromStream( meshDict_stream );
        configFile["mesh"] = &meshDict;
    }
    else
    {
        cerr << "generate_input_file_UNITTESTS.h: createUnitTestPorescaleConfig(): CRITICAL ERROR: Provided unit test ID unrecognized. Unit test ID: " << unit_test_number << "." << endl;
        exit(1);
    }



    // ===============================================================
    //   Define the stokes solver dictionary.
    // ===============================================================
    JSONDict stokesDict;

    
    if (unit_test_number == 0)
    {
        // Unit-test specs: 2D-to-1D, heterogeneous media (tutorial porous media geometry), 10x1 AR, traditional MoFA, Pe = \epsilon^{-1}
        vector<int> periodicity = {1, 0, 0};
        istringstream stokesDict_stream = createStokesSolverDict(project_dir,
            "porescale_u_sol.gf", "porescale_p_sol.gf", "porescale_mesh.mesh",
            2, 1, periodicity, 0.01,
            0, {337.84, 0.0, 0.0}, 0.0,
            10000, 1.0e-10, 0.0);
        stokesDict.loadFromStream( stokesDict_stream );
        configFile["stokes"] = &stokesDict;
    }
    else if (unit_test_number == 1)
    {
        // Unit-test specs: 2D-to-2D, heterogeneous media (.tif geometry), 10x10 AR with merge, localized mesh, generalized residual, Pe = 
        vector<int> periodicity = {0, 0, 0};
        istringstream stokesDict_stream = createStokesSolverDict(project_dir,
            "porescale_u_sol.gf", "porescale_p_sol.gf", "porescale_mesh.mesh",
            2, 1, periodicity, 0.01,
            1, {0.0, 0.0, 0.0}, 0.2320,
            10000, 1.0e-10, 0.0);
        stokesDict.loadFromStream( stokesDict_stream );
        configFile["stokes"] = &stokesDict;
    }
    else
    {
        cerr << "generate_input_file_UNITTESTS.h: createUnitTestPorescaleConfig(): CRITICAL ERROR: Provided unit test ID unrecognized. Unit test ID: " << unit_test_number << "." << endl;
        exit(1);
    }



    // ===============================================================
    //   Define the porescale transport solver dictionary.
    // ===============================================================
    JSONDict porescaleDict;


    if (unit_test_number == 0)
    {
        // Unit-test specs: 2D-to-1D, heterogeneous media (tutorial porous media geometry), 10x1 AR, traditional MoFA, Pe = \epsilon^{-1}
        vector<int> periodicity = {0, 0, 0};
        istringstream transportPorescaleDict_stream = createTransportPorescaleDict(project_dir,
            1, periodicity, 1, 1, 0,
            10.0, 1.0, {0.0}, 1.0,
            1000, 0.0001, 100);
        porescaleDict.loadFromStream( transportPorescaleDict_stream );
        configFile["scalar porescale"] = &porescaleDict;
    }
    else if (unit_test_number == 1)
    {
        // Unit-test specs: 2D-to-2D, heterogeneous media (.tif geometry), 10x10 AR with merge, localized mesh, generalized residual, Pe = 
        vector<int> periodicity = {0, 0, 0};
        istringstream transportPorescaleDict_stream = createTransportPorescaleDict(project_dir,
            1, periodicity, 1, 1, 0,
            10.0, 1.0, {0.0}, 1.0,
            10000, 0.00001, 100);
        porescaleDict.loadFromStream( transportPorescaleDict_stream );
        configFile["scalar porescale"] = &porescaleDict;
    }
    else
    {
        cerr << "generate_input_file_UNITTESTS.h: createUnitTestPorescaleConfig(): CRITICAL ERROR: Provided unit test ID unrecognized. Unit test ID: " << unit_test_number << "." << endl;
        exit(1);
    }

    

    // ===============================================================
    //   Save the config file.
    // ===============================================================
    configFile.saveToFile(project_dir + "config/" + config_file_name);
}

