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

#include <iostream>
#include <fstream>
#include "JSON_IO.h"



using namespace std;



// ===============================================================
// Define save/load functions for mesh attributes
// ===============================================================
// Define a function to save the custom type to a text file
void saveToFile(const CustomMap &data, const string &filename)
{
    // Create the output file
    ofstream file(filename);

    // Verify that the file was created
    if (!file)
    {
        cerr << "CRITICAL ERROR: JSON_IO.cpp: saveToFile(CustomMap): Error opening file for writing." << endl;
        return;
    }

    for (const auto &outer_pair : data)
    {
        file << outer_pair.first << ":\n";  // Write the outer key
        for (const auto &inner_pair : outer_pair.second)
        {
            file << "  " << inner_pair.first << ": " << inner_pair.second << "\n";  // Write inner key-value pairs
        }
        file << "\n";  // Blank line for separation between entries
    }

    file.close();
    cout << "Boundary attribute data has been saved to " << filename << endl;
}

// Define a function to load the custom type from a text file
CustomMap loadFromFile(const string &filename)
{
    // Find the previously generated output file
    ifstream file(filename);
    CustomMap data;
    
    if (!file)
    {
        cerr << "CRITICAL ERROR: JSON_IO.cpp: loadFromFile(CustomMap): Error opening file for reading." << endl;
        return data;
    }

    string line;
    string currentKey;
    
    while (getline(file, line))
    {
        // Ignore empty lines
        if (line.empty()) continue;

        // If the line ends with a colon, it's an outer key
        if (line.back() == ':')
        {
            currentKey = line.substr(0, line.size() - 1); // Remove the colon
        }
        else // Otherwise, it's a key-value pair (with an indentation of two spaces)
        {
            size_t pos = line.find(": ");
            if (pos != string::npos)
            {
                string innerKey = line.substr(2, pos - 2); // Extract inner key (skip "  ")
                string value_str = line.substr(pos + 2);  // Extract the value as a string
                int value_int = stoi(value_str);  // Turn the string into an int
                data[currentKey][innerKey] = value_int; // Record the value
            }
        }
    }

    file.close();
    return data;
}
