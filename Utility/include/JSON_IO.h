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
#include <fstream>
#include <memory>
#include <unordered_map>
#include <map>
#include <tuple>
#include <vector>
#include <any>
#include <variant>

#include <sstream>
#include <type_traits>
#include <typeinfo>




using namespace std;




// ===============================================================
// Declare custom type and save/load functions for mesh attributes
// ===============================================================
// Define an alias for the custom type
using CustomMap = unordered_map<string, unordered_map<string, int>>;

// Declare save function
void saveToFile(const CustomMap &data, const string &filename);

// Declare load function
CustomMap loadFromFile(const string &filename);









// Forward declaration for recursive type
struct JSONDict;
struct JSONVal;


using JSONARGType = variant< int, double, string, vector<JSONVal>, vector<vector<JSONVal>>, JSONDict* >;


// Define a basic structure for the JSON dictionary arguments (see JSONARGType for what types can be handled)
struct JSONVal : JSONARGType
{
    // Inherit constructors
    using JSONARGType::JSONARGType;

    // Allow assignment from the same base type (This has not been an issue yet...)
    //Value2& operator=(const Value2&) = default;
    //Value2& operator=(Value2&&) = default;
    
    // These are for writing JSONVal k = ...;
    //
    // Allows JSONVal variables to be assigned as a list of ints (e.g., JSONVal a = b, where b is a vector<int>)
    JSONVal(vector<int> &list)
    {
        vector<JSONVal> vec;
        for (int item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as a const list of ints (e.g., JSONVal a = b, where b is a const vector<int>)
    JSONVal(const vector<int> &list)
    {
        vector<JSONVal> vec;
        for (int item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as vector<vector<int>> (e.g., JSONVal a = b, where b is a vector<vector<int>>)
    JSONVal(vector<vector<int>> &list)
    {
        vector<vector<JSONVal>> vec;
        int count = 0;
        for (vector<int> item : list)
        {
            vec.emplace_back( vector<JSONVal>() );
            for (int item2 : item) { vec[count].emplace_back(item2); }
            count += 1;
        }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as a list of doubles (e.g., JSONVal a = b, where b is a vector<double>)
    JSONVal(vector<double> &list)
    {
        vector<JSONVal> vec;
        for (double item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as a const list of doubles (e.g., JSONVal a = b, where b is a const vector<double>)
    JSONVal(const vector<double> &list)
    {
        vector<JSONVal> vec;
        for (double item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as a list of strings (e.g., JSONVal a = b, where b is a vector<string>)
    JSONVal(vector<string> &list)
    {
        vector<JSONVal> vec;
        for (string item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    
    // Allows JSONVal variables to be assigned a list of any basic type (e.g., int, double, string) in the variant JSONARGType (e.g., JSONVal a = {1,2,3})
    template <typename T, typename = enable_if_t< is_convertible_v<T, JSONVal> >>
    JSONVal(initializer_list<T> list)
    {
        vector<JSONVal> vec;
        for (const auto& item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as a list of ints (needed because lists of ints and doubles are indistinguishable)
    JSONVal(initializer_list<int> list)
    {
        vector<JSONVal> vec;
        for (int item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned vector<vector<int>> (needed because lists of ints and doubles are indistinguishable) (e.g., JSONVal a = {{1,2},{3,4}})
    JSONVal(initializer_list<vector<int>> list)
    {
        vector<vector<JSONVal>> vec;
        int count = 0;
        for (vector<int> item : list)
        {
            vec.emplace_back( vector<JSONVal>() );
            for (int item2 : item) { vec[count].emplace_back(item2); }
            count += 1;
        }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as a list of doubles (needed because lists of ints and doubles are indistinguishable)
    JSONVal(initializer_list<double> list)
    {
        vector<JSONVal> vec;
        for (double item : list) { vec.emplace_back(item); }
        *this = vec;
    }
    // Allows JSONVal variables to be assigned as a list of mixed basic types (e.g., int, double, string) in the variant JSONARGType
    // (This allows for "python list"-like behavior; lists with heterogeneous types can be handled)
    template <typename... Args>
    JSONVal(Args&&... args)
    {
        vector<JSONVal> vec;
        (vec.emplace_back(JSONVal(forward<Args>(args))), ...);
        *this = std::move(vec);
    }

    // These are for writing int m = JSONDict k["somehting"];
    //
    explicit operator int() const {
        if (auto ptr = get_if<int>(this)) return *ptr; // If JSONVal holds an int, return it
        //if (auto ptr = get_if<double>(this)) return static_cast<int>(*ptr); // If JSONVal holds a double, return it as an int
        throw bad_variant_access();  // Error out
    }
    explicit operator double() const {
        if (auto ptr = get_if<double>(this)) return *ptr; // If JSONVal holds a double, return it
        //if (auto ptr = get_if<int>(this)) return static_cast<double>(*ptr); // If JSONVal holds an int, return it as a double
        throw bad_variant_access();  // Error out
    }
    operator string() const {
        if (auto ptr = get_if<string>(this)) return *ptr; // If JSONVal holds a string, return it
        throw bad_variant_access();  // Error out
    }
    operator JSONDict*() const {
        if (auto ptr = get_if<JSONDict*>(this)) return *ptr;
        throw bad_variant_access();  // Error out
    }
    /*
    operator vector<JSONVal>() const {
        if (auto ptr = get_if<vector<JSONVal>>(this)) return *ptr;
        throw bad_variant_access();  // Error out
    }
    */
    template <typename T, typename = enable_if_t< is_convertible_v<T, JSONVal> >>
    operator vector<T>() const {
        if (auto ptr = get_if<vector<JSONVal>>(this))
        {
            vector<T> ptr2;
            for (int i = 0; i < ptr->size(); i++) { ptr2.push_back( (T)(*ptr)[i] ); }
            return ptr2;
        }
        throw bad_variant_access();  // Error out
    }
    // For writting stuff like vector<vector<int>> var = vector<vector<JSONVal>> var2
    template <typename T, typename = enable_if_t< is_convertible_v<T, JSONVal> >>
    operator vector<vector<T>>() const {
        if (auto ptr = get_if<vector<vector<JSONVal>>>(this))
        {
            vector<vector<T>> ptr2;
            for (int i = 0; i < ptr->size(); i++)
            {
                ptr2.push_back( vector<T>() );
                for (int j = 0; j < (*ptr)[i].size(); j++) { ptr2[i].push_back( (T)(*ptr)[i][j] ); }
            }
            return ptr2;
        }
        throw bad_variant_access();  // Error out
    }



    // Function for converting JSONVal variables to a string
    string to_string(string indent = "    ", int level = 0) const
    {
        return visit([indent, level](const auto& val) -> string
        {
            using T = decay_t<decltype(val)>;
                
            if constexpr (is_same_v<T, string>)
            {
                return (string)"\"" + val + (string)"\"";
            }
            else if constexpr (is_same_v<T, JSONDict*>)
            {
                return val ? val->to_string(indent, level) : "<<ERROR IN TURNING JSONDICT TO STRING>>";
            }
            else if constexpr (is_same_v<T, vector<JSONVal>>)
            {
                ostringstream oss;
                oss << "[";
                for (int i = 0; i < val.size(); i++)
                {
                    oss << val[i].to_string(indent, level);
                    if (i + 1 < val.size()) { oss << ", "; }
                }
                oss << "]";
                return oss.str();
            }
            else if constexpr (is_same_v<T, vector<vector<JSONVal>>>)
            {
                ostringstream oss;
                oss << "[";
                for (int i = 0; i < val.size(); i++)
                {
                    oss << "[";
                    for (int j = 0; j < val[i].size(); j++)
                    {
                        oss << val[i][j].to_string(indent, level);
                        if (j + 1 < val[i].size()) { oss << ", "; }
                    }
                    oss << "]";
                    if (i + 1 < val.size()) { oss << ", "; }
                }
                oss << "]";
                return oss.str();
            }
            else if constexpr (is_same_v<T, int>)
            {
                ostringstream oss;
                oss << val;
                return oss.str();
            }
            else if constexpr (is_same_v<T, double>)
            {
                if ( (val < 1.0e+3 && val > 0.01) || (val > -1.0e+3 && val < -0.01) || (val == 0.0) )
                {
                    string val_str = std::to_string(val);
                    if (val_str.find('.') == string::npos) { val_str += '.'; }
                    ostringstream oss;
                    oss << val_str;
                    return oss.str();
                }
                else
                {
                    ostringstream oss;
                    oss << std::scientific << val;
                    return oss.str();
                }
            }
            else
            {
                cerr << "CRITICAL ERROR: JSON_IO.h: JSONVal::to_string(): Could not identify input." << endl;
                return "<<ERROR>>";
            }
        }, *this);
    }
};







// Define a basic structure for the JSON dictionaries
struct JSONDict : map<string, JSONVal>
{
    // Attach a vector of numbers as dtype "char" for loading JSONDict files 
    vector<char> numbers_char = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', '.'};



    // Function to get string values (as a unique_ptr<string>) from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, std::unique_ptr<string> &external_value)
    {
        // Check the strings
        if (*this->find(key) != *this->end())
        {
            string arg = (string)(*this)[key];
            external_value = std::make_unique<string>(arg);
            return;
        }
    }
    // Function to get string values (as a shared_ptr<string>) from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, std::shared_ptr<string> &external_value)
    {
        // Check the strings
        if (*this->find(key) != *this->end())
        {
            string arg = (string)(*this)[key];
            external_value.reset();
            external_value = std::make_shared<string>(arg);
            return;
        }
    }
    // Function to get string values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, string &external_value)
    {
        // Check the strings
        if (*this->find(key) != *this->end())
        {
            external_value = (string)(*this)[key];
            return;
        }
    }
    // Function to get vector<string> values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, vector<string> &external_value)
    {
        // Check the strings
        //if (*this->find(key) != *this->end())
        if ((*this).find(key) != (*this).end())
        {
            external_value = (vector<string>)(*this)[key];
            return;
        }
    }
    // Function to get int values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, int &external_value)
    {
        // Check the ints
        if (*this->find(key) != *this->end())
        {
            external_value = (int)(*this)[key];
            return;
        }
    }
    // Function to get vector<int> values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, vector<int> &external_value)
    {
        // Check the ints
        if (*this->find(key) != *this->end())
        {
            external_value = (vector<int>)(*this)[key];
            return;
        }
    }
    // Function to get vector<vector<int>> values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, vector<vector<int>> &external_value)
    {
        // Check the ints
        //if (*this->find(key) != *this->end())
        if ((*this).find(key) != (*this).end())
        {
            external_value = (vector<vector<int>>)(*this)[key];
            return;
        }
    }
    // Function to get double values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, double &external_value)
    {
        // Check the doubles
        if (*this->find(key) != *this->end())
        {
            external_value = (double)(*this)[key];
            return;
        }
    }
    // Function to get vector<double> values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, vector<double> &external_value)
    {
        // Check the doubles
        //if (*this->find(key) != *this->end())
        if ((*this).find(key) != (*this).end()) // Allows for returning nothing if key is not in "this"
        {
            external_value = (vector<double>)(*this)[key];
            return;
        }
    }
    // Function to get JSONDict values from a dictionary if the provided key is in the dictionary
    void getValue(const string &key, JSONDict &external_value)
    {
        // Check the JSONDicts
        //if (this->find(key) != this->end())
        if ((*this).find(key) != (*this).end())
        {
            external_value = *(*this)[key];
            return;
        }
    }
    void getValue(const string &key, JSONDict &external_value) const
    {
        // Check the JSONDicts
        //if (this->find(key) != this->end())
        if ((*this).find(key) != (*this).end())
        {
            external_value = *(*this).at(key);
            return;
        }
    }





    // Function for merging JSONDicts. Information in the same key path can be either
    // overwritten or not.
    void merge(const JSONDict &dict, bool overwrite = true)
    {
        // Get the keys from the current dictionary layer of "this" and the input JSONDict
        vector<string> this_keys, dict_keys;
        for (const auto& pair : (*this)) { this_keys.push_back(pair.first); }
        for (const auto& pair : dict) { dict_keys.push_back(pair.first); }

        // Look through the keys for a match; if a match is found,
        // go to the next layer in the JSONDict and repeat
        bool found_match;
        for (const string dict_key : dict_keys)
        {
            found_match = false;
            for (const string this_key : this_keys)
            {
                if (dict_key == this_key)
                {
                    found_match = true;
                    
                    // See if the object in each dictionary at the key are JSONDict
                    JSONDict check_dict_this, check_dict_dict;
                    this->getValue(this_key, check_dict_this); dict.getValue(dict_key, check_dict_dict);
                    if (check_dict_this.size() > 0 && check_dict_dict.size() > 0)
                    {
                        // If both objects are JSONDict, then merge them and replace the value in "this" with the merged JSONDict
                        check_dict_this.merge(check_dict_dict);
                        JSONDict* check_dict_this_pntr = new JSONDict(check_dict_this);
                        (*this)[this_key] = check_dict_this_pntr;
                    }
                    else
                    {
                        // If not, respect the overwrite option
                        if (overwrite) { (*this)[this_key] = dict.at(dict_key); }
                    }
                    break;
                }
            }

            // If no matching keys are found, merge the dict_key/associated information to this
            if (!found_match) { (*this)[dict_key] = dict.at(dict_key); }
        }
    }
    




    // Function for saving JSONDicts to a stream
    istringstream saveToStream()
    {
        istringstream stream( (*this).to_string() );
        return stream;
    }

    // Function for saving JSONDicts to a file
    void saveToFile(const string &path)
    {
        ofstream file(path);
        file << (*this).to_string();
        file.close();
    }

    // Function for converting JSONDict variables to a string
    string to_string(string indent = "    ", int level = 0) const
    {
        // Create output stream
        ostringstream oss;
        
        // Open the brackets
        oss << "{\n";

        // Iterate the indent level and create the indentation for the dictionary keys
        level += 1;
        string key_indent = "";
        for (int i = 0; i < level; i++) { key_indent += indent; }
        
        // Convert the arguments of all keys to strings
        size_t count = 0;
        for (const auto& [key, value] : *this)
        {
            oss << key_indent << "\"" << key << "\": " << value.to_string(indent, level);
            if (++count < this->size()) oss << ",\n";
        }

        // Create indentation for the closing brackets, and close the brackets
        string bracket_indent = "";
        for (int i = 0; i < level - 1; i++) { bracket_indent += indent; }
        oss << "\n" << bracket_indent << "}";

        // Return the string
        return oss.str();
    }




    // Function for loading JSONDicts and JSONVals from a file
    int loadFromStream(istringstream &stream)
    {
        string line;
        
        while (getline(stream, line))
        {
            // Get the position of the first non-space character
            size_t first_char_pos = GetFirstNonSpacePos(line);
            
            // Skip blank lines if there are any
            if (first_char_pos == string::npos) { continue; }

            // Skip c++ like comments (i.e., those that start with "//")
            if (line[first_char_pos] == '/' && line[first_char_pos + 1] == '/') { continue; }

            // Entering the structure
            if (line[first_char_pos] == '{')
            {
                JSONDict k = CreateDict(stream);
                *this = k;
                return 0;
            }
            else
            {
                cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.loadFromStream(): Unable to identify first char of the file in line: " << line << endl;
                return 1;
            }
        }
        cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.loadFromStream(): Unable to find first char of the file." << endl;
        return 1;
    }
    

    // Function for loading JSONDicts and JSONVals from a file
    int loadFromFile(const string &filename)
    {
        // Find the file
        ifstream file(filename);
        
        if (!file) { cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.loadFromFile(): Error opening file for reading." << endl; return 1; }

        string line;
        
        while (getline(file, line))
        {
            // Get the position of the first non-space character
            size_t first_char_pos = GetFirstNonSpacePos(line);
            
            // Skip blank lines if there are any
            if (first_char_pos == string::npos) { continue; }

            // Skip c++ like comments (i.e., those that start with "//")
            if (line[first_char_pos] == '/' && line[first_char_pos + 1] == '/') { continue; }

            // Entering the structure
            if (line[first_char_pos] == '{')
            {
                JSONDict k = CreateDict(file);
                file.close();
                *this = k;
                return 0;
            }
            else
            {
                cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.loadFromFile(): Unable to identify first char of the file in line: " << line << endl;
                file.close();
                return 1;
            }
        }
        cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.loadFromFile(): Unable to find first char of the file." << endl;
        file.close();
        return 1;
    }

    size_t GetFirstNonSpacePos(string &line, size_t start_pos = 0)
    {
        size_t pos = start_pos;
        while (pos < line.size() && line[pos] == ' ') { ++pos; }
        if (pos < line.size()) { return pos; }
        else { return string::npos; }
    }

    //JSONDict CreateDict(ifstream &file)
    JSONDict CreateDict(istream &file)
    {
        JSONDict k;
        string line;
            
        while (getline(file, line))
        {
            // Get the position of the first non-space character
            size_t pos = GetFirstNonSpacePos(line);
            
            // Skip blank lines if there are any
            if (pos == string::npos) { continue; }

            // Skip c++ like comments (i.e., those that start with "//")
            if (line[pos] == '/' && line[pos + 1] == '/') { continue; }

            // If the first character indicates a string...
            if (line[pos] == '\"')
            {
                // Read the string
                string key = ReadString(line, pos); // second argument is the pos of the first quote; it will be updated to the pos of the second quote
                pos += 3; // Skip the ": " to get to the first position of the argument
                GetAndSetArg(k, key, file, line, pos);
            }
            else if (line[pos] == '}')
            {
                return k;
            }
            else
            {
                cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.CreateDict(): Unable to identify the first char of the key in line: " << line << endl;
                string key = "<<ERROR IDENTIFYING KEY>>";
                k[key] = "<<ERROR IDENTIFYING ARGUMENT>>";
            }
        }
        cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.CreateDict(): The code shouldn't ever really get here. Line is: " << line << endl;
        return k; // Shouldn't ever get here
    }

    //void GetAndSetArg(JSONDict &k, const string &key, ifstream &file, const string &line, size_t pos)
    void GetAndSetArg(JSONDict &k, const string &key, istream &file, const string &line, size_t pos)
    {
        // Check if the first character is a number
        bool isNumber = isCharInVector(line[pos], numbers_char);

        if (isNumber) // If the first character indicates a number
        {
            bool isDouble;
            string number = GetNumberAsString(line, pos, isDouble);
            if (isDouble) { k[key] = stod(number); }
            else { k[key] = stoi(number); }
        }
        else if (line[pos] == '\"') // If the first character indicates a string
        {
            string str = ReadString(line, pos);
            k[key] = str;
        }
        else if (line[pos] == '[' && line[pos + 1] == '[') // If the first 2 characters indicate a nested vector
        {
            vector<vector<JSONVal>> nestedVec = GetNestedVector(line, pos);
            k[key] = nestedVec;
        }
        else if (line[pos] == '[') // If the first character indicates a vector
        {
            vector<JSONVal> vec = GetVector(line, pos);
            k[key] = vec;
        }
        else if (line[pos] == '{') // If the first character indicates a dictionary
        {
            JSONDict dict = CreateDict(file);
            JSONDict* kk = new JSONDict(dict);
            k[key] = kk;
        }
        else
        {
            cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.GetAndSetArg(): Unable to identify the first char of the argument in line: " << line << endl;
            k[key] = "<<ERROR IDENTIFYING ARGUMENT>>";
        }
    }

    string ReadString(const string &line, size_t &first_quote_pos)
    {
        size_t first_letter_pos = first_quote_pos + 1; // First letter pos is first quote pos + 1
        size_t last_quote_pos = line.find("\"", first_letter_pos);
        if (last_quote_pos == string::npos)
        {
            cerr << "CRITICAL ERROR: JSON_IO.h: JSONDict.ReadString(): Unable to find the second quote for string definition in line: " << line << endl;
            return "<<ERROR IDENTIFYING STRING>>";
        }
        first_quote_pos = last_quote_pos; // Update first_quote_pos to be last_quote_pos
        return line.substr(first_letter_pos, last_quote_pos - first_letter_pos);
    }

    string GetNumberAsString(const string &line, size_t &pos_first_char, bool &isDouble)
    {
        string number = "";
        isDouble = false;
        size_t pos = pos_first_char;

        vector<char> numbers_char_plus_e = numbers_char;
        numbers_char_plus_e.push_back( '+' );
        numbers_char_plus_e.push_back( 'e' );
        
        while(isCharInVector(line[pos], numbers_char_plus_e))
        {
            number += line[pos];
            if (line[pos] == '.') { isDouble = true; }
            pos += 1;
        }
        pos -= 1;
        pos_first_char = pos; // pos_first_char is updated with the pos of the last char
        return number;
    }

    vector<vector<JSONVal>> GetNestedVector(const string &line, size_t &pos_first_bracket)
    {
        vector<vector<JSONVal>> nestedVec;
        size_t pos = pos_first_bracket + 1;
        
        while(line[pos] != ']')
        {
            nestedVec.push_back( GetVector(line, pos) );
            pos += 1; // Get the pos just after the inner vector. Could be ']' on the outer vector, or ','
            if (line[pos] == ',') { pos += 2; } // If it is ',', add 2 to skip the space after the comma and land on '[' of the next inner vector
        }
        pos_first_bracket = pos; // Update pos_first_bracket to be at the position of ']' of the outer vector
        return nestedVec;
    }

    vector<JSONVal> GetVector(const string &line, size_t &pos_first_bracket)
    {
        vector<JSONVal> vec;
        size_t pos = pos_first_bracket + 1;
        
        while(line[pos] != ']')
        {
            if (line[pos] == '\"')
            {
                string str = ReadString(line, pos);
                vec.push_back(str);
            }
            else
            {
                bool isDouble;
                string number = GetNumberAsString(line, pos, isDouble);
                if (isDouble) { vec.push_back(stod(number)); }
                else { vec.push_back(stoi(number)); }
            }
            pos += 1; // Get the pos just after the vector argument. Could be ']' or ','
            if (line[pos] == ',') { pos += 2; } // If it is ',', add 2 to skip the space after the comma and land on the first char of the next vector arg
        }
        pos_first_bracket = pos; // Update pos_first_bracket to be at the position of ']'
        return vec;
    }
    
    bool isCharInVector(const char &test_char, const vector<char> &vec)
    {
        for (int i = 0; i < vec.size(); i++) { if (test_char == vec[i]) { return true; } }
        return false;
    }
};

