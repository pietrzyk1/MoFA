#!/bin/bash

#lL.lH.
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
#lL.lH.

#############################
#### Shell script set up ####
#############################

# Tell shell script to exit immediately if a command fails
set -e

# Get this shell script's name
SH_NAME=$(basename "$0")

# Get the directory that this shell script is in
SH_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"



##################################
#### Define utility functions ####
##################################

# This function extracts strings from the JSON-like dictionaries/files used to configure the MoFA
# simulations. It can be used in the following way:
#
# Given a bash array of strings, which represent some dictionary keys of a JSON-like dictionary/file, such as
#
#     input_keys=("mesh" "info path" "file name")
#
# this function can be used like
#
#     read -r var < <(extract_string config.txt "${input_keys[@]}")
#
# to get "var", the string stored at config.txt["mesh"]["info path"]["file name"].
#
extract_string()
{
    local file="$1"
    shift
    local -a keys=("$@")  # all remaining args are keys in order

    local depth=0
    local max_depth=${#keys[@]}
    local line
    local -a extracted_keys

    while IFS= read -r line; do
        # Check if current line contains expected key at current depth
        if echo "$line" | grep -q "\"${keys[depth]}\""; then
            ((++depth))
            # If we've reached last key, extract the array from this line
            if [[ $depth -eq $max_depth ]]; then
                # Extract array inside brackets
                local raw=$(echo "$line" | sed -E 's/.*:\s*(.*),?/\1/' | tr -d '" ')
                raw="${raw%,}"
                echo "$raw"
                return
            fi
            continue
        fi
    done < "$file"

    echo "$SH_NAME: extract_string(): Key path not found: ${keys[*]}" >&2
    return 1
}


# This function extracts vectors of strings from the JSON-like dictionaries/files used to configure the MoFA
# simulations. It can be used in the following way:
#
# Given a bash array of strings, which represent some dictionary keys of a JSON-like dictionary/file, such as
#
#     input_keys=("stokes closure" "residual info" "simulation keys")
#
# this function can be used like
#
#     mapfile -t keys < <(extract_string_vector config.txt "${input_keys[@]}")
#
# to get "keys", the vector of strings stored at config.txt["stokes closure"]["residual info"]["simulation keys"].
#
extract_string_vector()
{
    local file="$1"
    shift
    local -a keys=("$@")  # all remaining args are keys in order

    local depth=0
    local max_depth=${#keys[@]}
    local line
    local -a extracted_keys

    while IFS= read -r line; do
        # Check if current line contains expected key at current depth
        if echo "$line" | grep -q "\"${keys[depth]}\""; then
            ((++depth))
            # If we've reached last key, extract the array from this line
            if [[ $depth -eq $max_depth ]]; then
                # Extract array inside brackets
                local raw=$(echo "$line" | sed -E 's/.*\[(.*)\].*/\1/' | tr -d '" ')
                IFS=',' read -ra extracted_keys <<< "$raw"
                for key in "${extracted_keys[@]}"; do
                    echo "$key"
                done
                return
            fi
            continue
        fi
    done < "$file"

    echo "$SH_NAME: extract_string_vector(): Key path not found: ${keys[*]}" >&2
    return 1
}


# This function extracts numbers from the JSON-like dictionaries/files used to configure the MoFA
# simulations. It can be used in the following way:
#
# Given a bash array of strings, which represent some dictionary keys of a JSON-like dictionary/file, such as
#
#     MAX_ITER_KEYPATH=("stokes closure" "closure parameters" "max recursion iterations")
#
# this function can be used like
#
#     read -r var < <(extract_number config.txt "${MAX_ITER_KEYPATH[@]}")
#
# to get "var", the number stored at config.txt["stokes closure"]["closure parameters"]["max recursion iterations"].
#
extract_number()
{
    local file="$1"
    shift
    local -a keys=("$@")  # all remaining args are keys in order

    local depth=0
    local max_depth=${#keys[@]}
    local line
    local -a extracted_keys

    while IFS= read -r line; do
        # Check if current line contains expected key at current depth
        if echo "$line" | grep -q "\"${keys[depth]}\""; then
            ((++depth))
            # If we've reached last key, extract the array from this line
            if [[ $depth -eq $max_depth ]]; then
                # Extract number (can be a negative/positive integer/double, and in scientific notation)
                local raw=$(echo "$line" | sed -E 's/.*:\s*(-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?),?/\1/')
                echo "$raw"
                return
            fi
            continue
        fi
    done < "$file"

    echo "$SH_NAME: extract_number(): Key path not found: ${keys[*]}" >&2
    return 1
}


# This function extracts vectors of numbers from the JSON-like dictionaries/files used to configure the MoFA
# simulations. It can be used in the following way:
#
# Given a bash array of strings, which represent some dictionary keys of a JSON-like dictionary/file, such as
#
#     input_keys=("stokes closure" "closure parameters" "u closure variable components")
#
# this function can be used like
#
#     mapfile -t keys < <(extract_number_vector config.txt "${input_keys[@]}")
#
# to get "keys", the vector of numbers stored at config.txt["stokes closure"]["closure parameters"]["u closure variable components"].
#
extract_number_vector()
{
    local file="$1"
    shift
    local -a keys=("$@")  # all remaining args are keys in order

    local depth=0
    local max_depth=${#keys[@]}
    local line
    local -a extracted_keys

    while IFS= read -r line; do
        # Check if current line contains expected key at current depth
        if echo "$line" | grep -q "\"${keys[depth]}\""; then
            ((++depth))
            # If we've reached last key, extract the array from this line
            if [[ $depth -eq $max_depth ]]; then
                # Extract array inside brackets
                local raw=$(echo "$line" | sed -E 's/.*\[(.*)\].*/\1/' | tr -d ' ')
                IFS=',' read -ra extracted_keys <<< "$raw"
                for key in "${extracted_keys[@]}"; do
                    echo "$key"
                done
                return
            fi
            continue
        fi
    done < "$file"

    echo "$SH_NAME: extract_number_vector(): Key path not found: ${keys[*]}" >&2
    return 1
}


# This function writes a text file (meant to be temporary) that contains a JSON-like dictionary used to configure the MoFA
# simulations. It can be used in the following way:
#
# Given bash arrays like
#
#     KEYS=("key1" "key2")
#     VALUES=("[\"string1\", \"string2\"]" "[[1, 2, 3], [4, 5, 6]]")
#
# this function can be used like
#
#     write_temporary_JSONDict "output_file.txt" KEYS VALUES
#
# to write a short text file containing a JSON-like dictionary with the provided keys and values.
#
write_temporary_JSONDict()
{
    local file="$1"
    shift
    local -n keys=$1   # name reference (requires Bash 4.3+)
    shift
    local -n values=$1
    shift

    if [[ ${#keys[@]} -ne ${#values[@]} ]]; then
        echo "$SH_NAME: write_temporary_JSONDict(): Number of keys and values must match"
        return 1
    fi

    echo "{" > "$file"

    local last_index=$((${#keys[@]} - 1))
    for i in "${!keys[@]}"; do
        local key="${keys[i]}"
        local value="${values[i]}"

        # Detect if value looks like a JSON array or object (already quoted)
        if [[ "$value" =~ ^\[.*\]$ || "$value" =~ ^\{.*\}$ ]]; then
            line="    \"$key\": $value"
        else
            line="    \"$key\": \"$value\""
        fi

        # Add comma at the end unless it's the last item
        if [[ $i -lt $last_index ]]; then
            line="$line,"
        fi

        echo "$line" >> "$file"
    done

    echo "}" >> "$file"
    #echo "$SH_NAME: write_temporary_JSONDict(): Created $file"
}


# This function takes in an array of string arguments and converts it to a string of an array that contains the original
# string arguments. The final "string array" is formatted for the JSON-like dictionaries used to configure the MoFA
# simulations. This function can be used in the following way:
#
# Given a bash array of strings, such as
#
#     NUMBERS=("1" "2" "3")
#
# or
#
#     STRINGS=("key1" "key2" "key3")
#     
# this function can be used like
#
#     ARRAY1=$(array_to_string NUMBERS 0)
#
# or
#
#     ARRAY2=$(array_to_string STRINGS 1)
#
# to create strings of the original arrays (i.e., ARRAY1 = "[1, 2, 3]" and ARRAY2 = "["key1", "key2", "key3"]". Note that
# the last arguments (0 and 1) add or ignore (respectively) quotes on the components. If no number is provided, it
# defaults to 0 (no quotes).
#
array_to_string()
{
    local -n arr_ref=$1  # Reference the array by name
    local quote_toggle=${2:-0}
    local result="["

    for i in "${!arr_ref[@]}"; do
        if [ $quote_toggle -eq 1 ]; then
            result+="\"${arr_ref[i]}\""
        else
            result+="${arr_ref[i]}"
        fi
        if [[ $i -lt $((${#arr_ref[@]} - 1)) ]]; then
            result+=", "
        fi
    done

    result+="]"
    echo "$result"
}

