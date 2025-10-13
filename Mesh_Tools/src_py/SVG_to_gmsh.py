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

import numpy as np
import code


def main():

    # Provide the SVG directory and file name
    geometry_dir = '../../../Meshes/data/Single_Cylinder/11x1_D0p5/'
    geometry_file_name_prefix = 'Single_Cylinder'
    geometry_file_name_suffix = '.svg'

    # Provide the save file directory and file name
    save_dir = geometry_dir
    save_file_name_prefix = geometry_file_name_prefix
    save_file_name_suffix = '.txt'
    

    svgReader = SVGReader(geometry_dir + geometry_file_name_prefix + geometry_file_name_suffix)
    geometry_points = svgReader.GetPointsFromFile() # Format is geometry_points[surface][x or y][coord data]

    with open(save_dir + save_file_name_prefix + save_file_name_suffix, 'w') as file:
        file.write('{\n')
        
        # Write the geometry
        file.write('    "geometry": {\n')
        for i_surf in range(len(geometry_points)):
            file.write('        "surface_' + str(i_surf) + '": {\n')

            x_points = '['
            y_points = '['
            for i_pnt in range(len(geometry_points[i_surf][0]) - 1): # The data will produce the same point at the beginning and ending. We just get rid of this here.
                x_points += str(geometry_points[i_surf][0][i_pnt])
                y_points += str(geometry_points[i_surf][1][i_pnt])
                if i_pnt < len(geometry_points[i_surf][0]) - 2:
                    x_points += ', '
                    y_points += ', '
                else:
                    x_points += ']'
                    y_points += ']'
            
            file.write('            "x_coords": ' + x_points + ',\n')
            file.write('            "y_coords": ' + y_points + '\n')
            file.write('        }')
            if i_surf < len(geometry_points) - 1:
                file.write(',')
            file.write('\n')
        file.write('    },\n')

        # Write the additional information
        file.write('    "other": {\n')
        file.write('        "N_surfaces": ' + str(len(geometry_points)) + '\n')
        file.write('    }\n')

        file.write('}')
    
    print("Successfully transformed the SVG into point data for gmsh!")





class SVGReader:
    def __init__(self, geometry_file_name, ):

        self.geometry_file_name = geometry_file_name

        self.stringNumbers = ['0','1','2','3','4','5','6','7','8','9','-','.']


        self.f = open(geometry_file_name, 'r')
        self.content = self.f.readlines()



    def GetPointsFromFile(self):

        path_data = []
        for i_line in range(len(self.content)):
            
            if '<path' in self.content[i_line]:
                for i_d in range(len(self.content) - i_line):
                    
                    if 'd="' in self.content[i_line + i_d] and not 'id="' in self.content[i_line + i_d]:
                        line = self.content[i_line + i_d]
                        
                        # Find the first letter index
                        for i_letter in range(len(line)):
                            if line[i_letter] == '"':
                                letter_ind = i_letter + 1
                                break

                        terminate_symbol = line[letter_ind]
                        d_paths = 0
                        while terminate_symbol != '"':
                            
                            if terminate_symbol == 'm':
                                letter_ind += 2
                                path_data.append( [[], []] )
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][1].append( number )
                                if d_paths > 0:
                                    path_data[-1][0][-1] = path_data[-1][0][-1] + path_data[-2][0][-1]
                                    path_data[-1][1][-1] = path_data[-1][1][-1] + path_data[-2][1][-1]
                            
                            elif terminate_symbol == 'M': # IF THERE ARE ISSUES, CONSIDER THIS ELIF; IT IS NEW AND UNTESTED
                                letter_ind += 2
                                path_data.append( [[], []] )
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][1].append( number )
                        
                            elif terminate_symbol == 'l' or self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                if terminate_symbol == 'l':
                                    letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( path_data[-1][0][-1] + number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][1].append( path_data[-1][1][-1] + number )
                            
                            elif terminate_symbol == 'v':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( path_data[-1][0][-1] )
                                path_data[-1][1].append( path_data[-1][1][-1] + number )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'v'
                                    letter_ind -= 2

                            elif terminate_symbol == 'V':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( path_data[-1][0][-1] )
                                path_data[-1][1].append( number )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'V'
                                    letter_ind -= 2
                            
                            elif terminate_symbol == 'h':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( path_data[-1][0][-1] + number )
                                path_data[-1][1].append( path_data[-1][1][-1] )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'h'
                                    letter_ind -= 2

                            elif terminate_symbol == 'H':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( number )
                                path_data[-1][1].append( path_data[-1][1][-1] )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'H'
                                    letter_ind -= 2
                            
                            elif terminate_symbol == 'c':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                    
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( path_data[-1][0][-1] + number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][1].append( path_data[-1][1][-1] + number )
                                    
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'c'
                                    letter_ind -= 2

                            elif terminate_symbol == 'C':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                    
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][0].append( number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind)
                                path_data[-1][1].append( number )
                                    
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'C'
                                    letter_ind -= 2

                            elif terminate_symbol == 'z' or terminate_symbol == 'Z':
                                d_paths += 1
                                if line[letter_ind + 1] == '"':
                                    terminate_symbol = '"'
                                else:
                                    letter_ind += 2
                                    terminate_symbol = line[letter_ind]
                                
                        break

        return path_data 



    def ReadNumber(self, line, first_index):
        number = ''
        for i in range(first_index, len(line)):
            if line[i] == ',':
                terminate_symbol = line[i]
                break
            elif line[i] == ' ':
                i += 1
                terminate_symbol = line[i]
                break
            elif line[i] == '"':
                terminate_symbol = line[i]
                break
            else:
                number += line[i]
        
        return float(number), i, terminate_symbol

    def CompareNumbers(self, number, numberList):
        numberBool = False
        for i in range(len(numberList)):
            if numberList[i] == number:
                numberBool = True
                break
        return numberBool



if __name__ == "__main__":
    main()
