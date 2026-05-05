# lL.lH.
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
# lL.lH.

import numpy as np
import code


def main():

    # Provide the SVG directory and file name
    #geometry_dir = '../../../Meshes/data/Single_Cylinder/11x1_D0p5/'
    #geometry_dir = '../../../Meshes/data/Sintered_Iron_2D_Periodic_UnitCell/'
    geometry_dir = '../../../Meshes/data/Sand_Stone_1/'
    #geometry_file_name_prefix = 'Single_Cylinder'
    #geometry_file_name_prefix = 'Sintered_Iron_2D_Periodic_UnitCell'
    geometry_file_name_prefix = 'Sand_Stone_1'
    geometry_file_name_suffix = '.svg'

    # Provide the save file directory and file name
    save_dir = geometry_dir
    save_file_name_prefix = geometry_file_name_prefix
    save_file_name_suffix = '.txt'
    
    # Provide the number of places you would like after the decimal
    places_after_deci = 3


    svgReader = SVGReader(geometry_dir + geometry_file_name_prefix + geometry_file_name_suffix)
    geometry_points = svgReader.GetPointsFromFile(places_after_deci + 1) # Format is geometry_points[surface][x or y][coord data]
    with open(save_dir + save_file_name_prefix + save_file_name_suffix, 'w') as file:
        file.write('{\n')
        
        # Write the geometry
        file.write('    "geometry": {\n')
        for i_surf in range(len(geometry_points)):
            file.write('        "surface_' + str(i_surf) + '": {\n')

            x_points = '['
            y_points = '['
            for i_pnt in range(len(geometry_points[i_surf][0]) - 1): # The data will produce the same point at the beginning and ending. We just get rid of this here.
                
                if len(geometry_points[i_surf][0][i_pnt]) > geometry_points[i_surf][0][i_pnt].index(".") + places_after_deci:
                    geometry_points[i_surf][0][i_pnt] = geometry_points[i_surf][0][i_pnt][: geometry_points[i_surf][0][i_pnt].index(".") + places_after_deci + 1]
                if len(geometry_points[i_surf][1][i_pnt]) > geometry_points[i_surf][1][i_pnt].index(".") + places_after_deci:
                    geometry_points[i_surf][1][i_pnt] = geometry_points[i_surf][1][i_pnt][: geometry_points[i_surf][1][i_pnt].index(".") + places_after_deci + 1]

                x_points += geometry_points[i_surf][0][i_pnt]
                y_points += geometry_points[i_surf][1][i_pnt]
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



    def GetPointsFromFile(self, round_places):
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
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][1].append( number )
                                if d_paths > 0:
                                    path_data[-1][0][-1] = self.AddStrNumbers(path_data[-1][0][-1], path_data[-2][0][-1], round_places)
                                    path_data[-1][1][-1] = self.AddStrNumbers(path_data[-1][1][-1], path_data[-2][1][-1], round_places)
                            
                            elif terminate_symbol == 'M': # IF THERE ARE ISSUES, CONSIDER THIS ELIF; IT HAS NOT BEEN TESTED, BUT ASSUMED TO WORK
                                letter_ind += 2
                                path_data.append( [[], []] )
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][1].append( number )
                        
                            elif terminate_symbol == 'l' or self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                if terminate_symbol == 'l':
                                    letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( self.AddStrNumbers(path_data[-1][0][-1], number, round_places) )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][1].append( self.AddStrNumbers(path_data[-1][1][-1], number, round_places) )
                            
                            elif terminate_symbol == 'v':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( path_data[-1][0][-1] )
                                path_data[-1][1].append( self.AddStrNumbers(path_data[-1][1][-1], number, round_places) )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'v'
                                    letter_ind -= 2

                            elif terminate_symbol == 'V':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( path_data[-1][0][-1] )
                                path_data[-1][1].append( number )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'V'
                                    letter_ind -= 2
                            
                            elif terminate_symbol == 'h':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( self.AddStrNumbers(path_data[-1][0][-1], number, round_places) )
                                path_data[-1][1].append( path_data[-1][1][-1] )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'h'
                                    letter_ind -= 2

                            elif terminate_symbol == 'H':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( number )
                                path_data[-1][1].append( path_data[-1][1][-1] )
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'H'
                                    letter_ind -= 2
                            
                            elif terminate_symbol == 'c':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                    
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( self.AddStrNumbers(path_data[-1][0][-1], number, round_places) )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][1].append( self.AddStrNumbers(path_data[-1][1][-1], number, round_places) )
                                    
                                if self.CompareNumbers(terminate_symbol, self.stringNumbers):
                                    terminate_symbol = 'c'
                                    letter_ind -= 2

                            elif terminate_symbol == 'C':
                                letter_ind += 2
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                    
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
                                path_data[-1][0].append( number )
                                letter_ind += 1
                                number, letter_ind, terminate_symbol = self.ReadNumber(line, letter_ind, round_places)
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



    def ReadNumber(self, line, first_index, round_places):
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
        
        if len(number) == 0:
            return number, i, terminate_symbol

        '''
        if "e" in number:
            after_e = number[number.index("e") + 1]
            if number[number.index("e") + 1] != "-":
                return number, i, terminate_symbol
            if int(number[number.index("e") + 1 :]) <= -round_places - 2:
                return "0", i, terminate_symbol

            if int(number[number.index("e") + 1 :]) == -round_places:
                print("Problems")
                code.interact(local=locals())

            if int(number[number.index("e") + 1 :]) == -round_places - 1:
                if number[0] == "-":
                    neg = "-"
                    pre_decimal = number[1]
                else:
                    neg = ""
                    pre_decimal = number[0]

                if abs(int(neg + pre_decimal)) >= 5:
                    return neg + str(10**(-round_places)), i, terminate_symbol
                else:
                    return "0", i, terminate_symbol

        for j in range(len(number)):
            if number[j] == '.' and len(number) >= j + 1 + (round_places + 1):
                if int(number[j + (round_places + 1)]) >= int('5'):
                    if int(number[j + round_places]) + 1 == 10:
                        new_num = str(int(number[j + round_places - 1]) + 1)
                        number = number[0 : j + round_places - 1] + new_num + "0"
                    else:
                        number = number[0 : j + round_places] + str(int(number[j + round_places]) + 1)
                else:
                    number = number[0 : j + round_places + 1]
                break
        '''

        number = self.RoundStrNumber(number, round_places)
        
        return number, i, terminate_symbol
    


    def AddStrNumbers(self, a, b, round_places):        
        number = self.RoundStrNumber(str(float(a) + float(b)), round_places)
        return number



    def RoundStrNumber(self, number, round_places):
        # See if the number has an 'e' in it
        if "e" in number:
            # Get the exponent to 10 (i.e., 1.0eNUM, get NUM)
            exp_ind = number.index("e")
            exponent = int(number[exp_ind + 1 :])
            
            # If the exponent is positive, it is likely that the number is big and does not need to be rounded
            if exponent > 0:
                print("Case not considered yet")
                code.interact(local=locals())
            # If the exponent is too small, then return 0
            if exponent <= -round_places - 2:
                return "0.0"

            # Make sure there is a decimal
            if "." not in number:
                number = number[:exp_ind] + ".0" + number[exp_ind:]    
            deci_ind = number.index(".")

            before_deci = number[: deci_ind]
            neg = ""
            if before_deci[0] == "-":
                neg = "-"
                before_deci = before_deci[1:]
            after_deci = number[deci_ind + 1 : exp_ind]
            assert(round_places >= len(before_deci))
            
            number = neg + "0." + "0"*(round_places - len(before_deci)) + before_deci + after_deci
            

            '''
            # Hopefully this does not happen... it will require more code
            if exponent == -round_places:
                print("Case not considered yet")
                code.interact(local=locals())

            if exponent == -round_places - 1:
                if number[0] == "-":
                    neg = "-"
                    pre_decimal = number[1]
                else:
                    neg = ""
                    pre_decimal = number[0]

                if abs(int(neg + pre_decimal)) >= 5:
                    return neg + str(10**(-round_places))
                else:
                    return "0.0"
            '''
        
        foundDecimal = False
        for j in range(len(number)):
            if number[j] == '.':
                foundDecimal = True

            if number[j] == '.' and len(number) >= j + 1 + (round_places + 1):
                if int(number[j + (round_places + 1)]) >= int('5'):
                    new_after_deci = str(int(number[j + 1 : j + (round_places + 1)]) + 1)
                    if len(new_after_deci) <= round_places:
                        number = number[: j + 1] + "0"*(round_places - len(new_after_deci)) + new_after_deci
                    else:
                        if number[0] == "-":
                            first_num = str(int(number[: j]) - int(new_after_deci[0]))
                        else:
                            first_num = str(int(number[: j]) + int(new_after_deci[0]))
                        number = first_num + "." + new_after_deci[1:]
                    
                    '''
                    if int(number[j + (round_places)]) + 1 == 10:
                        new_num = str(int(number[j + round_places - 1]) + 1)
                        if int(new_num) == 10:
                            new_num = str(int(number[j + round_places - 2]) + 1)
                        
                        number = number[0 : j + round_places - 1] + new_num + "0"
                    else:
                        number = number[0 : j + round_places] + str(int(number[j + round_places]) + 1)
                    '''
                
                else:
                    number = number[0 : j + round_places + 1]
                break
        
        if not foundDecimal:
            number += ".0"

        return number



    def CompareNumbers(self, number, numberList):
        numberBool = False
        for i in range(len(numberList)):
            if numberList[i] == number:
                numberBool = True
                break
        return numberBool



if __name__ == "__main__":
    main()
