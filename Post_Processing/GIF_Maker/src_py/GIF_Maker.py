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

import sys
import os
import code
from PIL import Image


# Get input arguments (sys.argv[0] is the script name)
save_dir = sys.argv[1] # directory with the png files
PNG_file_prefix = sys.argv[2] # prefix of the files


# Navigate to the directory with the png files and get the names of the files in the directory
os.chdir(save_dir)
directory_file_names = os.listdir()


# Create the frames
frames = []
imgs = []
for i in range(len(directory_file_names)):
	PNG_file_name = PNG_file_prefix + str(i) + '.png'
	flag = 0
	for file in directory_file_names:
		if file == PNG_file_name:
			flag = 1
			imgs.append(file)
			break
	if flag == 0:
		break

for i in imgs:
	new_frame = Image.open(i)
	frames.append(new_frame)


# Save into a GIF file that loops forever
frames[0].save(PNG_file_prefix + 'GIF.gif', format='GIF', append_images=frames[1:], save_all=True, duration=60, loop=0)
