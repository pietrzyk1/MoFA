# The Method of Finite Averages (MoFA)



### Description :sparkles:

This repository contains code related to the implementation and validation of The Method of Finite Averages (MoFA), a rigorous upscaling methodology for modeling physical phenomena in heterogeneous porous media. Further details regarding MoFA can be found in the publications cited in the "Relevant Publications" section below.

The files of this repository are the existing code for implementing MoFA (as demonstrated in the publications cited in the
"Relevant Publications" section below) in MFEM. The objectives of this repository are:

1. to develop a parallelized code for implementing MoFA (both ensamble and solver parallelization).
2. to continue the development and extension of MoFA.
3. to provide an open-source repository that makes MoFA implementation more accessible.



### Relevant Publications :newspaper:

1. K. Pietrzyk. "The method of finite averages: a rigorous upscaling methodology for heterogeneous porous media", _Adv. Water Resour._, 188, 104689, 2024.
2. K. Pietrzyk. "The rigorous upscaling of advection-dominated transport in heterogeneous porous media via the method of finite averages", _Adv. Water Resour._, 200, 104958, 2025.
3. S. Xia, I. Battiato, & K. Pietrzyk. "Rigorous upscaling of the Navier-Stokes equations in heterogeneous porous media", (In preparation).



### Project Status :bar_chart:

The MoFA methodology is advancing on multiple fronts:

- Kyle Pietrzyk is extending the MoFA methodology to upscale systems with nonlinearities/nonlinear couplings.

- Thomas Roy is working to run the MFEM code on the LC and run the code with ensamble and solver parallelization.

- Shufan Xia is extending the MoFA methodology to the Navier-Stokes equations, applying MoFA to novel problems
in porous/fractured media and fluid mechanics, and extending the MoFA methodology to systems with coupled fluid flow and transport.




### Repository Contents

The codes in the repository are currently designed to: 

1. Help generate and label meshes (using Gmsh) from a provided geometry. This includes labeling the averaging regions of the mesh.
2. Solve the Stokes equations (for systems forced by inlet-outlet boundary conditions or body forces).
3. Solve the MoFA closure problems to obtain the closure solutions and residuals.
4. Solve the MoFA model.
5. Solve the corresponding pore-scale (i.e., fully-resolved) model.
6. Compare the MoFA and averaged-pore-scale solutions.



## Quick Setup :zap:



### Download MFEM and Gmsh; Clone the MoFA Repo :white_check_mark:

This MoFA repo requires the MFEM and Gmsh (C++) packages to run. It is highly recommended to also download GLVis to easily
view and plot the results. These packages can be downloaded following the steps on their respective webpages.

After downloading, clone the MoFA repo.



### Configure the MoFA Repo with Paths to MFEM and Gmsh :white_check_mark:

After cloning, edit ```Path/to/MFEM/source``` and ```Path/to/Gmsh/source``` in the following lines of
```cmake/Configuration_Paths.cmake``` to the correct paths:

```bash
# Define a cmake variable for the MFEM source directory path. This will be used for all cmake files in the MoFA package.
set(MFEM_SOURCE_DIR "Path/to/MFEM/source")

# Define a cmake variable for the Gmsh source directory path. This will be used for all cmake files in the MoFA package.
set(GMSH_SOURCE_DIR "Path/to/Gmsh/source")
```

<!-- 
The mesh generating codes, solver codes, and post-processing codes can now be built using the ```make``` command in the respective
```build``` directories; however, ```make``` will automatically be called when using the provided shell scripts to run the codes.
As such, we skip "```make```ing" the executables for now.
-->



### Run CMake to Prepare Build Files :white_check_mark:

After configuring the MoFA Repo with paths to MFEM and Gmsh, call ```run_cmake.sh``` in the ```cmake``` directory to create the CMake
build files for each executable in MoFA (if you are in the ```cmake``` directory, enter ```./run_cmake.sh```).



## Create Your First MoFA Model :pencil2:

After downloading MFEM and Gmsh, and configuring the MoFA repo, you are ready to setup your first MoFA model:



### Create a New Project Directory :open_file_folder:

> :exclamation:**TDLR**: From the top level of the MoFA repo, enter
>
> ```
> cp -R Projects/Example_Project/ ../My_First_MoFA_Model
> cd ../My_First_MoFA_Model
> ```

The easiest way to create a new MoFA project directory is to copy the template provided in this repo. From the top level of the
MoFA repo, navigate to the ```Projects``` directory by entering 

```bash
cd Projects
```

Inside, you will see the ```Example_Project``` directory. This directory is meant to help get you
started by containing

1. the directory hierarchy used to organize the configuration files, meshes, and results of a project,
2. shell scripts used to run the various codes provided in the MoFA repo.

Make a copy of this directory called ```My_First_MoFA_Model``` outside the MoFA repo by entering

```bash
cp -R Example_Project/ ../../My_First_MoFA_Model
```

Then, navigate to ```My_First_MoFA_Model``` by entering

```bash
cd ../../My_First_MoFA_Model
```

For the remainder of this tutorial, we will be working out of ```My_First_MoFA_Model```.



### Configure the Shell Scripts in the Project Directory :pencil:

> :exclamation:**TDLR**: in ```My_First_MoFA_Model/scripts/config.sh```, change ```MOFA_DIR="$PROJECT_DIR/../.."``` to
>
> ```
> MOFA_DIR="$PROJECT_DIR/../MoFA"
> ```

Inside of ```My_First_MoFA_Model```, navigate to the ```scripts``` directory:

```bash
cd scripts
```

Here, you will find

1. pre-written shell scripts for running the various MoFA codes,
2. a directory ```subscripts``` for "sub" shell scripts,
3. a file named ```config.sh```.

We will return to the shell scripts later; for now, find the line in ```config.sh``` that reads

```bash
MOFA_DIR="$PROJECT_DIR/../.."
```

Because we copied the ```Example_Project``` directory and moved it to ```My_First_MoFA_Model``` just outside the MoFA repo,
```MOFA_DIR``` is no longer defined as the top level of the MoFA repo. To correct this, change the ```MOFA_DIR``` definition to

```bash
MOFA_DIR="$PROJECT_DIR/../MoFA"
```

Then, save and close the file. The project shell scripts are now ready to run MoFA code.



### Create Input Files :books:

> :exclamation:**TDLR**: from ```My_First_MoFA_Model/scripts/```, enter
>
> ```
> ./make_and_run_input_file_generator.sh -t ut
> ./make_and_run_input_file_generator.sh -t pt
> ```

The last step is to create the project's input files for running the MoFA codes. Because such files can be long and tedious to
write, we have provided code (in ```Input_File_Generator/include```; ```generate_input_file.h``` and ```generate_input_file_TUTORIAL.h```) to automatically generate input files for the tutorial simulation. These "starter" input
files---as well as the provided code---can be modified for running other simulations in the future.

To generate the input files for the tutorial simulation, entering the following commands from the
```My_First_MoFA_Model/scripts``` directory (i.e., we ended the previous step in this directory):

```bash
./make_and_run_input_file_generator.sh -t ut
./make_and_run_input_file_generator.sh -t pt
```

These commands generate JSON-like files ```MoFA_config.txt``` and ```porescale_config.txt``` in the ```My_First_MoFA_Model/config```
directory. The MoFA codes take their commands from these files.



## Run Your First MoFA Model :computer:

> :exclamation:**TDLR**: from ```My_First_MoFA_Model/scripts/```, enter
>
> ```
> ./make_and_run_mesh_generator.sh -t u
> ./make_and_run_stokes_solver.sh -t u
> ./full_make_run_scalar_closure_solver.sh
> ./make_and_run_upscaled_model.sh
> ./make_and_run_mesh_generator.sh -t p
> ./make_and_run_stokes_solver.sh -t p
> ./make_and_run_porescale_solver.sh
> ./make_and_run_error_calc.sh
> ```
>
> The maximum absolute errors will appear in the terminal. They shoule be less than the _a priori_ error threshold, 0.1.

You are now ready to run the MoFA codes from the provided shell scripts. We will first run the shell scripts for formulating and
solving the MoFA model, and then the scripts for solving the fully-resolved model. Finally, we will compare the averaged solutions
from the models to verify that the absolute error is within the expected error threshold.

From the ```My_First_MoFA_Model/scripts``` directory, enter the following ordered commands to perform the associated tasks related
to formulating and solving the MoFA model:

1. Generate the upscaled mesh: ```./make_and_run_mesh_generator.sh -t u```
2. Solve the Stokes equations on the upscaled mesh: ```./make_and_run_stokes_solver.sh -t u```
3. Solve the closure problems: ```./full_make_run_scalar_closure_solver.sh```
4. Solve the upscaled model: ```./make_and_run_upscaled_model.sh```

Then, enter the following ordered commands to perform the associated tasks related to solving the fully-resolved model:

1. Generate the pore-scale mesh: ```./make_and_run_mesh_generator.sh -t p```
2. Solve the Stokes equations on the pore-scale mesh: ```./make_and_run_stokes_solver.sh -t p```
3. Solve the pore-scale model: ```./make_and_run_porescale_solver.sh```

Finally, enter the following command to compute the absolute error between the averaged solutions of the MoFA model and the
fully-resolved model:

1. Compute the absolute error: ```./make_and_run_error_calc.sh```

The maximum absolute error of each averaging region throughout time will appear in the terminal. These errors should all remain less
than the _a priori_ error threshold of 0.1.

:rainbow: Congrats! You have run your first MoFA model and verified the results against the fully-resolved simulation :checkered_flag:. If you have GLVis, you can view the simulation and closure problem results in the ```My_First_MoFA_Model/output``` directory. Feel free to look at
and modify the JSON-like input/configuration files ```MoFA_config.txt``` and ```porescale_config.txt```---as well as
```generate_input_file.h``` and ```generate_input_file_TUTORIAL.h```---to better understand how to apply the MoFA codes to your problems of interest.



## License :scroll:

MoFA is distributed under the terms of the BSD-3 license. All new contributions must be made under this license. See [LICENSE](LICENSE) and [NOTICE](NOTICE) for details.

LLNL Release Number:   LLNL-CODE-2006961

