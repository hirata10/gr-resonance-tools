
# gr-resonance-tools

Hello all! If you are reading this, you are here to use the gr-resonance-tools toolkit to help you along your own research projects. The README file below should give a good outline as to how you can use this code to replicate the results in our paper (arXiv number here) as well as how to modify it for your own needs. If you have any questions, you can reach the developers here:

hirata.10@osu.edu
k14masilv@gmail.com
h.g.blake-goszyk@vanderbilt.edu

Good luck and have fun!

________________________________________________________________________________________________________

Table of Contents:

Part 1: System Requirements/Installation
Part 2: Replicating (arxiv number)
Part 3: Physical Interpretations
Part 4: FAQ/Tips

________________________________________________________________________________________________________

PART 1:

SYSTEM REQUIREMENTS/INSTALLATION

***************************************************
Requirements:

This code requires that you have the lates versions of Python, Perl, Cython, and C installed. For those of you unfamiliar with it Cython, it can be described a wrapping language that allows you to use C functions in Python in one line. This menas that any of the C functions listed in the header file (CKerr.h) are fair game for Python users. To install Python and Cython, you can use pip:

pip install python
pip install cython

Nice that it rhymes, eh? If any issues arise with either Python or Cython, the documentation for both are accessible and easy to use. The C documentation is a bit denser, but is still nice to skim if needed:

http://docs.cython.org/en/latest/src/quickstart/install.html
https://www.python.org/downloads/
https://docs.python.org/3/
https://www.gnu.org/software/gnu-c-manual/gnu-c-manual.html
https://devdocs.io/c/


***************************************************
Running Cython code:

The way in which Cython is able to wrap C functions into Python is through three files: gr_wrapper.pxd, gr_wrapper.pyx, and setup.py. The first file acts as a sort of "header" file while the second file is actually doing the work of taking the inputs of the Python version of the function and putting them into the C version, and then back out again. The setup.py file takes all the C files with functions in the pxd file and "cythonizes" them. In order to initialize the setup.py file, type this into your terminal:

python setup.py build_ext --inplace

This will allow you to use all of the C functions in a python context. Note that the only function that is not wrapped is the rk4. If this statement is wrong, please let us know!


***************************************************
Running C code:

For those attempting to use this repository and being unfamiliar with C, do not fret! There are several C calculation routines that are already specified in the calling.c file. If you want to run any of them, type this command:





________________________________________________________________________________________________________

PART 2:

REPLICATING THE PAPER

***************************************************

This section provides a step-by-step guide on how to replicate the results from our paper, (). Note that this requires some knowledge of both Python and C, so be warned!

Step 1: To change the parameters of the EMRI systems you want to generate, go into pars.txt and change the values of the variables. Once you are satisfied, run the Perl script "pars.pl" (perl path/to/script/pars.pl). This should initialize a "globalpars.py" and "globalpars_c.h" file, which should allow you to use the designated variables in both coding languages.

Step 2: To initialize the Cython wrapper, type "python setup.py build_ext --inplace" into the terminal.

Step 3: The next step is to generate a bunch of potential EMRI systems with parameters as designated in the pars.txt file. Do this by running "action_angle_generator.py" into your terminal. You should get a text file with seven columns (system label, 6N action angle variables (J/mass)) named "action_angle_pairs.txt". If you wish to change the name of this file or it's location, you can do so by modifying line 34.

Step 4: Make the directory labeled "outputs_data" BEFORE compiling the following C code. Compile the C code "calling.c" using the following command "gcc -DIS_RK4_J_DOT calling.c resonance_find.c kerrtraj.c kerrmode.c kerrgwem.c J_dot.c J2J_dot.c -w -o J_evolve_single." This will generate an executable file with the name "J_evolve_single" that takes the following inputs: initial action variables (doubles), initial time (double), number of time steps (long), system label (long), system type (char; inner or outer). In order to run these executable files for each line in the files provided in Step 3, we use the python code "J_evolve_single.py." You can set the initial time, number of time steps, and data file name in lines 28-30 before executing the python code "J_evolve_single.py." The output will be the system label, the time, the action variables of that time, the corresponding frequencies (in the r, \theta, and \phi directions), and the time step. The output will go into a directory called "outputs_data" in txt format with file names "J_evolve_[system type]_[system label].txt"

Step 5: Upon completing all the rk4 runs, we find the resonances using "resfinder.py". This file does it by taking each pair of simulations (inner/outer) and checking them for potential resonances one-by-one. In order for this to happen, you must change lines 51/52 to read text from the directories you designated (outputs_data, Data/Inner_Body, etc). Once this is done, you can run resfinder.py from your terminal. You will get a text file names "potential_resonances.txt" that has a row for each potential resonance found within all of the simulations (resonance number, file number, resonance time, change in omega (inner/outer), inner body/SMBH mass ratio, outer body/SMBH mass ratio, gamma, n/k/m's, action variables for the inner and outer bodies (J/m_body), omegas for the inner and outer bodies, ancilliary data (inclination, periapse, apoapse)). If you wish to change the text file's name or location, you can do that in line 43.

Step 6: Make the directory labeled "Output_Delta_J" BEFORE compiling the following C code. Compile the C code "calling.c" using the following command "gcc calling.c -DIS_DELTA_J_SINGLE kerrtraj.c kerrmode.c kerrgwem.c resonance_find.c Gamma.c Delta_J.c J_dot.c -o Delta_J_single -w." This will generate an executable file with the name "Delta_J_single" that takes in the data from the resonance file and outputs the change in the action variables due to the resonance interaction. The python code "Delta_J_single_parallel.py" will run the executable file in parallel. On line 48 in the python file, input the name of the file with the resonance data (Step 3). Once the executable file is made, run the python code "Delta_J_single_parallel.py." The output file will have the resonance label (long; which resonance from the resonance file are we computing), the system label (long; in which system does this resonance occur), the integer mode labels for the inner and outer body (int; n_inner, k_inner, n_outer, k_outer, m_outer), the fundamental resonant angle (double; location on the torus), the total angular acceleration multiplied by the mass of each body (Gamma_outer - Gamma_inner), the frequencies of each body at that resonance resonance, the action variables for each body at that resonance, the final change in action variables for the inner body (double; J/(mass_inner)), and the ratio of the change in action variable to the action variable at resonance (this is for inner body only). The files will be in txt format in the directory "Output_Delta_J" with the naming scheme "Delta_J_log_{tot_chunk}_{chunk + 1}.txt," where tot_chunk is the total number of runs needed to cover the entire resonance data file (Step 3) and chunk+1 is the job in that run.

Step 7: [Figure out how we will classify resonances for astrophysical comparisons. Most likely compare the relative change in the action to the values of the action before the resonance. If this ratio is much smaller than 1, it is a perturbation in the same sense as our resonance calculation. If it is order unity/greater than one, the two bodies will undergo a transition to a non-linear interaction.]



________________________________________________________________________________________________________

PART 3:

PHYSICAL INTERPRETATIONS

***************************************************

This section will go through the physical interpretations of what the code is doing to get the results of (). More description can be found in the methods section.





________________________________________________________________________________________________________

PART 4:

FAQ/TIPS

***************************************************