READ ME for wpic1d MATLAB code used in AS 727 Plasmas, written by William Longley and extended by Kshitij Duraphe

This code was written for MATLAB R2016a, and has been tested on versions R2011a and 2018b

Components of the code:
run_pic.m   This is the script used to run the code. It sets the initial
	    conditions and calls wpic1d.m. The code outputs a set of .txt
plot_1dpic.m  This is the script used to plot the outputs of a wpic1d run.
	      Reads the .txt outputs of a run and plots the particle positions,
	      electric field, and phase space locations as a function of time.
wpic1d.m    This is the main code. It is not necessary to open or modify this file.
gather1d.m  This is a function called by wpic1d.m to gather the density onto the grid
poisson_invert.m  This is a function called by wpic1d.m that uses the gathered
		  density to solve Poisson's equation for the electric field.
		  There is no reason to modify gather1d.m or poisson_invert.m

INSTALLATION:
Make sure all 5 of the above files are in the same directory.

HOW TO USE:
Open run_pic.m. Click the run button. Select the desired initial conditions. Wait.
Open plot_1dpic.m. Click the run button. Follow the prompts.

INPUT CONFIGURATION: The run_pic.m file was deliberately written as a separate script
from wpic1d.m to allow initial conditions to be easily changed. Some default
configurations are present, and the user should look at what the initial values are
before using the code for the AS 727 assignments. The default values may be incorrect
for the desired simulation result or change year to year. The beginning of the script 
defines what all the variables are.

___________________________________
ADVANCED OPTIONS:
Some configurations to the code are listed below. Line references are approximate.

Window placement in plot_1dpic.m:
The first prompt in plot_1dpic.m asks for your monitor size. The two default options
are for a 23 inch monitor (looks great) and a ~14 inch laptop screen (may not look
great). If neither option is appealing, the plot window sizes and locations can be
changed in plot_1dpic.m at lines ~100 through ~125. Ideally you should only need to 
change 'basepos' (the top left corner of the first window) and 'dx' and 'dy' (width
and  height of each plot window). 'offset' is the array of window postions for the 
plots, and assumes 4 windows in each row. This may need to be changed to have only 2 
or 3 windows per row, which is done by removing the last entry in each row of the 
array (for example, the '3*dx 0 0 0' entry in the first row).

Window flashing in plot_1dpic.n: The annoying flash of the plot windows is an OS problem.

Output configuring: 
By default the run outputs are dumped in the same directory that run_1dpic.m is in. 
The major downside of this is that by default, new runs overwrite the output of old 
runs. This can be changed in wpic1d.m at line ~255. The variable 'prefix' can be 
changed so that the run is named, or the run is saved to a subdirectory (must exist 
before running the code). 3 examples are given as comments. plot_1dpic.m will also 
need to be modified at line ~5 to change 'prefix' to whateverit was changed to in 
wpic1d.m

Command line spam: 
By default the code prints to the command line every 100 timesteps. If this is too 
much or too little, at line ~150 is a statement 'if mod(t,100) == 0'. Changing the 
100 to another number will change how often the code prints the current time step. 
'mod(t,1)' will print every time step

Animation Speed: 
By default, plot_1dpic.m has 4 preset animation speeds. The 1x option plots every 
timestep of the simulation (slow), the 5x option plots every 5th timestep and so 
forth. This can be changed in lines ~10 through ~20 by changing the 'speed' variable.
My suggestion is to change the value of the 50x speed.

Initial particles:
The initial particle positions are output in the initial_particles.mat file. This is
generally not needed, but it does allow you to check the initial conditions to see
if they're as expected.

Units:
The code was written to use natural units of everything = 1. This can be changed some
by putting in real masses, charges, etc in the run_pic.m file, but it may not come out
consistent in the end. The natural unit choice also means the output plots of
potential, energy, etc are in weird units as well.

Memory issues: 
This code was developed and run on a machine with 8 GB of memory. I do not guarantee 
the code will be able to run without issues on machines with less memory. If you have
less than 8 GB of memory and the code appears to be taking forever to run, you're 
probably caching to virtual memory on your hard drive at every time step which is 
extremely slow. The easiest way to handle this is to reduce the number of particles
in the simulation, which is the variable np in run_pic.m. You can also reduce the
number of particles which are actually saved to the output files if plot_1dpic.n is slow.
This is done at line ~120 in wpic1d.m, by changing 'samp' to a smaller number. Other ways
of dealing with the memory issue and speed are reducing overall run size 
(np, nx, or nmax in run_pic.m), or modernizing your computer.
