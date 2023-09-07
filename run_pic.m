%%% MATLAB program for running wpic1d.m with an electron and ion fluid.

% run_pic.m can run any set of initial conditions for electrons and ions
% Built into this program are presets for acoustic waves and for a two
% stream instability. Not yet tested for more than 2 particle species

clear all;

%%%%%%  Format for initial conditions:
%%%%%%  particle0 = [np, q, m, kT, v0, den0, vel_perturb_mode, vel_perturb_mag]
%     np is the number of PIC macroparticles
%     q is the particle charge
%     m is the particle mass
%     kT is the particle distribution's initinal thermal temperature in energy units
%     v0 is the particle distribution's initial average velocity
%     den0 is the particle distribution's initial average density (the scaling between np and true density)
%     vel_perturb_mode is the oscillation mode of an initial perturbation.
%                       Don't set this to 0. vel_perturb_mode = 1 will disturb the initial
%                       with 1 sine wavelength, vel_perturb_mode = 2 will be 2 sine wavelengths...
%     vel_perturb_mag is the magnitude of the perturbation. Set this to 0 for no perturbation
%%%%%% particle1 = ... follows the same form as particle0 and allows the code to run with multiple species/populations.
%%%%%%  grid = [nx, dx, eps]
%     nx is the number of grid points.
%     dx is the spacing of the grid (must be less than the Debye length for numerical stability)
%            Total grid length is nx*dx.
%     eps is the dielectric constant of space
%%%%%%  time = [dt, n_skip, nmax]
%     dt is the time step of the simulation (less than dx/(v_thermal))
%     n_output is the frequency of output (n_output = 4 means plot every 4th timestep)
%     nmax is the number of time steps to calculate before stopping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%% 1D Initial Conditions %%%%%%%%%%%%%%%%%%%%%%

choice = menu('Initial Conditions','Two Stream','Ion Acoustic',...
    'Cold Plasma', 'Electron Langmuir (Plasma) Waves','Grid Heating', 'Landau Damping', 'Custom');

if choice == 1          % Two stream instability
    
    % Electron beam moving to the right
    % particle0 = [np, q, m, kT, v0, den0, vel_perturb_mode, vel_perturb_mag]
    particle0 = [1024*64, -1, 1, 1, 3, 1, 1, 0.1];
    % Electron beam moving to the left
    particle1 = [1024*64, -1, 1, 1, -3, 1, 1, 0.1];
    
    particle_params = [particle0; particle1];
    grid = [4096, 1, 1];          % [nx, dx, eps]
    time = [0.1, 10, 300];        % [dt, n_skip, nmax] 
    phaseflag = 1;
    
elseif choice == 2      % Ion acoustic waves
    
    % Ions,100 times heavier than electrons
    % particle0 = [np, q, m, kT, v0, den0, vel_perturb_mode, vel_perturb_mag]
    particle0 = [409600, 1, 100, 1, 0, 1, 4, 0.1];
    
    % Electrons
    particle1 = [409600, -1, 1, 4, 0, 1, 4, 0.1];
    
    particle_params = [particle0; particle1];
    
    grid = [128, 0.5, 1];          % [nx, dx, eps]
    time = [0.05, 10, 1000];        % [dt, n_skip, nmax]
    phaseflag = 0;
    
elseif choice == 3      % Cold plasma oscillations
    
    % Cold plasma, ions are treated as too boring to code
    % particle0 = [np, q, m, kT, v0, den0, vel_perturb_mode, vel_perturb_mag]
    particle0 = [1024*64, -1, 1, 0.01, 0, 1, 2, .5];
    
    particle_params = particle0;
    
    grid = [128, 1, 1];
    time = [0.1, 1, 1000];
    phaseflag = 0;
    
elseif choice == 4      % Langmuir modes. Same as plasma_waves in the IDL code
    
    % Only electrons present. Higher kT and different perturbation than
    % cold plasma.
    
    particle0 = [1024*64, -1, 1, 1, 0, 1, 5, 1]; 
    
    particle_params = particle0;
    
    grid = [128, 1, 1];
    time = [0.1, 1, 1000];
    phaseflag = 0;
    
elseif choice == 5
    % Use this area to code in a new set of initial conditions
    
    % particle0 = [np, q, m, kT, v0, den0, vel_perturb_mode, vel_perturb_mag]
    particle0 = [10000, -1, 1, 1, 0, 100, 0, 0];
      
    particle_params = particle0;
    
    grid = [512, 1, 1];          % [nx, dx, eps]
    time = [0.02, 1, 750];        % [dt, n_skip, nmax]
    phaseflag = 0;

elseif choice == 6
    % Use this area to code in a second new set of initial conditions
    
    % particle0 = [np, q, m, kT, v0, den0, vel_perturb_mode, vel_perturb_mag]
    particle0 = [10000, -1, 1, 0.4*0.4, 0, 1, 1, 0.1];
    
    
    particle_params = particle0;
    
    grid = [64, pi/32, 1];          % [nx, dx, eps]
    time = [0.2, 1, 75];        % [dt, n_skip, nmax]
    phaseflag = 0;
    
else 
    % In case you want to set up your own system without coding it in
    % particle0
    prompt = {['Number of Macroparticles, np'],['Particle charge, q'],['Particle mass, m'],['Initial thermal energy, kT'],...
        ['Average initial velocity, v0'],['Average initial density, den0'],...
        ['Velocity Perturbation Mode (1-4)'],['Velocity Pertubation Magnitude']};
    dlg_title = 'Particle 1';
    num_lines = 1;
    def = {'50000','1','100','1','0','1','4','0.05'};

    inputanswer = inputdlg(prompt,dlg_title,num_lines,def);

    particle0 = [str2num(inputanswer{1}),str2num(inputanswer{2}),str2num(inputanswer{3}),str2num(inputanswer{4})...
        str2num(inputanswer{5}),str2num(inputanswer{6}),str2num(inputanswer{7}),str2num(inputanswer{8})];
        
    % particle1
    prompt2 = {['Number of Macroparticles, np'],['Particle charge, q'],['Particle mass, m'],['Initial thermal energy, kT'],...
        ['Average initial velocity, v0'],['Average initial density, den0'],...
        ['Velocity Perturbation Mode (1-4)'],['Velocity Pertubation Magnitude']};
    dlg_title = 'Particle 2';
    num_lines = 1;
    def = {'50000','-1','1','4','0','1','4','0.05'};

    inputanswer2 = inputdlg(prompt2,dlg_title,num_lines,def);
    
    particle1 = [str2num(inputanswer2{1}),str2num(inputanswer2{2}),str2num(inputanswer2{3}),str2num(inputanswer2{4})...
        str2num(inputanswer2{5}),str2num(inputanswer2{6}),str2num(inputanswer2{7}),str2num(inputanswer2{8})];
    
    particle_params = [particle0; particle1];
    
    % Grid and time step
    % particle1
    prompt3 = {['Grid size, nx'],['Grid spacing, dx'],['Dielectric constant, eps'],['Time step, dt'],...
        ['Output skip, n_skip'],['Total time steps, nmax']};
    dlg_title = 'Grid, Time Step';
    num_lines = 1;
    def = {'128','0.5','1','0.05','20','10000'};

    inputanswer3 = inputdlg(prompt3,dlg_title,num_lines,def);
   
    grid = [str2num(inputanswer3{1}),str2num(inputanswer3{2}),str2num(inputanswer3{3})];
    time = [str2num(inputanswer3{4}),str2num(inputanswer3{5}),str2num(inputanswer3{6})];
    phaseflag = 1;
end
        
    
x = wpic1d(particle_params, grid, time, phaseflag);

