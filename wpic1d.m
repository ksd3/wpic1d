%%% MATLAB program to perform 1D plasma PIC simulations. Derived from
%%% Oppenheim's IDL program used for AS 727

% This is a function that simulates the plasma. This function is called by
% the following scripts, which mainly vary in terms of initial conditions:
%
% run_pic.m       Preset conditions for plasma waves and instabilities, 
%                   or customizable inputs for any intial conditions 
%                   instabily, 1 species cold plasma waves, 1 species 

% This function calls the following functions:
% 
% gather1d.m            Gathers particle data after time step
% poisson_invert.m      Inverts the Poisson equation for periodic boundary conditions



function [x] = wpic1d(particle_params, grid, time, phaseflag)

tic
% Input:
% particle_param must be an array of structures containing the following
% components describing each particle distribution's characteristics:
% 
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
%     
% grid is a structure containing:
%     nx is the number of grid points.
%     dx is the spacing of the grid (must be less than the Debye length for numerical stability)
%            Total grid length is nx*dx.
%     eps is the dielectric constant of space
%     
% time is a structure containing:
% 
%     dt is the time step of the simulation (less than dx/(v_thermal))
%     n_output is the frequency of output (n_output = 4 means plot every 4th timestep)
%     nmax is the number of time steps to calculate before stopping



%%%%% Initialize particle positions and velocities %%%%%

nspecies = size(particle_params,1);
Len = grid(1)*grid(2);              % Len = nx*dx  is the length of the grid
particles = cell(1,nspecies);

% Syntax for particles cell array:
% particles{1} contains a 2 by np array where particles{1}(1,m) is the
% position of the mth particle in fluid 1, and particles{1}(2,m) is the
% velocity of the mth particle in fluid 1
% particles{2} and so forth follow the same sytax
% This Cell Array method was chosen over a giant matrix so that np can be
% different for each fluid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% It would be nice to save the position and velocity of each
%%%%%% particle, at every time step, but each of the matrices below has about
%%%%%% 1 trillion entries in it. This overflowed 8 GB of memory, and forced
%%%%%% the computer to swap 1 GB of memory onto the hard drive. You can
%%%%%% probably save x and vx at all times for all particles simulated
%%%%%% if you have 16+ GB of memory (for 65,000 particles and 10,000 time
%%%%%% steps)

x_save = cell(1,nspecies);
vx_save = cell(1,nspecies);
%rden = cell(1,nspecies);
%lden = cell(1,nspecies);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Define intitial particle positions %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:nspecies
    
    % If you want a random distribution
    % particles{i}(1,:)= Len*random('Uniform',0,1,[1,particle_params(i,1)];

    % If you want a unifrom distribution
    particles{s} = [1:particle_params(s,1)]*Len/particle_params(s,1);
    
    % Gaussian velocity distribution.
    % Expression below is Random_number*sqrt(kT/m) + v0
    particles{s}(2,:) = random('Normal',0,1,[1,particle_params(s,1)])*...
        sqrt(particle_params(s,4)/particle_params(s,3)) + particle_params(s,5);
    
    % Add perturbations to the initial velocities
    % Expression below is v = v + vel_perturb_mag*Sin(2*pi*vel_perturb_mode*x/Len)
    %particles{s}(2,:) = particles{s}(2,:) + particle_params(s,8)*...
    %    sin(2*pi*particle_params(s,7)*particles{s}(1,:)/Len);
    
    % Add perturbations to initial particle positions
    gamma = 3;
    % Equation is sqrt(kT/m) where kT and m are summed over species (m = m1 + m2...)
    % The term vel_perturb_mag*0.7 is included for instances when kT is
    % near 0. 0.7 is the RMS value of sine.
    Cs = sqrt(gamma*sum(particle_params(:,4),1)/sum(particle_params(:,3)))...
        + particle_params(s,8)/sqrt(2);
%    Cs = 1;
    % Expression is x = x + (vel_perturb_mag/Cs)*Sin(2*pi*vel_perturb_mode*x/Len)
    particles{s}(1,:) = particles{s}(1,:) + (particle_params(s,8)/Cs)*...
        sin(2*pi*particle_params(s,7)*particles{s}(1,:)/Len);
    
    %%% Preallocate array for saving x and vx for each species
    %%% Array sizes are floor(nmax/n_output), where nmax is the number of
    %%% time steps taken, and n_output sets the frequency of data
    %%% recording. For simplicity, the code will always sample ~2500
    %%% particles. You will probably get errors if you use np < 2500. The
    %%% resulting arrays should have between 1 and 10 million entries,
    %%% depending on the length of the simulation.
    samp = 2500;
    ind = floor(particle_params(s,1)/samp);
    xsize = length(1:ind:length(particles{s}(2,:)));
    x_save{s} = zeros(floor(time(3)/time(2)),xsize);
    vx_save{s} = zeros(floor(time(3)/time(2)),xsize);
    
    %%% Preallocate arrays called by gather
    rden{s} = zeros(2,particle_params(s,1));
    lden{s} = rden{s};
    
end
save('initial_particles.mat','particles')
clear s
fprintf('Initial conditions determined \n')
toc

%%% Initial conditions have been set. Time to set the system in motion

% Preallocate matrices for saving data:
charge_den_save = zeros(floor(time(3)/time(2)),grid(1));
phi_save = zeros(floor(time(3)/time(2)),grid(1));
Ex_save = zeros(floor(time(3)/time(2)),(grid(1)+2));
den_save = zeros(floor(time(3)/time(2)),grid(1),nspecies);
KE = zeros(1,nspecies);
KE_save = zeros(floor(time(3)/time(2)), nspecies);      % Save kinetic energy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Begin time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:time(3)
    
    if mod(t,100) == 0  % Change the number in mod to change command line output cadence
        sprintf(['On the ' num2str(t) 'th time step at time \n'])
    end
    % First, calculate density and charge distribution:
    
    den = zeros(nspecies,grid(1));       % Mass density
    charge_den = zeros(1,grid(1));         % Charge density

    for s = 1:nspecies
        % Expression is density = gather(x/dx, nx)*den0*nx/(np)
        % gather1d.m is the function used to gather the particles onto the
        % grid and calculate density
        den(s,:) = gather1d(particles{s}(1,:)/grid(2), grid(1), rden{s}, lden{s})*...
            particle_params(s,6)*grid(1)/(particle_params(s,1));
        
        if mod(t,time(2)) == 0
            % Save the density
            den_save(t/time(2),:,s) = den(s,:);
        end        
        charge_den(1,:) = charge_den(1,:) + den(s,:)*particle_params(s,2);
    
        % While we're looping over particle species, let's calculate the
        % average kinetic energy of the species at the start of this
        % timestep. KE = 0.5*m*den0*<v^2>
        KE(1,s) = 0.5*particle_params(s,3)*particle_params(s,6)*...
            mean(particles{s}(2,:).^2);
        
    end
    
 
    
    % Calculate the electric potential from the charge density using a 
    % spectral technique in the poisson_invert.m function. Expression is
    % poisson_invert(charge_den*(-1/eps),dx, dy=1, dz=1)
    phi = poisson_invert(charge_den*(-1/grid(3)), 1, grid(2), 1, 1);
    
    
    % Calculate the electric field by finding the gradient of phi. grid(2) = dx
    Ex = [(circshift(phi,[0 -1]) - circshift(phi,[0 1]))/(-2*grid(2)), 0];
    
    % Enforce periodic boundary conditions
    Ex(length(Ex)) = Ex(1);     % For particles above grid point nx
    Ex = [Ex(length(Ex)-1) Ex];                % For particles below grid point 1
    
    %%% Save Ex, phi, KE, and charge_den
    if mod(t,time(2)) == 0
        Ex_save(t/time(2),:) = Ex;
        phi_save(t/time(2),:) = phi;
        charge_den_save(t/time(2),:) = charge_den(1,:);
        KE_save(t/time(2),:) = KE;
    end
    
        
    %%%% Great, now we can finally accelerate the particles
    clear s
    for s = 1:nspecies
        % vf = vo + at where a = E*q/m and the time step is dt = time(1)
        particles{s}(2,:) = particles{s}(2,:) + time(1)*...
            (particle_params(s,2)/particle_params(s,3))*...
            interp1((0:(length(Ex)-1)),Ex,((particles{s}(1,:)/grid(2))+1));       % Interpolate Ex
        
        if mod(t,time(2)) == 0
            ind = floor(particle_params(s,1)/samp);
            %%% Save velocities for every 2500 particles. The sampling
            %%% takes a particle every np/2500 columns in particles{s} 
            vx_save{s}((t/time(2)),:) = particles{s}(2, 1:ind:length(particles{s}(2,:)));
            %%%
        end
        
    end

    %%%%%save('particles.mat','particles')
    % Advance particle positions. This could probably be included in the
    % above loop
    for s = 1:nspecies
        % xf = xi + vt where v has just been updated. The expression is 
        % mod(xf + Len, Len) where Len is added to xf in order to avoid
        % taking the modulus of a negative number and getting an incorrect
        % answer. mod enforces periodic boundary conditions.
        particles{s}(1,:) = mod((particles{s}(1,:) + time(1)*particles{s}(2,:) + Len), Len);
        
        if mod(t,time(2)) == 0
            ind = floor(particle_params(s,1)/samp);
            %%% Save positions. Same scheme as the vx_save from above.
            x_save{s}((t/time(2)),:) = particles{s}(1, 1:ind:length(particles{s}(1,:)));
            %%%
        end
        
    end

    if mod(t, 200) ==0
        sprintf([num2str(t) ' time steps have completed \n'])
        toc
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Time loop finished %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprintf('Time loop complete, now saving \n')
toc

prefix = '';   % Default
% prefix = 'runnamehere_'; Save the run with a specific name, ie 2stream_position.txt
% prefix = 'Saved_Results/'; Save to a separate folder (folder must exist).
% prefix = 'Saved_Results/runnamehere/'; Save specific runs to separate folders

% Save particle inputs
dlmwrite([prefix 'particle_params.txt'], particle_params);
% Save grid, time, and phaseflag (phaseflag may not be needed?)
dlmwrite([prefix 'gridtime.txt'], [grid; time; phaseflag phaseflag phaseflag]);

% Save charge density
dlmwrite([prefix 'charge_den.txt'], charge_den_save);

% Save potential, phi
dlmwrite([prefix 'potential.txt'], phi_save);

% Save electric field, Ex
dlmwrite([prefix 'efield.txt'], Ex_save);

% x, vx and den are species dependent
for p = 1:nspecies
    %%% Save position for species p
    dlmwrite([prefix 'position' num2str(p) '.txt'], x_save{p});
    %%% Save velocity for species p
    dlmwrite([prefix 'velocity' num2str(p) '.txt'], vx_save{p});
    %%% Save density for species p
    dlmwrite([prefix 'density' num2str(p) '.txt'], den_save(:,:,p));
end
    
sprintf('Code done \n')
toc

x = 1;

end











