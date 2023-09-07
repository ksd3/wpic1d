%%% plot_1dpic.m is used to plot the data produced by 1D runs of run_pic.m
clear all
%close all

prefix = '';    % Default
% prefix = 'Saved_Results/runnamehere/'; % If outputpath was changed in wpic1d.m

% Plot window locations:
monitors = menu('Choose screen size', '23 inch monitor', 'Small laptop');

% Determine animation speed or turn animation off.
speedflag = menu('Choose Animation Speed','1x', '5x', '10x', '50x');

if speedflag == 1 % 1x speed
    speed = 1;
elseif speedflag == 2 % 2x speed
    speed = 5;
elseif speedflag == 3 % 3x speed
    speed = 10;
else
    speed = 50; % 50x or custom speed
end

%%% Start by loading the original input parameters of run_pic
tic
particle_params = csvread([prefix 'particle_params.txt']);
gridtime = csvread([prefix 'gridtime.txt']);

nspecies = size(particle_params,1);     % Check index

%%%%% Modify this if you need to plot 6 or more particle species for some reason
if nspecies > 5
    nspecies = 5;
    sprintf('The plot functions have not been set up for more than 5 particle species \n');
    springf('Only showing plots for first 5 species \n');
end

%%% Load the quantities dependent only on grid size

charge_den = csvread([prefix 'charge_den.txt']);
efield = csvread([prefix 'efield.txt']);
potential = csvread([prefix 'potential.txt']);

%%% Load the quantities dependent on particle species

x = cell(1, nspecies);
vx = cell(1, nspecies);
for s = 1:nspecies

    density(:,:,s) = csvread([prefix 'density' num2str(s) '.txt']);
    x{s} = csvread([prefix 'position' num2str(s) '.txt']);
    vx{s} = csvread([prefix 'velocity' num2str(s) '.txt']);
    
    vxmin(s) = min(min(vx{s}));     % To set vx plot range
    vxmax(s) = max(max(vx{s}));
    density(:,:,s) = density(:,:,s)/particle_params(s,6);
end
    vxmin = min(vxmin);
    vxmax = max(vxmax);

clear s
sprintf('Data reloaded \n')
toc


        
% Determine Electric field energy: E_en = E^2*eps/(2*Len)
E_energy = sum(efield(:,(2:(gridtime(1,1)-1))).^2,2)'*gridtime(1,3)/(2*gridtime(1,1)*gridtime(1,2));
        
% Determine particle energy:
for i = 1:nspecies
    % Kinetic energy of particle is 0.5*m*den0*<v^2>
    KE(i,:) = 0.5*particle_params(i,3)*particle_params(i,6)*...
        mean(vx{i}(:,:).^2,2);
    maxKE(i) = max(KE(i,:));
end
maxKE = max(maxKE);     % To fix the axis on the energy plot
maxE = max(max(E_energy),maxKE);

% Vector of the times n_output plots at
tplot = (0:gridtime(2,2):(gridtime(2,3)-1))*gridtime(2,1)';
% Vector of the spatial grid coordinates
xplot = (0:(gridtime(1,1)-1))*gridtime(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Plot stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot 1 is Energy vs Time (electric and kinetic energy
%%% Plot 2 is Electric potential at each grid point
%%% Plots 3-n are the Density at each grid point
%%% Plots n+ are phase space, Vx versus X for each particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set figure locations
fignum = 2 + 2*nspecies;
fig = zeros(fignum,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Window Placement %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basepos, offset, dx, and dy use numbers based on a 23 in widescreen monitor.
% Tweak these values if your plots show up stacked on each other.

%%% Left: Position (in pixels) of left figure edge from left edge of monitor
%%% Bottom: Position of bottom figure edge from bottom edge of monitor
%%% Width: How wide eachthe figure is in pixels
%%% Height: How tall each figure is in pixels

if monitors == 1 %%% Set for 23 inch monitor
    basepos = [20, 740, 500, 280];      % [Left, Top, Width, Height]
    dx = 500;       % Shift successive images to the right by dx pixels
    dy = 370;       % Shift successive images down by dy pixels
    % Each row of offset corresponds to the nth plot being modified by
    % [Right shift, Left shift, Width change, Height change]
    offset = [0 0 0 0; dx 0 0 0; 2*dx 0 0 0; 3*dx 0 0 0;...
        0 -dy 0 0; dx -dy 0 0; 2*dx -dy 0 0; 3*dx -dy 0 0;...
        0 -2*dy 0 0; dx -2*dy 0 0; 2*dx -2*dy 0 0; 3*dx -2*dy 0 0];
else %%% This should work for a laptop display, but tops out at 6 windows
    basepos = [20, 370, 400, 280];      % [Left, Top, Width, Height]
    dx = 400;       % Shift successive images to the right by dx pixels
    dy = 350;       % Shift successive images down by dy pixels
    % Each row of offset corresponds to the nth plot being modified by
    % [Right shift, Left shift, Width change, Height change]
    offset = [0 0 0 0; dx 0 0 0; 2*dx 0 0 0;...
        0 -dy 0 0; dx -dy 0 0; 2*dx -dy 0 0];
end
    
%%% Error catch if the code wants to plot more figures than offset has positions
if size(offset,2) < fignum
        offset(fignum, 4) = 0; % This will stack extra windows on top of Figure 1. Not ideal
end

for m = 1:fignum
    fig(m) = figure('position', basepos + offset(m,:));
end
    fig(fignum+1) = figure('position', basepos + offset(end,:));

% Plot the power spectrum
    Y = 10*log10((abs(fft2(density))));
    Y(1,1) = 0;
    Y = 0.5*Y(1:end/2, 2:end/2, 1) + 0.5*[Y(1,2:end/2,1); Y(end:-1:(end/2+2), 2:end/2,1)];
    Yw = [0:size(density,1)/2-1]*2*pi/(gridtime(2,1)*gridtime(2,3));
    Yk = [1:size(density,2)/2-1]*2*pi/(gridtime(1,1)*gridtime(1,2));
    figure(fignum+1); clf; hold on
    imagesc(Y(:,:,1));
    imagesc('XData',Yk,'YData',Yw,'CData',Y);
    colormap jet
    colorbar
    title('Power Spectrum (dB)');
    xlabel('Wavenumber, k'); ylabel('Frequency,w');
    axis([0 Yw(end) 0 Yk(end)]);    
    

for t = gridtime(2,2):gridtime(2,2):gridtime(2,3)
    
    % Plot only on the cadence of n_output.
    num = t/gridtime(2,2);
      
    % Only plot according to speedflag
    if mod(num, speed) == 0
        
        figure(1); clf;
        plot(tplot(1:num),E_energy(1:num),':',tplot(1:num),KE(:,1:num));
        title('Energy vs Time');
        xlabel('Time');
        ylabel('Energy');
        legend('E^2 density', 'Particle 1 KE');
        axis([0 max(tplot) 0 maxE*1.05]);
            

        figure(2); clf;
        plot(xplot,potential(num,:));
        %plot(xplot,efield(num,2:(size(efield,2)-1)));
        title('Potential at each grid point');
        xlabel('Position, x');
        ylabel('Potential');
        axis([0 max(xplot) min(0.8*min(min(potential)), 1.1*min(min(potential))) 1.2*max(max(potential))]);
        

        for s = 1:nspecies

            %%% Note: Density and velocity are plotted on the same scale for
            %%% each particle in order aid comparison. If one particle has a
            %%% density or velocity that is orders of magnitude higher, than
            %%% issues may arise and you should change the axis commands below
            figure(2*s+1); clf;
            plot(xplot,(density(num,:,s)-1))
            title(['Spatial Density of Particle' num2str(s)]);
            xlabel('Position, x');
            ylabel('Density - 1');
            dmax = max(max(max(density(:,:,:)))) - 1;
            dmin = min(min(min(density(:,:,:)))) - 1;
            axis([0 max(xplot) 1.1*dmin 1.1*dmax]);
            

            figure(2*s+2); clf; hold on
            scatter(x{s}(num,:),vx{s}(num,:),9,'filled','k')
            % The following highlights a few particles in phase space
            %scatter(x{s}(num,[700 1400 2100]),vx{s}(num,[700 1400 2100]),25,'filled','r')
            hold off
            title(['Vx versus X for Particle' num2str(s)]);
            xlabel('Position, x');
            ylabel('Velocity, vx');
            axis([0 max(xplot) vxmin vxmax]);
            
        end

        
        drawnow

    end

end
