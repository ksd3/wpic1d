function [density] = gather1d(xg, nx, rden, lden)
% This function is called by wpic1d.m to gather the particle data on a grid
% by linearly interpolating from the particle positions
%
% Inputs:
% xg is an array of particle positions normalized to grid values (from 1 to nx)
% nx is the number of elements in the grid


% Linear PIC interpolation (Macroparticles are tent functions)
% Add the component of the particle to the cell with an index less than x
%sprintf('Gather starting at')
%toc


rden(2,:) = floor(xg);    % The index below the grid point, xg
lden(2,:) = ceil(xg);

% How close xg is to floor(xg), weighted such that f(xg) = xg gives wxg = 1. This
% means if f(xg) = xg, the whole particle's density is at grid point f(xg).
rden(1,:) = 1 - (xg - rden(2,:)); 
lden(1,:) = 1 - (lden(2,:) - xg);

% rden(1,n) is the density weight of particle n to the grid point on the
% particle's right side, specified by rden(2,n). lden(1,n) is the density
% weight of particle n to the grid point on the particle's left side,
% specified by lden(2,n). These two lines basically wrap the grid.

rden(2,rden(2,:) == nx) = 0;
lden(2,lden(2,:) == nx) = 0;


%%%%%% If the particle is exactly on a grid point (usually only happens for
%%%%%% t = 0), then rden and lden are both claiming that particle with a
%%%%%% weight of 1. Thus the particle is double counted. Below we adjust
%%%%%% these weights to 0.5 so that the particle is no longer double counted.
rden(1, (rden(1,:) == 1)) = 0.5;
lden(1, (lden(1,:) == 1)) = 0.5;

% accumarray is now used on each density matrix to compile the total
% density at each grid point. The first line initializes the total density
% matrix with length nx. No clue why accumarray forces the matrices to be
% transposed. The rden(2,:) + 1 is used to shift grid point 0 to array
% index 1, since Matlab starts indices at 1.

density(1,nx+1) = 0;
density = accumarray((rden(2,:)+1)', rden(1,:)')';
%save('density_save.mat','density')
%save('lden_save.mat','lden')
density = density + accumarray((lden(2,:)+1)',lden(1,:)')';

% Wrap the grid around
%density(1)=density(1) + density(nx+1);

%density = density(1:nx);

%sprintf('Gather complete at')
%toc


end

