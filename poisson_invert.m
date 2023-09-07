function phi = poisson_invert(b, dim, dx, dy, dz)
% This is a function to calculate the potential phi in Poisson's equation,
% del^2(phi) = b, where b is a given quantity. Assumes periodic boundary
% conditions. Returns a complex array.

% Function must be called with a dy and dz


Nx = size(b,2);
hnx = floor(Nx/2);

cx = zeros(1,Nx);
cx(1:(hnx+1)) = [0:hnx]*2*pi/(dx*Nx);
cx((hnx+2):Nx) = -(hnx - 1 - [0:(hnx-2)])*2*pi/(dx*Nx);

P = -1./cx.^2;
P(1) = 1;

    
btrans = fft2(b);       % Fourier transform
btrans = btrans/(Nx); %%%%%%%%      % IDL's fft automatically normalizes by 1/N
phitrans = btrans.*P;   % phitrans is the Fourier transform of phi
phitrans(1) = 0;        % Set the DC term to 0
phi = ifft2(phitrans);  % We came here for phi
phi = real(phi)*Nx; %%%%%%%%        % Still not on Reimmann sphere


end
