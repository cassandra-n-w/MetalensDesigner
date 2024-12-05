% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


function [A, kx, ky] = FourierEtoA(E, dx, dy)
%FourierEtoA Convert an E-field [ E(x,y) ] to a vector potential A(kx, ky)
%   Detailed explanation goes here
    % calculate vector potential
    % note the conjugate is taken here to account for phase convention
    A = fft2(conj(E));
    
    dims = flip(size(E) - [1, 1]);
    %dims = flip(size(E));
    
    % calculate kx and ky parts of the k vector
    %[kx, ky] = meshgrid(linspace(-dims(1)/2,dims(1)/2, dims(1)+1), linspace(-dims(2)/2,dims(2)/2, dims(2)+1));
    
    xn = 0:dims(1);
    yn = 0:dims(2);
    
    xn = xn - (xn > length(xn)/2) * length(xn);
    yn = yn - (yn > length(yn)/2) * length(yn);
    
    [kx,ky] = meshgrid(xn, yn);
    
    % this reorganizes the k vector components to align with how FFT is
    % calculated
    %kx = fftshift(fftshift(kx, 1), 2);
    %ky = fftshift(fftshift(ky, 1), 2);
    
    kx = 2*pi*kx / (dims(2)*dx);
    ky = 2*pi*ky / (dims(1)*dy);
end

