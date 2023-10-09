% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


function [E] = FourierAtoE(A, z, kx, ky, k0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
%    zscale(1,1,:) = z;

%     kz = repmat(sqrt(k0^2 - kx.^2 - ky.^2), 1, 1, length(z));
%     
%     dims = size(A);
%     zscale = repmat(zscale, dims(1), dims(2));
%     A = repmat(A, 1, 1, length(z));
%     
%     Az = A .* exp(kz .* zscale * 1.0i);
    
    kz = sqrt(k0^2 - kx.^2 - ky.^2);
    zscale = permute(z, [3 1 2]);
    Az = (A .* exp(kz .* zscale * 1.0i));    

    % note the conjugate is taken to take into account phase convention
    E = conj(ifft2(Az));
end

