% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


function [EField] = HornField(a1, b1, p, a, b, k, x, y)
%UNTITLED Summary of this function goes here
%   Input parameters align with what is shown in Fig 13.16 of 
%   Antenna Theory: Analaysis and Design by Balanis, 4th Ed.
    
    % for a pyramidal horn, we assume p_e = p_h = p
    % we then use properties of similar triangles to find the
    % center-lengths
    rho1 = p * b1 / (b1 - b);
    rho2 = p * a1 / (a1 - a);
    
    mask = (abs(x) <= a1/2) .* (abs(y) <= b1/2);

    EField = mask .* cos(pi*x/a1) .* exp(-1.0i * k * (x.^2/rho2 + y.^2/rho1) / 2);
    
end

