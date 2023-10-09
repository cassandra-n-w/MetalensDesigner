% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


classdef DiagonalHornClass
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Based on Johansson 1992, The diagonal horn as a sub-mmwave antenna
    
    properties
        a
        L
    end
    
    methods
        function obj = DiagonalHornClass(a, L)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.a = a;
            obj.L = L;
            
        end
        
        % we use eta and xi here because they align with the polarization
        % directions. eta is polarized axis and xi is cross-pol
        function out = EField(obj, eta, xi, f)
            [k, lambda] = Simulation.f_to_k_lambda(f);

            % transform to x/y coordinates
            x = (eta + xi)/sqrt(2);
            y = (eta - xi)/sqrt(2);

            lambda = reshape(lambda, 1,1,[]);

            mask = (abs(x) <= obj.a) .* (abs(y) <= obj.a);
                       
            k_del = (2*pi ./ lambda) .* ((2* (obj.a).^2 - x.^2 - y.^2))/ (2*obj.L);

            E_eta_mag = 1/(sqrt(2)) * (cos(pi*y / (2*obj.a)) + cos(pi*x/(2*obj.a)));

            out = mask .* E_eta_mag .* exp(1.0i * k_del);
        end
    end
end

