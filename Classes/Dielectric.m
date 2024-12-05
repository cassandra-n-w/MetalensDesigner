% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


classdef Dielectric
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Eps_r
        Tandelta
        Freqs
        abcd_grid
        unit
        Min
        Max
    end
    
    methods
        function obj = Dielectric(epsr, tand, freqs, minsize, maxsize)
            %Dielectric Create a new dielectric material using its relative
            %   permittivity and loss tangent.
            
            obj.Eps_r = epsr;
            obj.Tandelta = tand;
            obj.Freqs = freqs;
            obj.Min = minsize;
            obj.Max = maxsize;
            
            Z0 = 376.730314; % ohms; impedance of free space
            obj.unit = logspace(log10(minsize), log10(maxsize), log10(maxsize/minsize) * 60); % in mm; unit calculation length is 1um
            c = 299792458000; % speed of light (mm/s);
            
            
             
            % calculate abcd matrices
            abcd_matrix = zeros(2,2,numel(freqs),numel(obj.unit));
            
            for i=1:numel(freqs)
                % note: this is probably reliant on good dielectric
                % approximation. can possibly be extended to include all
                % mu = 1.0 materials
                k = sqrt(epsr(i)) .* 2*pi*freqs(i)/c;
                beta = k;
                alpha = k .* tand / 2;
                gamma = alpha + 1.0i * beta;
                prop = exp(-gamma * obj.unit);
                
                s_matrix = zeros(2,2,numel(obj.unit));
                s_matrix(2,1,:) = prop;
                s_matrix(1,2,:) = prop;
                %s_matrix(4,3,:) = prop;
                %s_matrix(3,4,:) = prop;
                
                %s_matrix = [0 prop(i) 0 0; prop(i) 0 0 0; 0 0 0 prop(i); 0 0 prop(i) 0];
                impedance = Z0 / sqrt(epsr(i));
                
                abcd = s2abcd(s_matrix, impedance);
                szabcd = size(abcd);
                
                for j = 1:szabcd(3)
                    abcd_matrix(:, :, i, j) = abcd(:,:,j);
                end
                
            end     
            
            obj.abcd_grid = griddedInterpolant({1:2, 1:2, obj.Freqs, obj.unit}, abcd_matrix, 'linear', 'none');
            
        end
        
        function out = abcd_interp(obj, freq, z)
            % interpolate between the dielectrics's ABCD parameters in both
            % length and frequency dimensions
            out = obj.abcd_grid({1:2, 1:2, freq, z});
            
%             out = zeros(4,4,numel(freq), numel(z));
%             for i = 1:numel(freq)
%                 for j = 1:numel(z)
%                     
%                     out(:,:,i,j) = (abcds(i))^(length(j)/obj.unit);
%                 end
%             end
            
        end
    end
end

