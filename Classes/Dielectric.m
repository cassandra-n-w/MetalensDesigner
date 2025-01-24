% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


classdef Dielectric < TwoPort_1Param
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Eps_r
        Eps_r_grid
        Tandelta
        unit
        Min
        Max
    end
    
    methods
        function obj = Dielectric(epsr, tand, freqs, minsize, maxsize)
            %Dielectric Create a new dielectric material using its relative
            %   permittivity and loss tangent.
            
            if isscalar(epsr)
                epsr = epsr * ones(size(freqs));
            end
            obj.Eps_r = epsr;
            obj.Eps_r_grid = griddedInterpolant(freqs, obj.Eps_r, 'linear','nearest');


            obj.Tandelta = tand;
            obj.freqs = freqs;
            obj.Min = minsize;
            obj.Max = maxsize;
            
            %Z0 = 376.730314; % ohms; impedance of free space
            obj.unit = logspace(log10(minsize), log10(maxsize), log10(maxsize/minsize) * 1000); % in mm; unit calculation length is 1um
            %c = 299792458000; % speed of light (mm/s);
            
            
             
            % calculate abcd matrices
            abcd_matrix = zeros(2,2,numel(freqs),numel(obj.unit));

            % convert things into matrices for vectorization
            epsr_m = reshape(epsr, 1, 1, []);
            freqs_m = reshape(freqs, 1, 1, []);
            unit_m = reshape(obj.unit, 1, 1, 1, []);

            s_matrix = zeros(2,2,length(freqs_m), length(unit_m));
            
            % note: this is probably reliant on good dielectric
            % approximation. can possibly be extended to include all
            % mu = 1.0 materials
            k = sqrt(epsr_m) .* 2.*pi.*freqs_m./Constants.c;
            beta = k;
            alpha = k .* tand / 2;
            gamma = alpha + 1.0i * beta;
            prop = exp(-gamma .* unit_m);
            
            s_matrix(2,1,:,:) = prop;
            s_matrix(1,2,:,:) = prop;
            %s_matrix(4,3,:) = prop;
            %s_matrix(3,4,:) = prop;
            
            %s_matrix = [0 prop(i) 0 0; prop(i) 0 0 0; 0 0 0 prop(i); 0 0 prop(i) 0];
            impedance = ones(size(unit_m)) .* Constants.Z0 ./ sqrt(epsr_m);

            s_matrix_flat = reshape(s_matrix, 2, 2, length(unit_m) * length(freqs_m));
            impedance_flat = reshape(impedance, 1, 1, length(unit_m) * length(freqs_m));
            
            abcd = s2abcd(s_matrix_flat, impedance_flat);

            % create a new S parameter matrix normalized to vacuum
            % characteristic impedance
            obj.impedance = impedance_flat(1);
            s_matrix_vacuum = s2s(s_matrix_flat, impedance_flat, obj.impedance);
            
            s_matrix_vacuum = reshape(s_matrix_vacuum, 2, 2, length(freqs_m), length(unit_m));
                
            abcd_matrix = reshape(abcd, 2, 2, length(freqs_m), length(unit_m));
            %abcd_matrix = [cosh(gamma .* unit_m), impedance.*sinh(gamma .* unit_m); sinh(gamma .* unit_m)./impedance, cosh(gamma .* unit_m)];

            obj.abcd_grid = griddedInterpolant({1:2, 1:2, obj.freqs, obj.unit}, abcd_matrix, 'linear', 'none');
            obj.s_grid = griddedInterpolant({1:2,1:2, obj.freqs, obj.unit}, s_matrix_vacuum, 'linear', 'none');
        end

        function out = abcd_interp(obj, freq, lengths)
            % interpolate between the element's ABCD parameters in both size
            % and frequency dimensions
            lengths = squeeze(lengths);
            %out = obj.abcd_grid({1:2, 1:2, freq, size});
            [s, impedance] = obj.calc_s(freq, lengths);

            s_matrix_flat = reshape(s, 2, 2, length(freq) * length(lengths));
            abcd = s2abcd(s_matrix_flat, impedance);
            out = reshape(abcd, 2, 2, length(freq), length(lengths));
        end

        function [s, impedance] = calc_s(obj, freqs, lengths)
            % calculate abcd matrices
            %abcd_matrix = zeros(2,2,numel(freqs),numel(obj.unit));

            % convert things into matrices for vectorization
            epsr = obj.Eps_r_grid(freqs);
            epsr_m = reshape(epsr, 1, 1, []);
            freqs_m = reshape(freqs, 1, 1, []);
            unit_m = reshape(lengths, 1, 1, 1, []);

            s_matrix = zeros(2,2,length(freqs_m), length(unit_m));
            
            % note: this is probably reliant on good dielectric
            % approximation. can possibly be extended to include all
            % mu = 1.0 materials
            k = sqrt(epsr_m) .* 2.*pi.*freqs_m./Constants.c;
            beta = k;
            alpha = k .* obj.Tandelta / 2;
            gamma = alpha + 1.0i * beta;
            prop = exp(-gamma .* unit_m);
            
            s_matrix(2,1,:,:) = prop;
            s_matrix(1,2,:,:) = prop;
            %s_matrix(4,3,:) = prop;
            %s_matrix(3,4,:) = prop;
            
            %s_matrix = [0 prop(i) 0 0; prop(i) 0 0 0; 0 0 0 prop(i); 0 0 prop(i) 0];
            impedance = ones(size(unit_m)) .* Constants.Z0 ./ sqrt(epsr_m);

            s_matrix_flat = reshape(s_matrix, 2, 2, length(unit_m) * length(freqs_m));
            impedance_flat = reshape(impedance, 1, 1, length(unit_m) * length(freqs_m));

            % create a new S parameter matrix normalized to vacuum
            % characteristic impedance
            impedance = impedance_flat(1);
            s = s2s(s_matrix_flat, impedance_flat, obj.impedance);
        end
    end
end

