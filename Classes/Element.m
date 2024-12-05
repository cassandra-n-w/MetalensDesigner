% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


classdef Element
    %Element: Contains S-parameters for multiple different sizes of a
    %single metal square element within metamaterial
    %   Detailed explanation goes here
    
    %Z_freespace = 376.73;
    
    properties
        Sparams
        sizes
        freqs
        Sparams_matrix
        impedance
        abcd_matrix
        abcd_grid
        s_grid
    end
    
    methods
        function obj = Element(Path)
            %Element Construct an instance of this class
            %   Using a path to touchstone files of various size
            %   parameters, create s parameters
            
            % list the files
            files = dir(fullfile(Path, '*.s4p'));
                       
            num_files = numel(files);
            
            %obj.Sparams = cell([num_files 1]);
            % loop through the files
            for i = 1:num_files
                %read in the s params
                %print(files(i))
                sparamtemp = sparameters(fullfile(files(i).folder, files(i).name));
                
                if (i > 1)
                    freqschange = obj.freqs ~= sparamtemp.Frequencies;
                    
                    if sum(freqschange) > 0
                        warning('Warning: different Element files are using different frequency vectors');
                    end    
                    
                    if obj.impedance ~= sparamtemp.Impedance
                        warning('Warning: different Element files are using different impedance values');
                    end
                end
                
                obj.freqs = sparamtemp.Frequencies;
                obj.impedance = sparamtemp.Impedance;
                
                obj.Sparams = [obj.Sparams sparamtemp];
                
                temp_mat = sparamtemp.Parameters([1 3], [1 3], :);
                
                obj.Sparams_matrix = cat(4, obj.Sparams_matrix, temp_mat);
                
                obj.abcd_matrix = cat(4, obj.abcd_matrix, s2abcd(temp_mat, obj.impedance));
                
                % read in the element size (note: this is based on the
                % filename and therefore maybe not be 100% robust)
                obj.sizes(i) = str2double(files(i).name(1:end-4));
            end
            
            
            % this sets up objects for supposedly-efficient future
            % interpolation. Hopefully they are as efficient as advertised.
            obj.s_grid = griddedInterpolant({1:2, 1:2, obj.freqs, obj.sizes}, obj.Sparams_matrix, 'linear', 'none');
            obj.abcd_grid = griddedInterpolant({1:2, 1:2, obj.freqs, obj.sizes}, obj.abcd_matrix, 'linear', 'none');
            
            return;
        end
        
        function out = Sparam_interp(obj, freq, size)
            % interpolate between the element's sparameters in both size
            % and frequency dimensions
            
            out = obj.s_grid({1:2, 1:2, freq, size});
            
        end
        
        function out = abcd_interp(obj, freq, size)
            % interpolate between the element's ABCD parameters in both size
            % and frequency dimensions
            size = squeeze(size);
            out = obj.abcd_grid({1:2, 1:2, freq, size});
        end
    end
end

