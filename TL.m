% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


classdef TL
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        MyProtoTL
        sizes
        thicknesses
    end
    
    methods
        function obj = TL(element_sizes, dielectric_thicknesses, proto)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0
            	obj.thicknesses = dielectric_thicknesses;
                obj.sizes = element_sizes;
                obj.MyProtoTL = proto;
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

