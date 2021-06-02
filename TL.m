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

