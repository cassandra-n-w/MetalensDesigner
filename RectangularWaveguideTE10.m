classdef RectangularWaveguideTE10
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Based on Pozar
    
    properties
        a
        b
    end
    
    methods
        function obj = RectangularWaveguideTE10(a,b)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.a = a;
            obj.b = b;
            
        end
        
        % we use eta and xi here because they align with the polarization
        % directions. eta is polarized axis and xi is cross-pol
        function out = EField(obj, x, y, f)

            mask = (abs(x) <= obj.a/2) .* (abs(y) <= obj.b/2);

            out = mask .* cos(pi*x/obj.a);
        end
    end
end

