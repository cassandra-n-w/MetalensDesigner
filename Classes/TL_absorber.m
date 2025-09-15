classdef TL_absorber
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        function obj = TL_absorber()
            
        end
        
        function out = SParam(obj, frequencies)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            out = zeros([2,2,length(frequencies)]);
        end

        function out = SParamModify(obj, frequencies, thickness_mult)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            out = obj.SParam(frequencies);  
        end
    end
end

