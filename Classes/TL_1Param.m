classdef TL_1Param
    %TL_1PARAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        my_proto_tl
        sizes
    end
    
    methods
        function obj = TL_1Param(sizes, proto)
            %TL_1PARAM Construct an instance of this class
            %   Detailed explanation goes here
            obj.sizes = sizes;
            obj.my_proto_tl = proto;
        end
        
        function out = SParam(obj, frequencies)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            out = obj.my_proto_tl.SParamTL(frequencies, obj);
        end

        function out = SParamModify(obj, frequencies, thickness_mult)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            out = obj.SParam(frequencies);
            warning("SParamModify not yet implemented for TL_1Param. Defaulting to standard SParam behaviour");
            
        end
    end
end

