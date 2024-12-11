classdef TL_interface
    %TL_INTERFACE Summary of this class goes here
    %   Detailed explanation goes here

    methods (Abstract)
        SParam(obj, frequencies)
        SParamModify(obj, frequencies, thickness_mult)
    end
end

