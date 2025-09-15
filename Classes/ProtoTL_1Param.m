classdef ProtoTL_1Param < handle
    %PROTOTL_GENERAL A more general version of ProtoTL
    %   Can handle multiple dielectrics, differently ordered dielectrics and elements,
    %   etc.
    
    properties
        TwoPortArray
    end
    
    methods
        function obj = ProtoTL_1Param(TwoPortArray)
            %PROTOTL_GENERAL Construct an instance of this class
            %   Detailed explanation goes here
            obj.TwoPortArray = TwoPortArray;
        end
        
        function out = SParam(obj, frequencies, sizes)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            Z0 = 376.730314; % ohms; impedance of free space

            prod = obj.TwoPortArray(1).abcd_interp(frequencies, sizes(1));

            for i = 2:length(sizes)
                prod = pagemtimes(prod, obj.TwoPortArray(i).abcd_interp(frequencies, sizes(i)));
            end

            out = abcd2s(prod, Z0);
        end

        function out = SParamTL(obj, frequencies, TL)
            out = obj.SParam(frequencies, TL.sizes);
        end
    end
end

