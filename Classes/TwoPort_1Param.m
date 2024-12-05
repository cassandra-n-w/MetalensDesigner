classdef TwoPort_1Param
    %TWOPORT_1PARAM A class containing one two-port network varied by a
    %single parameter. 
    %   (E.G. varying length for a dielectric or varying
    %square size for an element). Not intended to be instantiated as an
    %object, only to act as a superclass
    
    properties
        freqs
        abcd_grid
    end
    
    methods
        function out = abcd_interp(obj, freq, size)
            % interpolate between the element's ABCD parameters in both size
            % and frequency dimensions
            size = squeeze(size);
            out = obj.abcd_grid({1:2, 1:2, freq, size});
        end
    end
end

