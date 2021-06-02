classdef ProtoTL < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Dielectric
        Element
    end
    
    methods
        function obj = ProtoTL(element, dielectric)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Element = element;
            obj.Dielectric = dielectric;
        end
        
        function out = SParam(obj, frequencies, sizes, thicknesses)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            Z0 = 376.730314; % ohms; impedance of free space
            
            el_abcd = obj.Element.abcd_interp(frequencies, sizes);
            di_abcd = obj.Dielectric.abcd_interp(frequencies, thicknesses);
            
            szel = size(el_abcd);
            szdi = size(di_abcd);
            
            if szdi(4) ~= (szel(4) + 1)
                warning('Warning: expected there to be one more dielectic layer than element layers');
            end    
            
            % calculate the abcd parameters of the full TL by multiplying
            % together the abcd parameters of its individual elements.
            prod = di_abcd(:,:,:,1);
            
            for i = 1:szel(4)
                for j = 1:length(frequencies)
                    prod(:,:,j) = prod(:,:,j) * el_abcd(:,:,j,i) * di_abcd(:,:,j,i+1);
                end
            end
            
            out = abcd2s(prod, Z0);
        end
    end
end

