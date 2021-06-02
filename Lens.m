classdef Lens < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sim_frequency
        ideal_phase
        
        diameter
        layer_thicknesses
        focal_length
        grid_dimension
        
        TL_prototype
        TL_array
        
        S_array
    end
    
    methods
        function obj = Lens(diameter, thickness, focal, grid, proto)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.diameter = diameter;
            obj.layer_thicknesses = thickness;
            obj.focal_length = focal;
            obj.grid_dimension = grid;
            obj.TL_prototype = proto;
            
            obj.TL_array = TL.empty(0);
            %obj.TL_array(obj.index_dimension, obj.index_dimension) = TL();
        end
        
        function out = CalcSParam(obj,frequency)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            sz_array = size(obj.TL_array);
            
            obj.S_array = zeros(sz_array(1), sz_array(2), numel(frequency));
            obj.sim_frequency = frequency;
            
            for i = 1:sz_array(1)
                for j = 1:sz_array(2)
                    thisTL = obj.TL_array(i,j);
                    
                    S = obj.TL_prototype.SParam(frequency, thisTL.sizes, thisTL.thicknesses);
                    
                    obj.S_array(i,j,:) = S(2,1,:);
                end
            end
            
            out = obj.S_array;
        end
        
        
        function out = PadSParams(obj, pad, frequency)
            %PadSParams This method inserts the s-parameters of this lens
            %into a larger set of s-parameters contained in pad
            %Note: calculate only at a single frequency. if multiple
            %frequencies are given, the first will be used.
            
            idx = find(obj.sim_frequency == frequency(1), 1);
            
            if numel(idx) > 0
                s_param = obj.S_array(:,:,idx);
            else
                obj.CalcSParam(frequency);
                s_param = obj.S_array(:,:,1);
            end
            
            szpad = size(pad);
            szsim = size(obj.S_array);
%             
%             maxlim = ceil(szpad/2 + szsim/2);
%             minlim = ceil(1 + szpad/2 - szsim/2);
%             
%             pad(minlim(1):maxlim(1), minlim(2):maxlim(2)) = s_param;
            
            padamount = ceil((szpad - szsim)/2);

            temp = padarray(s_param, padamount, 0);
            sztemp = size(temp);
            
            szerror = sztemp - szpad;
            
            out = temp(1:end - szerror(1), 1:end - szerror(2));
        end
    end
end

