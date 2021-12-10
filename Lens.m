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
        S11_array
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
            obj.S11_array = zeros(sz_array(1), sz_array(2), numel(frequency));
            obj.sim_frequency = frequency;
            
            proto = obj.TL_prototype;
            Sout = zeros( [sz_array(1) sz_array(2) length(frequency)]);
            S11out = zeros( [sz_array(1) sz_array(2) length(frequency)]);
            
            imax = sz_array(1);
            jmax = sz_array(2);
            for i = 1:imax
                i
                for j = 1:jmax
                    thisTL = obj.TL_array(i,j);
                    
                    S = proto.SParam(frequency, thisTL.sizes, thisTL.thicknesses);
                    
                    Sout(i,j,:) = S(2,1,:);
                    S11out(i,j,:) = S(1,1,:);
                end
            end
            obj.S_array = Sout;
            obj.S11_array = S11out;
            out = obj.S_array;
        end
        
        function out = PadMultiFreq(obj, pad, frequencies)
            sz = size(pad);
            out = zeros([sz(1), sz(2), length(frequencies)]); 
            
            for i = 1:length(frequencies)
                out(:,:,i) = obj.PadSParams(pad, frequencies(i));
            end
        end
        
        function out = PadGeneral(obj, input, pad, frequencies)
            sz = size(pad);
            out = zeros([sz(1), sz(2), length(frequencies)]); 
            
            for i = 1:length(frequencies)
%                 idx = find(obj.sim_frequency == frequencies(i), 1);
%             
%                 if numel(idx) > 0
%                     s_param = input(:,:,idx);
%                 else
%                     obj.CalcSParam(frequency);
%                     s_param = input(:,:,1);
%                 end
                s_param = input(:,:,i);

                szpad = size(pad);
                szsim = size(s_param);
    %             
    %             maxlim = ceil(szpad/2 + szsim/2);
    %             minlim = ceil(1 + szpad/2 - szsim/2);
    %             
    %             pad(minlim(1):maxlim(1), minlim(2):maxlim(2)) = s_param;

                padamount = ceil((szpad - szsim)/2);

                temp = padarray(s_param, padamount, 0);
                sztemp = size(temp);

                szerror = sztemp - szpad;

                out(:,:,i) = temp(1:end - szerror(1), 1:end - szerror(2));
            end
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
            szsim = size(s_param);
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

