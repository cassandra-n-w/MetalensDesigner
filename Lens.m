classdef Lens < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sim_frequency
        ideal_phase
        
        diameter
        struct_diam
        layer_thicknesses
        focal_length
        grid_dimension
        
        TL_prototype
        TL_array
        
        S_array
        S11_array
        dcode
    end
    
    methods
        function obj = Lens(diameter, thickness, focal, grid, proto, struct_diam)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.diameter = diameter;
            obj.struct_diam = struct_diam;
            obj.layer_thicknesses = thickness;
            obj.focal_length = focal;
            obj.grid_dimension = grid;
            obj.TL_prototype = proto;
            obj.dcode = 10;
            
            obj.TL_array = TL.empty(0);
            %obj.TL_array(obj.index_dimension, obj.index_dimension) = TL();
        end
        
        function gerberfy(obj)
            
            
            
            t = datetime;
            t.Format = 'yyyy-MM-dd''T''HHmm';
            str = string(t);
            g = obj.grid_dimension;
            mkdir('gerber' , str);
            
            %folder = strcat('gerber\', str, '\');
            
            prefix = '480GHz_layer';
            
            size_lens = size(obj.TL_array);
            extentx = size_lens(1);
            extenty = size_lens(2);
            
            xstart = -(extentx - 1)*g/2;
            ystart = -(extenty - 1)*g/2;
            
            for i = 1:(length(obj.layer_thicknesses) - 1)
                filename = fullfile('gerber', str, strcat(prefix, num2str(i), '.gbr'));
                gerberfile = fopen(filename, 'w');
                obj.dcode = 10;

                % set decimal format
                fprintf(gerberfile, '%%FSLAX36Y36*%%\n');

                % set to mm measurements
                fprintf(gerberfile, '%%MOMM*%%\n');

                obj.writeDrillMacro(gerberfile);

                obj.drillHoles(gerberfile, 24);

                obj.clearCircle(gerberfile,obj.diameter+2)

                obj.writeCrosshairs(gerberfile, 0.1)

                obj.clearCircle(gerberfile, obj.diameter+0.1)

                for j = 1:extentx
                    x = xstart+j*g;
                    if mod(j, 30) == 0
                        x
                    end
                    for k = 1:extenty
                         y = ystart+k*g;
                         
                         TL = obj.TL_array(j, k);
                         b = TL.sizes(i)*g;
                         if (x^2 + y^2 < (obj.diameter/2)^2)
                             obj.writeRect(x, y, b, b, gerberfile);
                         end
                    end
                end
                fprintf(gerberfile, 'M02*');
                fclose(gerberfile);
            end
            
            
        end

        function writeAlignmentMacro(obj, gerberfile)
            
        end

        function writeDrillMacro(obj, gerberfile)
            fprintf(gerberfile, '%%AMHOLEMARKER*\n');
            % outermost circle
            fprintf(gerberfile, '1,1,10.0,0,0*\n');
            % outermost circle
            fprintf(gerberfile, '1,0,9.5,0,0*\n');

            % second circle
            fprintf(gerberfile, '1,1,6.6,0,0*\n');
            % second circle
            fprintf(gerberfile, '1,0,6.3,0,0*\n');

            % third circle
            fprintf(gerberfile, '1,1,3.6,0,0*\n');
            % third circle
            fprintf(gerberfile, '1,0,3.4,0,0*\n');

            % inner circle
            fprintf(gerberfile, '1,1,0.2,0,0*%%\n');

        end

        function writeCrosshairs(obj, gerberfile, thickness)

            obj.writeRect(0,0,thickness, obj.diameter + 4, gerberfile)
            obj.writeRect(0,0,obj.diameter + 4, thickness, gerberfile)

        end

        function drillHoles(obj, gerberfile, numholes)
            r = (obj.struct_diam/2) - 7;
            obj.createDrillAperture(gerberfile);
            for i = 1:numholes
                
                ang = i/numholes * 2 * pi;
            
                x = sin(ang) * r;
                y = cos(ang) * r;

                obj.flashApp(x, y, gerberfile);             

            end
        
        end

        function clearCircle(obj, gerberfile, diam)
            
             fprintf(gerberfile, '%%LPC*%%\n');

             obj.writeCircle(0,0,diam,0,gerberfile)

             fprintf(gerberfile, '%%LPD*%%\n');

        end

        function createDrillAperture(obj, gerberfile)
            % create aperture
            dc = obj.dcode;
            fprintf(gerberfile, '%%ADD%dHOLEMARKER*%%\n', dc);
            fprintf(gerberfile, 'D%d*\n', dc);
            obj.dcode = dc + 1;
        end

        function writeCircle(obj, xpos, ypos, d_outer, d_inner, gerberfile)
            % create aperture
            dc = obj.dcode;
            if d_inner > 0
                fprintf(gerberfile, '%%ADD%dC,%#.6fX%#.6f*%%\n', dc, d_outer, d_inner);
            else
                fprintf(gerberfile, '%%ADD%dC,%#.6f*%%\n', dc, d_outer);
            end

            fprintf(gerberfile, 'D%d*\n', dc);
            obj.dcode = dc + 1;

            %flash aperture at each relevant location
            fprintf(gerberfile, 'X%dY%dD03*\n', round(xpos*1e6), round(ypos*1e6));
            
        end


        function flashApp(obj, xpos, ypos, gerberfile)
            %flash aperture at each relevant location
            fprintf(gerberfile, 'X%dY%dD03*\n', round(xpos*1e6), round(ypos*1e6));
        end
        
        function writeRect(obj, xpos, ypos, xsize, ysize, gerberfile)
            % create aperture
            dc = obj.dcode;
            fprintf(gerberfile, '%%ADD%dR,%#.6fX%#.6f*%%\n', dc, xsize, ysize);
            fprintf(gerberfile, 'D%d*\n', dc);
            obj.dcode = dc + 1;

            %flash aperture at each relevant location
            fprintf(gerberfile, 'X%dY%dD03*\n', round(xpos*1e6), round(ypos*1e6));
            
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

