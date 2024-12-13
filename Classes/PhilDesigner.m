classdef PhilDesigner
    %PHILDESIGNER Design a lens using Phil's equations
    %   Detailed explanation goes here
    
    properties
        lens_model
        sim
    end
    
    methods
        function obj = PhilDesigner(lens_model, sim)
            %PHILDESIGNER Construct an instance of this class
            %   Detailed explanation goes here
            obj.lens_model = lens_model;
            obj.sim = sim;
        end
        
        function out = DesignDielectricGrin(obj)
            %METHOD1 Design lens using GRIN equation
            %   Assumes ideal varying dielectric (not realistic)
            
            % dist_field = sqrt(obj.lens_model.x.^2 + obj.lens_model.y.^2);
            % 
            % for i = 1:length(obj.lens_model.xvec)
            %     disp(i);
            %     for j = 1:length(obj.lens_model.yvec)
            %         dist = min(dist_field(i,j), obj.lens_model.diameter/2);
            %         obj.lens_model.TL_array(i,j) = obj.DesignDielectricGrinTL(dist);
            %     end
            % end
            
            dr = obj.lens_model.grid_dimension / 4;
            radii = 0:dr:(obj.lens_model.diameter/2);

            default_tl = TL_1Param(0,0);
            TL_collection = createArray([length(radii),1], FillValue=default_tl);

            for i=1:length(radii)
                if (mod(i,10) == 0)
                    disp(i);
                end
                TL_collection(i) = obj.DesignDielectricGrinTL(radii(i));
            end

            obj.lens_model.SetRadialDesign(TL_collection, radii);

        end

        function out = DesignDielectricGrinTL(obj, dist)
            %DesignDielectricGrinTL Given a dist from the lens center,
            %design a single dielectric grin TL
            t_tot = obj.lens_model.layer_thicknesses;
            
            n0 = 7; % maximum index of refraction for the lens

            lamb_x = Constants.c/obj.sim.designfrequency; % X-band (9 GHz) wavelength in mm
            
            f = obj.lens_model.focal_length; % focal length in millimeters
            
            t0 = t_tot - 2.*lamb_x/(4. * sqrt(n0)); % Thickness of the high index region in the center of the lens
            
            r = dist; %radius in millimeters
            
            b = (lamb_x^2/(4.*t_tot^2)) + 2*(n0*t0 - r^2/(2*f))/t_tot; %coefficient
            
            c = ((n0*t0 - r^2/(2*f))/t_tot)^2; % coefficient
            
            n_of_r = (b + sqrt(b^2 - 4.*c))/2.; % value of high index material as a function of radius
            
            t_of_r = (n0*t0 - r^2/(2*f))/n_of_r; % thickness of high index material as a function of radius
            
            phase = 180.*(t_tot - t_of_r)*sqrt(n_of_r)/lamb_x; % phase of AR coating (should be 90 degrees)
            
            ar_thick = (t_tot - t_of_r)/2; % thickness of AR coating as a function of radius (index is the square root of the high index part)
            
            n_thick = t_of_r; % thickness of the high index part again

            inner_dielectric = Dielectric(n_of_r.^2, 0, obj.sim.frequency, t_of_r./1.2, t_of_r.*1.2);

            % note we use n_inner as the AR dielectric because epsilson_ar is the
            % the square of n_ar, but n_ar is the square root of n_inner
            ar_dielectric = Dielectric(n_of_r, 0, obj.sim.frequency, ar_thick./1.2, ar_thick.*1.2);
    
            proto = ProtoTL_1Param([ar_dielectric, inner_dielectric, ar_dielectric]);
            TL = TL_1Param([ar_thick, n_thick, ar_thick], proto);

            out = TL;
        end
    end
end

