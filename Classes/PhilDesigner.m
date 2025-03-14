classdef PhilDesigner
    %PHILDESIGNER Design a lens using Phil's equations
    %   Detailed explanation goes here
    
    properties
        lens_model
        sim
        total_thickness
    end
    
    methods
        function obj = PhilDesigner(lens_model, sim)
            %PHILDESIGNER Construct an instance of this class
            %   Detailed explanation goes here
            obj.lens_model = lens_model;
            obj.sim = sim;
            obj.total_thickness = obj.lens_model.layer_thicknesses;
        end
        
        function DesignDielectricGrin(obj)
            obj.DesignDielectricGrinGeneral(0);
        end

        function DesignDielectricGrinGeneral(obj, layer_thickness)
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
            pathlength = obj.FindIdealPathLength(radii);

            default_tl = TL_1Param(0,0);
            TL_collection = createArray([length(radii),1], FillValue=default_tl);

            for i=1:length(radii)
                if (mod(i,10) == 0)
                    disp(num2str(i) + "/" + num2str(length(radii)));
                end
                if (layer_thickness == 0)
                    TL_collection(i) = obj.DesignDielectricGrinTL(pathlength(i));
                else
                    TL_collection(i) = obj.DesignDielectricGrinTLDiscretize(pathlength(i), layer_thickness);
                end
            end

            obj.lens_model.SetRadialDesign(TL_collection, radii);
        end

        function DesignDielectricGrinDiscretizeZ(obj, layer_thickness)
            obj.DesignDielectricGrinGeneral(layer_thickness);
        end

        function pathlength = FindIdealPathLength(obj, r)
            lamb_x = Constants.c/obj.sim.designfrequency;
            f = obj.sim.designfrequency;
            phasesim = Simulation(obj.lens_model, f, f);
            phasegrid = angle(phasesim.calc_ideal_phase());
            rgrid = sqrt(phasesim.x.^2 + phasesim.y.^2);
            rinterp = rgrid(round(phasesim.dims(1)/2), 1:round(phasesim.dims(2)/2));
            phaseinterp = phasegrid(round(phasesim.dims(1)/2), 1:round(phasesim.dims(2)/2));
            phaseinterp = unwrap(phaseinterp);

            phaser = interp1(rinterp, phaseinterp, r, 'spline');
            phase0 = interp1(rinterp, phaseinterp, 0, 'spline');

            pathlength = lamb_x * phaser / (2*pi);
            pathlength0 = lamb_x * phase0 / (2*pi);
            pathlength = pathlength - pathlength0;
        end

        function pathlength = FindNaivePathLength(obj, r)
            Rz = obj.sim.calc_phase_radius_of_curvature();
            pathlength = r.^2 / (2*Rz);
        end

        function r = idealCurvature(obj)
            dr = obj.lens_model.grid_dimension / 4;
            radii = -(obj.lens_model.diameter/2):dr:(obj.lens_model.diameter/2);

            pathlength = obj.FindIdealPathLength(abs(radii));

            p = polyfit(radii, pathlength, 2);
            r = 1/(2*p(1));
        end


        %% returns the index of the center part, thickness of the center part, and thickness of the AR coating part 
        function [n_of_r, t_of_r, ar_thick] = DesignIdealDielectric(obj, pathlength)
            t_tot = obj.total_thickness;
            
            n0 = 7; % maximum index of refraction for the lens

            lamb_x = Constants.c/obj.sim.designfrequency; % X-band (9 GHz) wavelength in mm
            
            f = obj.lens_model.focal_length; % focal length in millimeters

            % Rz = obj.sim.calc_phase_radius_of_curvature();

            % pathlength = r.^2/(2*Rz); % naive
            % pathlength = r.^2/(2*obj.idealCurvature()); % quadratic fit
            % pathlength = obj.FindIdealPathLength(dist);

            t0 = t_tot - 2.*lamb_x/(4. * sqrt(n0)); % Thickness of the high index region in the center of the lens
            
            %r = dist; %radius in millimeters
            
            b = (lamb_x.^2/(4.*t_tot^2)) + 2*(n0.*t0 - pathlength)./t_tot; %coefficient
            
            c = ((n0.*t0 - pathlength)./t_tot).^2; % coefficient
            
            n_of_r = (b + sqrt(b.^2 - 4.*c))/2.; % value of high index material as a function of radius
            
            t_of_r = (n0.*t0 - pathlength)./n_of_r; % thickness of high index material as a function of radius
            
            phase = 180.*(t_tot - t_of_r).*sqrt(n_of_r)/lamb_x; % phase of AR coating (should be 90 degrees)
            
            ar_thick = (t_tot - t_of_r)/2; % thickness of AR coating as a function of radius (index is the square root of the high index part)
            
            n_thick = t_of_r; % thickness of the high index part again
        end

        function out = DesignDielectricGrinTL(obj, pathlength)
            %DesignDielectricGrinTL Given a dist from the lens center,
            %design a single dielectric grin TL

            [n_of_r, n_thick, ar_thick] = obj.DesignIdealDielectric(pathlength);

            inner_dielectric = Dielectric(n_of_r.^2, 0, obj.sim.frequency, n_thick./1.2, n_thick.*1.2);

            % note we use n_inner as the AR dielectric because epsilson_ar is the
            % the square of n_ar, but n_ar is the square root of n_inner
            ar_dielectric = Dielectric(n_of_r, 0, obj.sim.frequency, ar_thick./1.2, ar_thick.*1.2);
    
            proto = ProtoTL_1Param([ar_dielectric, inner_dielectric, ar_dielectric]);
            TL = TL_1Param([ar_thick, n_thick, ar_thick], proto);

            out = TL;
        end

        function out = DesignDielectricGrinTLDiscretize(obj, pathlength, layer_thickness)
            %DesignDielectricGrinTL Given a dist from the lens center,
            %design a single dielectric grin TL

            [n_of_r, n_thick, ar_thick] = obj.DesignIdealDielectric(pathlength);

            ar_thick = layer_thickness * (round(ar_thick/layer_thickness));
            n_thick = obj.total_thickness - 2*ar_thick;

            assert(n_thick >= 0);
            assert(ar_thick >= 0);
            assert(n_thick + 2*ar_thick == obj.total_thickness);

            inner_dielectric = Dielectric(n_of_r.^2, 0, obj.sim.frequency, n_thick./1.2, n_thick.*1.2);

            % note we use n_inner as the AR dielectric because epsilson_ar is the
            % the square of n_ar, but n_ar is the square root of n_inner
            ar_dielectric = Dielectric(n_of_r, 0, obj.sim.frequency, ar_thick./1.2, ar_thick.*1.2);
    
            proto = ProtoTL_1Param([ar_dielectric, inner_dielectric, ar_dielectric]);
            TL = TL_1Param([ar_thick, n_thick, ar_thick], proto);

            out = TL;
        end
        
        % size is the size of meta-atom grid squares. epsr is the is
        % the relative permittivity of the 
        function out = GetIndexFromSpacing(obj, size, epsr, s, alpha, beta)
            g = obj.lens_model.grid_dimension;
            a = g - size; % the space between squares assuming constant size
            
            lambda = Constants.c/obj.sim.designfrequency;
            % capacitance of the meta-atom square
            C = (epsr * 4 * g / (2 * pi * Constants.c * Constants.Z0)) * ln(csc(a*pi / (2*g)));

            Zeff = sqrt( Constants.mu_0 / (Constants.eps_0 * epsr + (alpha * C / (s + (beta * lambda * ...
                a / g)))));

            out = Constants.Z0 / Zeff;
        end
    end
end

