classdef PhilDesigner < handle
    %PHILDESIGNER Design a lens using Phil's equations
    %   Detailed explanation goes here
    
    properties
        lens_model
        sim
        total_thickness
        n0

        rinterp
        phaseinterp

        zone_bounds
        zone_jump

        layer_boundaries
        proto_tl

        tl_collection
        radii_collection

        cachedGetIndex

        is_metamaterial
    end
    
    methods
        function obj = PhilDesigner(lens_model, sim)
            %PHILDESIGNER Construct an instance of this class
            %   Detailed explanation goes here
            obj.lens_model = lens_model;
            obj.sim = sim;
            obj.total_thickness = obj.lens_model.layer_thicknesses;
            obj.n0 = 7;
            obj.is_metamaterial = false;

            obj.SimIdealPathLength();
        end
        
        function DesignDielectricGrin(obj)
            obj.DesignDielectricGrinGeneral(0);
        end

        function DesignDielectricGrinGeneral(obj, layer_count)
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
            pathlength = obj.PathLengthInterp(radii);

            default_tl = TL_1Param(0,0);
            TL_collection = createArray([length(radii),1], FillValue=default_tl);
            
            if isscalar(layer_count)
                layer_boundaries = (obj.total_thickness / layer_count) * (1:layer_count);
            elseif length(layer_count) == 2
                % treat first number as matching zone layer count, second
                % number as phase zone layer count
                matching_count = layer_count(1);
                phase_count = layer_count(2);
                [~, ~, matching_zone_thickness, ~] = obj.DesignIdealDielectric(pathlength(end));
                layer_boundaries = (matching_zone_thickness / matching_count) * (1:matching_count);
                phase_thickness = obj.total_thickness - 2*matching_zone_thickness;
                layer_boundaries = [layer_boundaries linspace(matching_zone_thickness + phase_thickness/phase_count, matching_zone_thickness + phase_thickness, phase_count)];
                layer_boundaries = [layer_boundaries layer_boundaries(end) + (matching_zone_thickness / matching_count) * (1:(matching_count - 1))];
            end

            obj.layer_boundaries = layer_boundaries;

            for i=1:length(radii)
                if (mod(i,10) == 0)
                    disp(num2str(i) + "/" + num2str(length(radii)));
                end
                if (length(layer_boundaries) == 0)
                    TL_collection(i) = obj.DesignDielectricGrinTL(pathlength(i));
                else
                    TL_collection(i) = obj.DesignDielectricGrinTLDiscretize(pathlength(i), layer_boundaries);
                end
            end
            obj.tl_collection = TL_collection;
            obj.radii_collection = radii;
            obj.lens_model.SetRadialDesign(TL_collection, radii);
        end

        function DesignDielectricGrinDiscretizeZ(obj, layer_count)
            obj.DesignDielectricGrinGeneral(layer_count);
        end

        function [index] = GetIndexFromPhaseSlope(obj, bgratio, quasi_layer_thickness)
            testfreqs = obj.sim.designfrequency * linspace(0.5, 1.0, 200);
            metal_thickness = 200e-6; %200nm converted to mm
            
            layer_count = 150;
            layer_thicknesses = (quasi_layer_thickness - metal_thickness) * ones([1, layer_count - 1]);
            layer_thicknesses = [(quasi_layer_thickness - metal_thickness)/2 layer_thicknesses (quasi_layer_thickness - metal_thickness)/2];

            TLlength = sum(layer_thicknesses);

            for i = 1:length(bgratio)
                squaresize = bgratio(i);
                metamaterialTL = TL(squaresize * ones([1, layer_count]), layer_thicknesses, obj.proto_tl);

                Sparams(:,:,:,i) = metamaterialTL.SParam(testfreqs);
                S21 = Sparams(2,1,:,i);
                phase = angle(S21);
                phase = unwrap(phase);
                
                delta_phase = phase(end) - phase(1);
                delta_freq = testfreqs(end) - testfreqs(1);

                index(i) = -(delta_phase / delta_freq) * Constants.c / (2 * pi * TLlength);
            end
            
            % iterate
            for i = 1:length(bgratio)
                Sparam_temp = s2s(Sparams(:,:,:,i), Constants.Z0, Constants.Z0/index(i));
                S21 = Sparam_temp(2,1,:);
                phase = angle(S21);
                phase = unwrap(phase);
                
                delta_phase = phase(end) - phase(1);
                delta_freq = testfreqs(end) - testfreqs(1);

                index(i) = -(delta_phase / delta_freq) * Constants.c / (2 * pi * TLlength);
            end 
        end


        function metamaterialTL = ConvertDielectricTLToMetamaterial(obj, dielectric_tl, tracespace)
            metal_thickness = 200e-6;
            % find thicknesses
            % convert boundaries to thicknesses and locations
            boundaries_extended = [0 obj.layer_boundaries obj.total_thickness];
            quasi_layer_thicknesses = diff(boundaries_extended);
            square_locations = (boundaries_extended(1:end-1) + boundaries_extended(2:end))/2;
            locations_extended = [0 square_locations obj.total_thickness];
            layer_thicknesses = diff(locations_extended) - metal_thickness;

            unique_quasi_thicknesses = unique(quasi_layer_thicknesses);

            % set up inverse square sizes
            g = obj.lens_model.grid_dimension;
            squaresize = linspace(tracespace, g-tracespace, 100);

            epsr = obj.proto_tl.Dielectric.Eps_r_grid(obj.sim.designfrequency);
            alpha = 1.15;
            beta = 0.07;
            
            for j = 1:length(unique_quasi_thicknesses)
                spacing = unique_quasi_thicknesses(j);
                %index(j, :) = obj.GetIndexFromSpacing(squaresize, epsr, spacing, alpha, beta);
                indexj = obj.cachedGetIndex(obj, squaresize/g, spacing);
                index(j, :) = indexj;
            end
            
            current_dielectric_index = 1;

            % create list of artificial dielectric boundaries
            for k = 1:length(dielectric_tl.my_proto_tl.TwoPortArray)
                artificial_dielectric = dielectric_tl.my_proto_tl.TwoPortArray(k);
                epsr_artificial_list(k) = artificial_dielectric.Eps_r_grid(obj.sim.designfrequency);
                artificial_thickness_list(k) = dielectric_tl.sizes(k);
            end

            artificial_boundaries_list = cumsum(artificial_thickness_list);
             
            assert(abs(artificial_boundaries_list(end) - obj.total_thickness) < 1e-6); % demand nm accuracy

            % translate list into an index for each square location
            % kinda janky.
            % basically, this should set an index of 1 for meta-atoms in
            % the first artificial dielectric zone, an index of 2 for meta-atoms
            % in the second artificial dielectric zone, etc.
            % this could probably be vectorized but idk how.
            metaatom_identity_index = ones(size(square_locations));
            for k = 1:length(artificial_boundaries_list)
                metaatom_identity_index = metaatom_identity_index + (square_locations > artificial_boundaries_list(k));
            end

            metaatom_size_list = zeros(size(square_locations));
            % find square sizes
            for i = 1:length(square_locations)
                % find desired artificial dielectric constant
                epsr_artificial = epsr_artificial_list(metaatom_identity_index(i));
                quasi_thickness_index = find(quasi_layer_thicknesses(i) == unique_quasi_thicknesses,1);
                index_artificial = sqrt(epsr_artificial);
                metaatom_size_list(i) = interp1(index(quasi_thickness_index,:), squaresize, index_artificial, 'linear', 'extrap');
                if metaatom_size_list(i) > (g - tracespace)
                    metaatom_size_list(i) = (g - tracespace);
                end
                if metaatom_size_list(i) < tracespace
                    metaatom_size_list(i) = tracespace;
                end
            end
            % need to convert to fraction of g for metaatom size
            metamaterialTL = TL(metaatom_size_list/obj.lens_model.grid_dimension, layer_thicknesses, obj.proto_tl);
            
        end
        
        function ConvertLensToMetamaterial(obj, proto_tl, tracespace)
            obj.is_metamaterial = true;
            default_tl = TL(0,0, proto_tl);
            TL_collection = createArray([length(obj.radii_collection),1], FillValue=default_tl);
            obj.proto_tl = proto_tl;
            
            function_handle = @GetIndexFromPhaseSlope;
            obj.cachedGetIndex = memoize(function_handle);

            %obj.tl_collection

            % convert each ideal dielectric tranmission line to a
            % metamaterial transmission line
            for i = 1:length(obj.radii_collection)
                if (mod(i,10) == 0)
                    disp(num2str(i) + "/" + num2str(length(obj.radii_collection)));
                end
                TL_collection(i) = obj.ConvertDielectricTLToMetamaterial(obj.tl_collection(i), tracespace);
            end

            obj.lens_model.SetRadialDesign(TL_collection, obj.radii_collection);
        end

        % function pathlength = FindIdealPathLength(obj, r)
        %     lamb_x = Constants.c/obj.sim.designfrequency;
        %     f = obj.sim.designfrequency;
        %     phasesim = Simulation(obj.lens_model, f, f);
        %     phasegrid = angle(phasesim.calc_ideal_phase());
        %     rgrid = sqrt(phasesim.x.^2 + phasesim.y.^2);
        %     rinterp = rgrid(round(phasesim.dims(1)/2), 1:round(phasesim.dims(2)/2));
        %     phaseinterp = phasegrid(round(phasesim.dims(1)/2), 1:round(phasesim.dims(2)/2));
        %     phaseinterp = unwrap(phaseinterp);
        %     obj.rinterp = rinterp;
        %     obj.phaseinterp = phaseinterp;
        % 
        %     phaser = interp1(rinterp, phaseinterp, r, 'spline');
        %     phase0 = interp1(rinterp, phaseinterp, 0, 'spline');
        % 
        %     pathlength = lamb_x * phaser / (2*pi);
        %     pathlength0 = lamb_x * phase0 / (2*pi);
        %     pathlength = pathlength - pathlength0;
        % end

        
        function SimIdealPathLength(obj)
            % use the fourier optics simulator to determine the ideal path
            % length from the phase output of an ideal gaussian beam. Store
            % the ideal path length for later use
            f = obj.sim.designfrequency;
            phasesim = Simulation(obj.lens_model, f, f);
            phasegrid = angle(phasesim.calc_ideal_phase());
            rgrid = sqrt(phasesim.x.^2 + phasesim.y.^2);
            rinterp = rgrid(round(phasesim.dims(1)/2), 1:round(phasesim.dims(2)/2));
            phaseinterp = phasegrid(round(phasesim.dims(1)/2), 1:round(phasesim.dims(2)/2));
            phaseinterp = unwrap(phaseinterp);
            obj.rinterp = rinterp;
            obj.phaseinterp = phaseinterp;
        end

        function pathlength = PathLengthInterp(obj, r)
            % Use the stored path length information to interpolate the
            % path length for a given radius range. Include fresnel zone
            % jumps.
            pathlength = obj.PathLengthInterpRaw(r);

            %lamb_x = Constants.c/obj.sim.designfrequency;

            if length(obj.zone_bounds) > 0
                pathmat = repmat(pathlength', 1, length(obj.zone_bounds));
                boundmat = repmat(obj.zone_bounds, length(r), 1);
    
                zone_boolean = pathmat > boundmat;
                correction = squeeze(zone_boolean * obj.zone_jump');
            else
                correction = 0;
            end
            
            pathlength = pathlength - correction';
        end

        function pathlength = PathLengthInterpRaw(obj, r)
            % Use the stored path length information to interpolate the
            % path length for a given radius range.
            lamb_x = Constants.c/obj.sim.designfrequency;
            phaser = interp1(obj.rinterp, obj.phaseinterp, r, 'spline');
            phase0 = interp1(obj.rinterp, obj.phaseinterp, 0, 'spline');

            pathlength = lamb_x * phaser / (2*pi);
            pathlength0 = lamb_x * phase0 / (2*pi);
            pathlength = pathlength - pathlength0;
        end
        
        function SetZoneBounds(obj, bounds, jumpsizes)
            % bounds are the locations (in mm pathlength) of the jumps, jumpsizes
            % is how many mm to jump.
            obj.zone_bounds = bounds;
            obj.zone_jump = jumpsizes;
        end

        function ZoneBoundsEveryTwoPi(obj)
            % example of zone testing, places z bound every two pi
            pathmax = obj.PathLengthInterpRaw(obj.lens_model.diameter/2);
            lamb_x = Constants.c/obj.sim.designfrequency;

            jumpcount = floor(pathmax / (lamb_x));
            
            for i = 1:jumpcount
                bounds(i) = i * lamb_x;
                jump(i) = lamb_x;
            end

            obj.SetZoneBounds(bounds, jump);
        end

        function pathlength = FindNaivePathLength(obj, r)
            Rz = obj.sim.calc_phase_radius_of_curvature();
            pathlength = r.^2 / (2*Rz);
        end

        function r = idealCurvature(obj)
            dr = obj.lens_model.grid_dimension / 4;
            radii = -(obj.lens_model.diameter/2):dr:(obj.lens_model.diameter/2);

            pathlength = obj.PathLengthInterpRaw(abs(radii));

            p = polyfit(radii, pathlength, 2);
            r = 1/(2*p(1));
        end

        
        % returns the index of the center part, thickness of the center part, and thickness of the AR coating part 
        function [n_of_r, t_of_r, ar_thick, t0] = DesignIdealDielectric(obj, pathlength)
            t_tot = obj.total_thickness;
            
            %n0 = 7; % maximum index of refraction for the lens

            lamb_x = Constants.c/obj.sim.designfrequency; % X-band (9 GHz) wavelength in mm
            
            f = obj.lens_model.focal_length; % focal length in millimeters

            % Rz = obj.sim.calc_phase_radius_of_curvature();

            % pathlength = r.^2/(2*Rz); % naive
            % pathlength = r.^2/(2*obj.idealCurvature()); % quadratic fit
            % pathlength = obj.FindIdealPathLength(dist);

            t0 = t_tot - 2.*lamb_x/(4 .* sqrt(obj.n0)); % Thickness of the high index region in the center of the lens
            
            %r = dist; %radius in millimeters
            
            b = (lamb_x.^2/(4.*t_tot^2)) + 2*(obj.n0.*t0 - pathlength)./t_tot; %coefficient
            
            c = ((obj.n0.*t0 - pathlength)./t_tot).^2; % coefficient
            
            n_of_r = (b + sqrt(b.^2 - 4.*c))/2.; % value of high index material as a function of radius
            
            t_of_r = (obj.n0.*t0 - pathlength)./n_of_r; % thickness of high index material as a function of radius
            
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

        function out = DesignDielectricGrinTLDiscretize(obj, pathlength, layer_boundaries)
            %DesignDielectricGrinTL Given a dist from the lens center,
            %design a single dielectric grin TL

            [n_of_r, n_thick, ar_thick, t0] = obj.DesignIdealDielectric(pathlength);

            pathlength_central = n_thick .* n_of_r;

            % round to nearest layer thicknesses
            [error, roundindex] = min(abs(ar_thick - layer_boundaries));
            ar_thick = layer_boundaries(roundindex);
            n_thick = obj.total_thickness - 2*ar_thick;            

            assert(n_thick >= 0);
            assert(ar_thick >= 0);
            assert(n_thick + 2*ar_thick == obj.total_thickness);

            % correct for discretization of phase section
            n_of_r = pathlength_central ./ n_thick;

            % correct for discretization of matching section
            lamb_x = Constants.c/obj.sim.designfrequency;
            n_ar = lamb_x/(4*ar_thick);

            inner_dielectric = Dielectric(n_of_r.^2, 0, obj.sim.frequency, n_thick./1.2, n_thick.*1.2);

            % note we use n_inner as the AR dielectric because epsilson_ar is the
            % the square of n_ar, but n_ar is the square root of n_inner
            ar_dielectric = Dielectric(n_ar.^2, 0, obj.sim.frequency, ar_thick./1.2, ar_thick.*1.2);
    
            proto = ProtoTL_1Param([ar_dielectric, inner_dielectric, ar_dielectric]);
            TL = TL_1Param([ar_thick, n_thick, ar_thick], proto);

            out = TL;
        end
        
        % size is the size of meta-atom grid squares. epsr is the is
        % the relative permittivity of the actual bulk dielectric
        function out = GetIndexFromSpacing(obj, size, epsr, s, alpha, beta)
            g = obj.lens_model.grid_dimension;
            a = g - size; % the space between squares assuming constant size
            
            lambda = Constants.c/obj.sim.designfrequency;
            % capacitance of the meta-atom square
            C = (epsr * 4 * g / (2 * pi * Constants.c * Constants.Z0)) * log(csc(a*pi / (2*g)));

            Zeff = sqrt( Constants.mu_0 ./ (Constants.eps_0 .* epsr + (alpha .* C ./ (s + (beta * lambda * ...
                a ./ g)))));

            out = Constants.Z0 ./ Zeff;
        end
    end
end

