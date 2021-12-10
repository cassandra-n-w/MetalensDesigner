classdef Simulation < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c
        frequency
        designfrequency
        lens_model
        dx
        dy
        x
        y
        xvec
        yvec
        dims
        E_saved
        E_current
        designidx
        lensdims
        lensbounds_small
        lensbounds_big
    end
    
    methods
        function obj = Simulation(lens_model, frequencies, designfrequency)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.c = 299792458000;
            obj.lens_model = lens_model;
            obj.frequency = frequencies;
            obj.designfrequency = designfrequency;
            [d, ix] = min( abs( frequencies - designfrequency) );
            obj.designidx = ix;
            obj.dx = lens_model.grid_dimension;
            obj.dy = lens_model.grid_dimension;
            diam = lens_model.diameter;
            dim_mult = 1.5;
            
            % set up an x-y grid with the center at 0,0
            dims(1) = round(diam*dim_mult / (obj.dx));
            dims(2) = round(diam*dim_mult / (obj.dy));
            
            obj.lensdims(1) = round(diam/obj.dx);
            obj.lensdims(2) = round(diam/obj.dy);
            
            obj.dims = dims;
            
            xbounds = [-(dims(1)-1)/2, (dims(1)-1)/2];
            ybounds = [-(dims(2)-1)/2, (dims(2)-1)/2];
            
            
            obj.lensbounds_small = round((obj.dims - obj.lensdims) / 2) + 1;
            obj.lensbounds_big = obj.lensbounds_small + obj.lensdims - 1;

            [x, y] = meshgrid(linspace(xbounds(1), xbounds(2), dims(1)), linspace(ybounds(1), ybounds(2), dims(2)));

            obj.x = x * obj.dx;
            obj.y = y * obj.dy;
            
            obj.xvec = obj.x(1,:);
            obj.yvec = obj.y(:,1);

            obj.E_current = zeros(obj.dims(1), obj.dims(2), length(obj.frequency));
        end
        
        function phasepattern = calc_ideal_phase(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            

            quicksim = Simulation(obj.lens_model, obj.frequency, obj.designfrequency);

            w0 = quicksim.calc_waist();
            gaussfunc = @(x, y, f) exp(-(x.^2+y.^2)/w0^2);

            quicksim.initialize_E_field(gaussfunc);

            r = sqrt(obj.x.^2 + obj.y.^2);
            mag = r < (obj.lens_model.diameter/2);

            quicksim.propagate(obj.lens_model.focal_length);
            E = quicksim.E_current;
            theta = angle(E);

            phasepattern = mag .*exp(-1i*theta);

        end
        
        function phasepattern = calc_gaussian_phase(obj)
            lambda = obj.c / (obj.designfrequency);
            r = sqrt(obj.x.^2 + obj.y.^2);
            mag = r < (obj.lens_model.diameter/2);
            
            w0 = obj.calc_waist(); % presumed beam waist of receiver in mm

            zr = pi * w0^2 / lambda; %confocal distance

            zf = obj.lens_model.focal_length; % focal distance of lens

            Rz = zf * (1 + (zr/zf)^2); %beam radius of curvature at focal distance

            phasetrans = exp(+1.0i * pi * r.^2 / (lambda * Rz));
            phasepattern = mag .* phasetrans;
        end
        
        function waist = calc_waist(obj)
            % approximate the desired far-field beam angle of the receiver
            % horn as z_focus / (2 * diameter of lens)
            theta0 = (obj.lens_model.diameter) / (2*obj.lens_model.focal_length);
            [k, lambda] = Simulation.f_to_k_lambda(obj.designfrequency);
            waist =  lambda/(pi*theta0);
        end
        
        function initialize_E_field(obj, E_func)
            
            for i = 1:length(obj.frequency)
                f = obj.frequency(i);
                efield = E_func(obj.x, obj.y, f);

                obj.E_current(:,:,i) = Simulation.normalize(efield);
            end         
            
        end
        
        function propagate(obj, z)
            [ks, lambdas] = Simulation.f_to_k_lambda(obj.frequency);
            
            % saved_num is the number of E_field slices saved previously
            saved_num = size(obj.E_saved);
            if (length(saved_num) < 3)
                saved_num = 1;
            else
                saved_num = saved_num(3);
            end
            new_num = length(z);
            for i = 1:length(obj.frequency)
                k = ks(i);
                [save, final] = obj.PropagateField(obj.E_current(:,:,i), z, k);
                obj.E_current(:,:,i) = final;
                obj.E_saved(:,:,(saved_num+1) : (saved_num+new_num), i) = save;
            end      
            
        end
        
        % breaks down loss into reflection, dielectric, and optical
        function [reflect, dielectric, optical] = CalcLoss(obj)
            

            S11 = obj.lens_model.PadGeneral(obj.lens_model.S11_array, zeros(obj.dims), obj.frequency);
            S21 = obj.lens_model.PadMultiFreq(zeros(obj.dims), obj.frequency);

            efield = @(x,y,f) (sqrt(x.^2 + y.^2) < (obj.lens_model.diameter/2)) * exp(0);
            sim = Simulation(obj.lens_model, obj.frequency, obj.designfrequency);
            
            % initialize incoming plane wave
            sim.initialize_E_field(efield);

            % transform through reflection of lens
            sim.transform(S11);
            
            reflect = squeeze(1 - sim.calculate_power());
            sim = Simulation(obj.lens_model, obj.frequency, obj.designfrequency);
            sim.initialize_E_field(efield);
            sim.transform(S21);
            
            dielectric = squeeze(sim.calculate_power()) ./ reflect;
            
            optical = squeeze(sim.StrehlRatio());
            
            optical = optical ./ (dielectric .* reflect);
            
        end
        
        % single frequency!
        function scan = receiver_scan(obj, receiver_func, xscan, yscan)
            scan = zeros(length(xscan), length(yscan));
            for i = 1:length(xscan)
                for j = 1:length(yscan)
                    x_offset = xscan(i);
                    y_offset = yscan(j);
                    
                    scan(i,j) = obj.calculate_coupling(receiver_func(obj.x+x_offset, obj.y+y_offset, obj.designfrequency));
                end
            end
        end
        
        %single frequency!
        function coupling = calculate_coupling(obj, receiverfield)
            integrand = (obj.E_current .* (Simulation.normalize(receiverfield)));
    
            coupling = trapz(trapz(integrand));
        end
        
        function power = calculate_power(obj)
            power = Simulation.calc_power(obj.E_current);
        end
        
        function power = calculate_power_masked(obj, mask)
             power = Simulation.calc_power(obj.E_current .* mask);
        end
        
        % takes the current E field and transforms it through a field
        % defined by S21_field. This field is a function of x, y, frequency
        function transform(obj, S21_field)
            for i = 1:length(obj.frequency)
                obj.E_current(:,:,i) = obj.E_current(:,:,i) .* S21_field(:,:,i);
            end    
        end
        
        function lensTransform(obj)
            
            %obj.lens_model.CalcSParam(obj.frequency);
            S = obj.lens_model.PadMultiFreq(zeros(obj.dims), obj.frequency);
            
            obj.transform(S);
        end
        
        function [save, final] = PropagateField(obj, E, z, k)
            % note: designed for single frequency input, vectorized z input
            
            %final = zeros([size(E) size(z)]);
            
%             if (length(z) > 5)
%                 ztemp = z(6:end);
%                 z = z(1:5);
%                 
%                 recur_answer = obj.PropagateField(E, ztemp, k);
%                 save(:,:,6:end) = recur_answer;
%                 
%                 [A, kx, ky] = obj.FourierEtoA(E);
%                 save(:,:,5) = FourierAtoE(A, z, kx, ky, k);
%                 final = save;
%             else
            [A, kx, ky] = obj.FourierEtoA(E);
            save = FourierAtoE(A, z, kx, ky, k);
            final = save(:,:,end);

        end
        
        function [waist] = FitGauss(obj)
            waist = zeros( 'like', obj.frequency);
            for i = 1:length(obj.frequency)
                E = obj.E_current(:,:,i);
                middle = round(obj.dims(2)/2);
                Efit = (abs(E(:,middle)));
                gaussfit = fit(obj.xvec.', Efit, 'gauss1', 'StartPoint', [0.01, 0, 1.0]);
                waist(i) = gaussfit.c1;
            end
        end
        
        function [A, kx, ky] = FourierEtoA(obj, E)
        %FourierEtoA Convert an E-field [ E(x,y) ] to a vector potential A(kx, ky)
            %   Detailed explanation goes here
            % calculate vector potential
            % note the conjugate is taken here to account for phase convention
            dx = obj.dx;
            dy = obj.dy;
            
            A = fft2(conj(E));

            dims = flip(size(E) - [1, 1]);
            %dims = flip(size(E));

            % calculate kx and ky parts of the k vector
            %[kx, ky] = meshgrid(linspace(-dims(1)/2,dims(1)/2, dims(1)+1), linspace(-dims(2)/2,dims(2)/2, dims(2)+1));

            xn = 0:dims(1);
            yn = 0:dims(2);

            xn = xn - (xn > length(xn)/2) * length(xn);
            yn = yn - (yn > length(yn)/2) * length(yn);

            [kx,ky] = meshgrid(xn, yn);

            % this reorganizes the k vector components to align with how FFT is
            % calculated
            %kx = fftshift(fftshift(kx, 1), 2);
            %ky = fftshift(fftshift(ky, 1), 2);

            kx = 2*pi*kx / (dims(2)*dx);
            ky = 2*pi*ky / (dims(1)*dy);
        end
        
        
        function [E] = Unpad(obj, Ein)
            
            E = Ein(obj.lensbounds_small(1):obj.lensbounds_big(1), obj.lensbounds_small(2):obj.lensbounds_big(2),:);
            
        end
        
        function [strehl] = StrehlRatio(obj)
            
            ideallens = obj.calc_ideal_phase();
            actuallens = obj.lens_model.PadMultiFreq(zeros(obj.dims), obj.frequency);
            

            efield = @(x,y,f) (sqrt(x.^2 + y.^2) < (obj.lens_model.diameter/2)) * exp(0);
            sim = Simulation(obj.lens_model, obj.frequency, obj.designfrequency);
            
            % initialize incoming plane wave
            sim.initialize_E_field(efield);

            % transform through ideal lens pattern
            sim.transform(ideallens);

            % propagate
            sim.propagate(obj.lens_model.focal_length);
            
            idealfocus = sim.E_current;
            
            % initialize incoming plane wave
            sim.initialize_E_field(efield);

            % transform through ideal lens pattern
            sim.transform(actuallens);

            % propagate
            sim.propagate(obj.lens_model.focal_length);
            
            actualfocus = sim.E_current;

            idealpeak = max(abs(idealfocus), [], [1 2]);
            actualpeak = max(abs(actualfocus), [], [1 2]);
            
            strehl = (actualpeak./idealpeak).^2;
            
        end
    end
    
    methods(Static)
        function [k, lambda] = f_to_k_lambda(f)
            c = 299792458000;
            lambda = c./f;
            k = 2*pi./lambda;
        end
        
        function power = calc_power(efield)
            power = trapz(trapz(abs(efield .* conj(efield))));
        end
        
        function [normalized] = normalize(efield)
             power = abs(efield .* conj(efield));
             totalpower = trapz(trapz(power));
                
             normalized = efield ./ sqrt(totalpower);
        end
        
        function [E] = FourierAtoE(A, z, kx, ky, k0)
        %FourierAtoE Summary of this function goes here
        %   Detailed explanation goes here

            kz = sqrt(k0^2 - kx.^2 - ky.^2);
            zscale = permute(z, [3 1 2]);
            Az = (A .* exp(kz .* zscale * 1.0i));    

            % note the conjugate is taken to take into account phase convention
            E = conj(ifft2(Az));
        end
        
    end
end

