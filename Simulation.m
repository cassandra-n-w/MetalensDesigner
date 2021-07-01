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
        dims
        E_saved
        E_current
    end
    
    methods
        function obj = Simulation(lens_model, frequencies, designfrequency)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.c = 299792458000;
            obj.lens_model = lens_model;
            obj.frequency = frequencies;
            obj.designfrequency = designfrequency;
            obj.dx = lens_model.grid_dimension;
            obj.dy = lens_model.grid_dimension;
            diam = lens_model.diameter;
            
            % set up an x-y grid with the center at 0,0
            dims = diam*4 / (obj.dx);
            obj.dims = dims;
            
            xbounds = [-(dims(1)-1)/2, (dims(1)-1)/2];
            ybounds = [-(dims(2)-1)/2, (dims(2)-1)/2];

            [x, y] = meshgrid(linspace(xbounds(1), xbounds(2), dims(1)), linspace(ybounds(1), ybounds(2), dims(2)));

            obj.x = x * dx;
            obj.y = y * dy;

            obj.E_current = zeros(obj.dims(1), obj.dims(2), obj.frequency);
        end
        
        function phasepattern = calc_ideal_phase(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            quicksim = Simulation(obj.lens_model, obj.designfrequency, obj.designfrequency);
            
            quicksim.setup();
            
            w0 = quicksim.calc_waist();
            gaussfunc = @(x, y, f) exp(-(x.^2+y.^2)/w0^2);
            
            quicksim.initialize_E_field(gaussfunc);
            
            r = sqrt(obj.x.^2 + obj.y.^2);
            mag = r < (obj.lens_model.diameter/2);
            
            quicksim.propagate(obj.lens_model.focal_length);
            E = quicksim.E_current;
            
            phasepattern = mag .* ang(E);
        end
        
        function phasepattern = calc_gaussian_phase(obj)
            lambda = obj.c / (obj.designfrequency);
            r = sqrt(obj.x.^2 + obj.y.^2);
            mag = r < (obj.lens_model.diameter/2);
            
            w0 = obj.calc_waist(); % presumed beam waist of receiver in mm

            zr = pi * w0^2 / lambda; %confocal distance

            zf = obj.lens_model.focal_length; % focal distance of lens

            Rz = zf * (1 + (zr/zf)^2); %beam radius of curvature at focal distance

            phasetrans = exp(-1.0i * pi * r.^2 / (lambda * Rz));
            phasepattern = mag .* phasetrans;
        end
        
        function waist = calc_waist(obj)
            % approximate the desired far-field beam angle of the receiver
            % horn as z_focus / (2 * diameter of lens)
            theta0 = obj.lens_model.focal_length / (2*obj.lens_model.diameter);
            
            waist = pi*theta0 / obj.designfrequency;
        end
        
        function initialize_E_field(obj, E_func)
            
            for i = 1:length(obj.frequencies)
                f = obj.frequencies(i);
                obj.E_current(:,:,i) = E_func(obj.x, obj.y, f);
            end         
            
        end
        
        function propagate(obj, z)
            [ks, lambdas] = Simulation.f_to_k_lambda(obj.frequencies);
            
            % saved_num is the number of E_field slices saved previously
            saved_num = size(obj.E_saved);
            saved_num = saved_num(3);
            new_num = length(z);
            for i = 1:length(obj.frequencies)
                k = ks(i);
                [save, final] = obj.PropagateField(obj.E_current, z, k);
                obj.E_current(:,:,i) = final;
                obj.E_saved(:,:,(saved_num+1) : (saved_num+new_num), i) = save;
            end      
            
        end
        
        % takes the current E field and transforms it through a field
        % defined by S21_field. This field is a function of x, y, frequency
        function transform(obj, S21_field)
            for i = 1:length(obj.frequencies)
                obj.E_current(:,:,i) = obj.E_current(:,:,i) .* S21_field(:,:,i);
            end    
        end
        
        function [save, final] = PropagateField(obj, E, z, k)
            % note: designed for single frequency input, vectorized z input
            [A, kx, ky] = obj.FourierEtoA(E);
            save = FourierAtoE(A, z, kx, ky, k);
            final = save(:,:,end);
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
    end
    
    methods(Static)
        function [k, lambda] = f_to_k_lambda(f)
            c = 299792458000;
            lambda = c/f;
            k = 2*pi/lambda;
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

