classdef Simulation
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
            
            xbounds = [-(dims(1)-1)/2, (dims(1)-1)/2];
            ybounds = [-(dims(2)-1)/2, (dims(2)-1)/2];

            [x, y] = meshgrid(linspace(xbounds(1), xbounds(2), dims(1)), linspace(ybounds(1), ybounds(2), dims(2)));

            obj.x = x * dx;
            obj.y = y * dy;

            
        end
        
        function phasepattern = calc_ideal_phase(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        function phasepattern = calc_gaussian_phase(obj)
            lambda = obj.c / (obj.designfrequency);
            r = sqrt(obj.x.^2 + obj.y.^2);
            mag = r < radius;
            
            w0 = ; % presumed beam waist of receiver in mm

            zr = pi * w0^2 / lambda; %confocal distance

            zf = obj.lens_model.focal_length; % focal distance of lens

            Rz = zf * (1 + (zr/zf)^2); %beam radius of curvature at focal distance

            phasetrans = exp(-1.0i * pi * r.^2 / (lambda * Rz));
            phasepattern = mag .* phasetrans;
        end
        
        function waist = calc_waist(obj)
            
        end
        
        function setup(obj)
            
        end
    end
    
    methods(Static)
        function [out] = Propagate(E, z, k)
            % note: designed for single frequency input, vectorized z input
            [A, kx, ky] = FourierEtoA(E, dx, dy);
            out = FourierAtoE(A, z, kx, ky, k);
        end
        
        function [A, kx, ky] = FourierEtoA(E, dx, dy)
        %FourierEtoA Convert an E-field [ E(x,y) ] to a vector potential A(kx, ky)
            %   Detailed explanation goes here
            % calculate vector potential
            % note the conjugate is taken here to account for phase convention
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

