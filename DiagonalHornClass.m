classdef DiagonalHornClass
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Based on Johansson 1992, The diagonal horn as a sub-mmwave antenna
    
    properties
        a
        L
    end
    
    methods
        function obj = HornClass(a, L)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.a = a;
            obj.L = L;
            
        end
        
        % we use eta and xi here because they align with the polarization
        % directions. eta is polarized axis and xi is cross-pol
        function out = EField(obj, eta, xi, f)
            [k, lambda] = Simulation.f_to_k_lambda(f);

            % transform to x/y coordinates
            x = (eta + xi)/sqrt(2);
            y = (eta - xi)/sqrt(2);

            mask = (abs(x) <= obj.a) .* (abs(y) <= obj.a);
                       
            k_del = (2*pi ./ lambda) * ((2* (obj.a).^2 - x.^2 - y.^2))/ (2*obj.L);

            E_eta_mag = 1/(sqrt(2)) * (cos(pi*y / (2*obj.a)) + cos(pi*x/(2*obj.a)));

            out = mask .* E_eta_mag .* exp(1.0i * k_del);
        end
    end
end

