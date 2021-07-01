classdef HornClass
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Input parameters align with what is shown in Fig 13.16 of 
    %   Antenna Theory: Analaysis and Design by Balanis, 4th Ed.
    
    properties
        a1
        b1 
        p
        a
        b 
    end
    
    methods
        function obj = HornClass(a1, b1, p, a, b)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.a1 = a1;
            obj.b1 = b1;
            obj.p = p;
            obj.a = a;
            obj.b = b;
            
        end
        
        function out = EField(obj, x, y, f)
            [k, lambda] = Simulation.f_to_k_lambda(f);
                       
            rho1 = obj.p * obj.b1 / (obj.b1 - obj.b);
            rho2 = obj.p * obj.a1 / (obj.a1 - obj.a);
    
            mask = (abs(x) <= obj.a1/2) .* (abs(y) <= obj.b1/2);

            out = mask .* cos(pi*x/obj.a1) .* exp(-1.0i * k * (x.^2/rho2 + y.^2/rho1) / 2);
        end
    end
end

