function [coupling] = CalculateCoupling(EField1, EField2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    integrand = EField1 .* (EField2);
    
    coupling = trapz(trapz(integrand));
end

