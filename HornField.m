function [EField] = HornField(a1, b1, p, a, b, k, x, y)
%UNTITLED Summary of this function goes here
%   Input parameters align with what is shown in Fig 13.16 of 
%   Antenna Theory: Analaysis and Design by Balanis, 4th Ed.
    
    % for a pyramidal horn, we assume p_e = p_h = p
    % we then use properties of similar triangles to find the
    % center-lengths
    rho1 = p * b1 / (b1 - b);
    rho2 = p * a1 / (a1 - a);
    
    mask = (abs(x) <= a1/2) .* (abs(y) <= b1/2);

    EField = mask .* cos(pi*x/a1) .* exp(-1.0i * k * (x.^2/rho2 + y.^2/rho1) / 2);
    
end

