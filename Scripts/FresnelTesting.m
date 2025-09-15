lambdas= linspace(1-0.12, 1.12, 1000);

for i = 1:length(lambdas)

    numpoints = 100;
    OPDr = linspace(0, 19.9, numpoints);
    x = 1:numpoints;
    
    lambda0 = 1;
    lambda = lambdas(i);
    
    OPDr_mod = mod(OPDr, 10);
    
    %plot(x, OPDr, x, OPDr_mod);
    
    phase_desired = mod(2*pi * OPDr / lambda, 2*pi);
    phase_actual = mod(2 * pi * OPDr_mod / lambda, 2*pi);
    phase_error = abs(phase_desired - phase_actual);
    phase_error = min(phase_error, 2*pi - phase_error);

    phase_error_l(i) = max(phase_error);

    %plot(x, phase_desired, x, phase_actual, x, phase_error);
end

plot(lambdas, phase_error_l, lambdas, ones(size(lambdas)) * pi/4);

%%
numpoints = 1000;
    OPDr = linspace(0, 20, numpoints);
    x = 1:numpoints;
    
    lambda0 = 1;
    lambda = 0.98;

    OPDr_mod = mod(OPDr, 10);
    
    plot(x, OPDr, x, OPDr_mod);
    
    phase_desired = mod(2*pi * OPDr / lambda, 2*pi);
    phase_actual = mod(2 * pi * OPDr_mod / lambda, 2*pi);
    phase_error = mod(phase_desired - phase_actual + pi, 2*pi) - pi;
    %phase_error = min(phase_error, 2*pi - phase_error);

    phase_error_l(i) = max(phase_error);

    plot(x, phase_desired, x, phase_actual, x, phase_error);
    legend(["Desired phase", "Actual Phase", "Phase Error"]);
    xlabel("Arbitrary units");
    ylabel("Phase (rad)");