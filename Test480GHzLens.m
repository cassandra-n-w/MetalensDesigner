% designed to be run AFTER Create480GHzLens.m

freqs = linspace(f*0.97, f*1.03, 7);

resim = false;
if (resim)
    sim = Simulation(lensmodel, freqs, f);

    % calculate the gaussian coupling of the lens

    efield = @(x,y,f) (sqrt(x.^2 + y.^2) < (lensmodel.diameter/2)) * exp(0);
    % initialize incoming plane wave
    sim.initialize_E_field(efield);

    lensmodel.CalcSParam(freqs);
    % transform through ideal lens pattern
    sim.lensTransform();

    % propagate
    sim.propagate(focal_length);
end


strehl_ratio = sim.StrehlRatio();

waists = sim.FitGauss();

w0 = waists(sim.designidx);
gaussfunc = @(x, y, f) exp(-(x.^2+y.^2)/w0^2);
gaussfield = gaussfunc(sim.x, sim.y, f);
coupling = sim.calculate_coupling(gaussfield);

% plot the ideal phase of the lens
S = lensmodel.S_array(:,:,1);
imagesc(sim.xvec, sim.yvec, 180+angle(S)*180/pi);
%imagesc(sim.xvec, sim.yvec, abs(S));
xlabel("X position (mm)");
ylabel("Y position (mm)");
colorbar;
title("Phase pattern over lens surface (degrees)");
xlim([-75, 75]);
ylim([-75, 75]);
