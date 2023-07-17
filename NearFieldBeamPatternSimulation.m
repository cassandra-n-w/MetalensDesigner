% designed to be run AFTER Create480GHzLens.m

% replace lens material with aluminum
lensmodel.TL_prototype = proto_al;
freqs = 460:5:495;
freqs = freqs * 1e9;

plotfreq = f;

[minval, freqidx] = min(abs(freqs - plotfreq));

% from https://www.vadiodes.com/images/AppNotes/VDI_Feedhorn_Summary_2020.05.04.pdf
% 325-500 horn is L = 36.0mm and aperture diameter (2a) = 3.6mm
horn = DiagonalHornClass(3.6/2, 36);
hornfunc = @(x,y,f) EField(horn, x,y,f);

sim = Simulation(lensmodel, freqs, f);

% initialize incoming plane wave
sim.initialize_E_field(hornfunc);
horn_dist = focal_length-5.7;
sim.propagate(horn_dist);
resim = false;
if (resim)
    lensmodel.CalcSParam(freqs);
    % transform through ideal lens pattern

end

sim.lensTransform();

%% 
phaseweight = ((sim.xvec.^2 + sim.yvec.^2) < (lensmodel.diameter^2/4)) .* abs(sim.E_current(:,:,freqidx));
avgphase = sum(angle(sim.E_current(:,:,freqidx)).*phaseweight)/sum(phaseweight);

avg_phase_error = sum( abs(angle(sim.E_current(:,:,freqidx)) - avgphase) .* phaseweight) / sum(phaseweight);

subplot(1,2,1);
imagesc(sim.xvec, sim.yvec, mag2db(abs(sim.E_current(:,:,freqidx))));
colorbar;
subplot(1,2,2);
imagesc(sim.xvec, sim.yvec, mod(pi - avgphase + ((sim.xvec.^2 + sim.yvec.^2) < (lensmodel.diameter^2/4)) .* angle(sim.E_current(:,:,freqidx)), 2*pi));
colorbar;
sgtitle("Power (dB) and Phase (rad) of near-field scan for 480GHz lens at " + plotfreq/1e9 + "GHz, with horn at " + horn_dist + "mm distance")







