% designed to be run AFTER Create480GHzLens.m

% replace lens material with aluminum
lensmodel.TL_prototype = proto_al;

% from https://www.vadiodes.com/images/AppNotes/VDI_Feedhorn_Summary_2020.05.04.pdf
% 325-500 horn is L = 36.0mm and aperture diameter (2a) = 3.6mm
horn = DiagonalHornClass(3.6/2, 36);
hornfunc = @(x,y,f) EField(horn, x,y,f);

freqs = linspace(f*0.97, f*1.03, 7);



sim = Simulation(lensmodel, freqs, f);

% calculate the gaussian coupling of the lens

% initialize incoming plane wave
sim.initialize_E_field(hornfunc);

sim.propagate(focal_length*2);
resim = true;
if (resim)
    lensmodel.CalcSParam(freqs);
    % transform through ideal lens pattern

end

sim.lensTransform();

% propagate
sim.propagate(focal_length*2);

xscan = -10:0.2:10;
yscan = -10:0.2:10;

xcut = sim.receiver_scan(hornfunc, xscan, 0);
ycut = sim.receiver_scan(hornfunc, 0, yscan);

%% 
meas_dB = mag2db(abs(xcut(:,:,4)));
horn_loss = 1 - 0.0947;
horn_loss_dB = 10*log10(horn_loss);
meas_dB = meas_dB + 2*horn_loss_dB;
plot(xscan, meas_dB, yscan, mag2db(abs(ycut(:,:,4))) - 3);




