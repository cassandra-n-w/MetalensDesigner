% frequency in hz
freq = 20e9;

% speed of light in m/s
c = 3e8;

% speed of light in mm/s
c = c * 1000;

% wavelength in mm
lambda = c / freq;

VerticalScanhand1 = readtable('C:\\Users\\Cassie\\Google Drive\\Dropbox\\Documents\\My Work\\Research\\metamaterial_lens\\Measurements\\VerticalScan_hand.xlsx');

height = VerticalScanhand1.ReceiverHeight_mm__1;

powerdB = VerticalScanhand1.Power__dBm__3;

powerdB = powerdB - max(powerdB);

powerlin = 10 .^ (powerdB./10);

% calculate "center of mass" of the len's measurement
com = sum(height .* powerlin) / sum(powerlin);
height = height - com;



% half-power beamwidth of horn antennas in degrees
hpbw_24 = 33;

% half-power beamwidth of horn antenna at 20GHz, calculated as hpbw_24 * sqrt( 10 ^ ((14.2 - 12.9)/10)) 
% where 14.2dBi is the 24GHz directivity and 12.9dBi is the 20GHz directivity
hpbw_20 = 38.3;

% (gaussian) divergence angle in radians
theta_24 = hpbw_24/1.18 * pi / 180;

% waist in mm, calculated at 24.125GHz because this is the nominal
% frequency of the horn. This should give 8.10mm at 24GHz 
lambda_24 = c / 24.15e9;
w_24 = lambda_24 /(theta_24 * pi);

% this should give a waist of 8.42mm at 20GHz 
theta_20 = hpbw_20/1.18 * pi/180;
w_20 = lambda / (theta_20 * pi);

xheights = -30:.1:30;
% gaussian shape of single horn beam
hornbeam = exp(-xheights.^2/w_20^2);

% gaussian shape of convolved horn beams
convbeam = exp(-xheights.^2 / (2 * w_20^2));

% fit by matlab
fitgauss = fit(height, powerlin, 'gauss1');

% fit by hand
param = 1.22;
wfit = param * w_20;
wfit = fitgauss.c1;
offset = fitgauss.b1;
fitbeam = exp(-(xheights).^2 / (2 * wfit^2));

[coupling, heights] = TestLensHornReal(lensmodel, f0);
%temp = find(coupling==max(coupling));
coupling = abs(coupling);
couplingoffset = sum(coupling .* heights) / sum(coupling);
couplingnorm = coupling/(max(coupling));


subplot(1,1,1)
p = plot(xheights, 20*log10(convbeam), height-offset, powerdB, xheights-offset, 10*log10(fitgauss(xheights)/fitgauss.a1), heights-couplingoffset, 20*log10(abs(couplingnorm)));
p(2).Marker = 'o';
p(2).LineStyle = 'none';
p(1).LineStyle = '--';
title("Image-plane Power Measurement");
xlabel("Offset from Lens Axis (mm)");
ylabel("Power (dB relative to peak)");
legend(["Ideal Focus", "Measurement", "Best Fit to Meas.", "EM Simulation"]);
ylim([-20,0]);
%plot(height, power);
