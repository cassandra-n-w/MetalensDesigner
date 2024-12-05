% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


f = 480e9; % 480GHz central frequency
tracewidth = 0.005; %trace/space is 5 microns
numfreqs=31;
freqs = linspace(f*0.88, f*1.12, numfreqs);


%el = Element('C:\\Users\\Cassie\\OneDrive\\Drive\\Documents\\My Work\\Research\\metamaterial_lens\\code\\HFSS Models\\Polyimide 3_4 g0p2 v99 touchstone');
el = Element('..\\..\\HFSS Models\\Polyimide 3_4 g0p2 v99 touchstone');
di = Dielectric(3.50 * ones(size(freqs)), 0.015, freqs, 0.005, 0.015); %layer thickness 5 - 15um, dielectric 3.5, loss tangent 0.015


% elements with aluminum instead of copper
el_al = Element('..\\..\\HFSS Models\\Polyimide 3_4 Aluminum g0p2 v99 touchstone');

% silicon nitride dielectric
SiN = Dielectric(6.0 * ones(size(freqs)), 0.002, freqs, 150e-6, 450e-6);
di2_thickness = 300e-6;

proto = ProtoTL(el, di);
proto_al = ProtoTL(el_al, di);

%gridding size
g=0.11992; % 0.11992mm ; this should be equal to lambda(500GHz, free space)/5

layer_thickness = 0.010; %10 um layer thickness

% diameter of active lens area and structural lens area, respectively
diameter = 124;
struct_diam = 150;

focal_length = 150;

lensmodel = Lens(diameter, layer_thickness*ones(11, 1),  focal_length, g, proto, struct_diam);

sim = Simulation(lensmodel, freqs, f);

horn_waist = sim.calc_waist();

ideal = sim.calc_ideal_phase();
ideal = ideal(:,:,sim.designidx);

gaussapprox = sim.calc_gaussian_phase();

% plot the ideal phase of the lens
imagesc(sim.xvec, sim.yvec, 180+angle(ideal)*180/pi);
xlabel("X position (mm)");
ylabel("Y position (mm)");
colorbar;
title("Phase pattern over lens surface (degrees)");
xlim([-75, 75]);
ylim([-75, 75]);

% calculate the gaussian coupling of the lens

efield = @(x,y,f) (sqrt(x.^2 + y.^2) < (lensmodel.diameter/2)) * exp(0);
% initialize incoming plane wave
sim.initialize_E_field(efield);

% transform through ideal lens pattern
sim.transform(repmat(ideal,1,1,numfreqs));

% propagate
sim.propagate(focal_length);


waists = sim.FitGauss();

w0 = waists(sim.designidx);
gaussfunc = @(x, y, f) exp(-(x.^2+y.^2)/w0^2);
gaussfield = gaussfunc(sim.x, sim.y, f);
coupling = sim.calculate_coupling(gaussfield);

optim = Optimizer(lensmodel, tracewidth, f);
optim.Optim360();
optim.CreateGoalSparam();
optim.AssignOptimalTLs();

%% 

lensmodel.gerberfy()

%% 

for i = 1:length(optim.TLlist)
    TL = optim.TLlist(i);

    S21_1 = proto.SParamTL(f, TL);
    S21_2 = proto.SParamTL_Di(f, TL, SiN, di2_thickness);
    trans1(i) = abs(S21_1(2, 1));
    phase1(i) = angle(S21_1(2, 1));

    trans2(i) = abs(S21_2(2, 1));
    phase2(i) = angle(S21_2(2, 1));

end

subplot(2,1,1);
hold off
plot(trans1);
hold on
plot(trans2);
title("Transmittance with Silicon Nitride layers")
xlabel("TL index");
ylabel("Transmittance (unitless)");
legend(["As designed", "With SiN layers"]);


subplot(2,1,2);
hold off
plot(wrapToPi(phase1 - phase2) * 180/pi);
hold on
title("Phase error with Silicon Nitride layers")
xlabel("TL index");
ylabel("Phase Error (degrees)");

%% Plot transmission lines by phase

% first, single-frequency for all transmission lines optimized
S21_singfreq_trans1 = [];
S21_singfreq_phase1 = [];
for i = 1:length(optim.TLlist)-1
    TL = optim.TLlist(i);
    S21_1 = proto.SParamTL(f, TL);
    S21_singfreq_trans1(i) = abs(S21_1(2, 1));
    S21_singfreq_phase1(i) = angle(S21_1(2, 1));
end
hold off
subplot(1,1,1);
phase = ((1:length(optim.TLlist)-1) - 1) * 360/60;
yyaxis left;
plot(phase, mag2db(S21_singfreq_trans1));
ylabel("Optimized TL Transmission Magnitude (dB)")
yyaxis right;
plot(phase, (180/pi) * wrapTo2Pi(S21_singfreq_phase1));
hold on;
plot(phase,phase, '--');
xlim([0,360]);
xlabel("Nominal TL Phase (deg)");
ylim([0,360]);
ylabel("Optimized TL Output Phase (deg)");
legend(["TL Transmittance", "TL Phase", "Nominal Phase"]);
title("Optimized TL Performance at 480 GHz");

%% Plot transmission lines by radius

% first, single-frequency for all transmission lines optimized
S21_singfreq_trans1 = [];
S21_singfreq_phase1 = [];
radial_phase_ideal = [];

center_index = floor(size(lensmodel.TL_array) / 2);
center_index = center_index(1);


for i = 1:center_index
    TL = lensmodel.TL_array(center_index, i);
    S21_1 = proto.SParamTL(f, TL);
    S21_singfreq_trans1(i) = abs(S21_1(2, 1));
    S21_singfreq_phase1(i) = angle(S21_1(2, 1));
    radial_phase_ideal(i) = 180/pi * wrapTo2Pi(angle(optim.goalSparam(center_index, i)));
end

hold off
subplot(1,1,1);
phase = ((1:length(optim.TLlist)-1) - 1) * 360/60;
radius = (fliplr(1:center_index)-1)*lensmodel.diameter/(2*center_index);

yyaxis left;
plot(radius, mag2db(S21_singfreq_trans1));
ylabel("Optimized TL Transmission Magnitude (dB)")

yyaxis right;
plot(radius, (180/pi) * wrapTo2Pi(S21_singfreq_phase1));
hold on;

plot(radius,radial_phase_ideal, '--');

xlim([0,62]);
xlabel("Distance from Lens Center (mm)");
ylim([0,360]);
ylabel("Optimized TL Output Phase (deg)");
legend(["TL Transmittance", "TL Phase", "Nominal Phase"]);
title("Optimized TL Performance at 480 GHz");

%% Plot arbitrary transmission lines, freq-dependent

index = 16;
nominal_phase = (index - 1) * 360/60;
plotfreq = (451:1:499)*1e9;
TL_1 = optim.TLlist(index);
S_1 = proto.SParamTL(plotfreq, TL_1);

S21_1_trans = abs(S_1(2, 1,:));
S21_1_phase = 180/pi*angle(S_1(2, 1,:));

S11_1_trans = abs(S_1(1, 1,:));
S11_1_phase = 180/pi*angle(S_1(1, 1,:));

clf
yyaxis left
plot(plotfreq/1e9, mag2db(squeeze(S21_1_trans)));
hold on
plot(plotfreq/1e9, mag2db(squeeze(S11_1_trans)));
ylabel("Magnitude (dB)");
ylim([-35,0]);

yyaxis right
plot(plotfreq/1e9, wrapTo360(squeeze(S21_1_phase)));
plot(plotfreq/1e9, wrapTo360(squeeze(S11_1_phase)));
plot(plotfreq/1e9, ones(size(plotfreq)) * wrapTo360(nominal_phase));
ylabel("Phase (deg)");
ylim([0,360]);

xlabel("Frequency (GHz)");

legend(["S21 Magnitude", "S11 Magnitude", "S21 Phase", "S11 Phase", "Nominal Phase"]);
title(strcat("S-parameters for arbitrary TL (nominal phase ", num2str(nominal_phase), " deg)"));

%% Plot arbitrary meta-atoms

el_size = 0.95;
S_1 = el.Sparam_interp(plotfreq, el_size);


S21_1_trans = abs(S_1(2, 1,:));
S21_1_phase = 180/pi*angle(S_1(2, 1,:));

S11_1_trans = abs(S_1(1, 1,:));
S11_1_phase = 180/pi*angle(S_1(1, 1,:));

clf
yyaxis left
plot(plotfreq/1e9, mag2db(squeeze(S21_1_trans)));
hold on
plot(plotfreq/1e9, mag2db(squeeze(S11_1_trans)));
ylabel("Magnitude (dB)");
ylim([-35,0]);

yyaxis right
plot(plotfreq/1e9, wrapTo360(squeeze(S21_1_phase)));
plot(plotfreq/1e9, wrapTo360(squeeze(S11_1_phase)));
%plot(plotfreq/1e9, ones(size(plotfreq)) * wrapTo360(nominal_phase));
ylabel("Phase (deg)");
ylim([0,360]);

xlabel("Frequency (GHz)");
legend(["S21 Magnitude", "S11 Magnitude", "S21 Phase", "S11 Phase"]);
title(strcat("S-parameters for arbitrary meta-atom (s = ", num2str(size), " of pixel grid spacing)"));

%% plot sizes of optimized TLs by phase

s1 = [];
s2 = [];
s3 = [];
s4 = [];
s5 = [];

for i = 1:length(optim.TLlist)-1
    TL = optim.TLlist(i);
    
    s1(i) = TL.sizes(1);
    s2(i) = TL.sizes(2);
    s3(i) = TL.sizes(3);
    s4(i) = TL.sizes(4);
    s5(i) = TL.sizes(5);
end

clf
hold on;
plot(phase, s1*100);
plot(phase, s2*100);
plot(phase, s3*100);
plot(phase, s4*100);
plot(phase, s5*100);

xlabel("Nominal TL Phase (deg)");
xlim([0,360]);
ylabel("s:p ratio (percent)");
title("Variation of meta-atom sizes by layer for optimized TLs");

legend(["Layer 1 (outermost)", "Layer 2", "Layer 3", "Layer 4", "Layer 5 (innermost)"]);

%% plot sizes of meta-atoms by radius

s1 = [];
s2 = [];
s3 = [];
s4 = [];
s5 = [];

center_index = floor(size(lensmodel.TL_array) / 2);
center_index = center_index(1);


for i = 1:center_index
    TL = lensmodel.TL_array(center_index, i);
    
    s1(i) = TL.sizes(1);
    s2(i) = TL.sizes(2);
    s3(i) = TL.sizes(3);
    s4(i) = TL.sizes(4);
    s5(i) = TL.sizes(5);
end

radius = (fliplr(1:center_index)-1)*lensmodel.diameter/(2*center_index);

% for i = 1:length(optim.TLlist)-1
%     TL = optim.TLlist(i);
% 
%     s1(i) = TL.sizes(1);
%     s2(i) = TL.sizes(2);
%     s3(i) = TL.sizes(3);
%     s4(i) = TL.sizes(4);
%     s5(i) = TL.sizes(5);
% end

clf
hold on;
plot(radius, s1*100);
plot(radius, s2*100);
plot(radius, s3*100);
plot(radius, s4*100);
plot(radius, s5*100);

xlabel("Distance from Lens Center (mm)");
xlim([0,62]);
ylabel("s:p ratio (percent)");
title("Variation of meta-atom sizes by layer for optimized TLs");

legend(["Layer 1 (outermost)", "Layer 2", "Layer 3", "Layer 4", "Layer 5 (innermost)"]);