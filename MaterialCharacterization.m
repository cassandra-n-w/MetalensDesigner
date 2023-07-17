freqs = linspace(330e9, 500e9, 1000);
Z0 = 376.730314;

polyimide = Dielectric(3.40 * ones(size(freqs)), 0.015, freqs, 0.005, 1.0);
polyethylene = Dielectric(2.34 * ones(size(freqs)), 0.0004, freqs, 0.0001, 1.0);
polypropylene = Dielectric(2.25 * ones(size(freqs)), 0.0007, freqs, 0.005, 2.0);
ptfe = Dielectric(2.07 * ones(size(freqs)), 0.0006, freqs, 0.005, 1.0);
SiN = Dielectric(6.0 * ones(size(freqs)), 0.002, freqs, 150e-6, 450e-6);

% glued polyimide and polyethylene
abcd_total = pagemtimes(pagemtimes( ...
             polyethylene.abcd_interp(freqs, 0.425), ...
             polyimide.abcd_interp(freqs, 0.132)), ...
             polyethylene.abcd_interp(freqs, 0.425));

s_tot = abcd2s(abcd_total, Z0);

s11 = mag2db(abs(permute(s_tot(1,1,:), [3 2 1])));
subplot(1,3,1)
load("HDPEtestdata.mat");
plot(freqs/1e9, s11, f, mag_load);
xlabel("Frequency (GHz)");
ylabel("S11 (dB)");
ylim([min(s11)-5,0]);
xlim([min(freqs/1e9), max(freqs/1e9)]);
title("Polyimide core with glued HDPE");
legend(["Model", "Measurement"])

% glued polyimide and polypropylene
% initial nominal:
% polyimide 110um
% polypropylene 1.645mm

% much closer match:
% polyimide 147um
% polypropylene 1.616mm
abcd_total = pagemtimes(pagemtimes( ...
             polypropylene.abcd_interp(freqs, 1.646), ...
             polyimide.abcd_interp(freqs, 0.110)), ...
             polypropylene.abcd_interp(freqs, 1.646));

s_tot = abcd2s(abcd_total, Z0);

s11 = mag2db(abs(permute(s_tot(1,1,:), [3 2 1])));

abcd_total = pagemtimes(pagemtimes( ...
             polypropylene.abcd_interp(freqs, 1.616), ...
             polyimide.abcd_interp(freqs, 0.147)), ...
             polypropylene.abcd_interp(freqs, 1.616));

s_tot = abcd2s(abcd_total, Z0);

s11_2 = mag2db(abs(permute(s_tot(1,1,:), [3 2 1])));

subplot(1,3,2)
load("polypropylenetestdata.mat");
plot(freqs/1e9, s11, f, mag_load, freqs/1e9, s11_2);
xlabel("Frequency (GHz)");
ylabel("S11 (dB)");
ylim([min(s11)-5,0]);
xlim([min(freqs/1e9), max(freqs/1e9)]);
title("Polyimide core with glued Polypropylene");
legend(["Model (nominal)", "Measurement", "Model (best fit)"]);

% glued polyimide and polyethylene
abcd_total = pagemtimes(pagemtimes( ...
             polyethylene.abcd_interp(freqs, 0.26), ...
             polyimide.abcd_interp(freqs, 0.110)), ...
             polyethylene.abcd_interp(freqs, 0.26));

s_tot = abcd2s(abcd_total, Z0);

% s11 = mag2db(abs(permute(s_tot(1,1,:), [3 2 1])));
% subplot(2,2,3)
% plot(freqs/1e9, s11);
% xlabel("Frequency (GHz)");
% ylabel("S11 (dB)");
% ylim([min(s11)-5,0]);
% xlim([min(freqs/1e9), max(freqs/1e9)]);
% title("Polyimide core with glued Polyethylene (destroyed)");

% glued polyimide and polyethylene
p_thick = 0.180;
SiN_thick_layer = 300e-6;
p_thick_layer = (p_thick - SiN_thick_layer * 10) / 11;
abcd_total = polyimide.abcd_interp(freqs, p_thick_layer);

for i = 1:10
    abcd_total = pagemtimes(abcd_total, SiN.abcd_interp(freqs, SiN_thick_layer));
    abcd_total = pagemtimes(abcd_total, polyimide.abcd_interp(freqs, p_thick_layer));
end
    

s_tot = abcd2s(abcd_total, Z0);

s11 = mag2db(abs(permute(s_tot(1,1,:), [3 2 1])));

abcd_total = polyimide.abcd_interp(freqs, p_thick);

    

s_tot = abcd2s(abcd_total, Z0);

s11_2 = mag2db(abs(permute(s_tot(1,1,:), [3 2 1])));
%subplot(1,3,3)
subplot(1,1,1)
load("sintestdata.mat");
plot(freqs/1e9, s11_2, freqs/1e9, s11);
%plot(freqs/1e9, s11, f, mag_load, freqs/1e9, s11_2);
legend(["180um polyimide only", "180um with 10 x 300nm interlayers"])
xlabel("Frequency (GHz)");
ylabel("S11 (dB)");
ylim([min(s11_2)-5,0]);
xlim([min(freqs/1e9), max(freqs/1e9)]);
title("Polyimide core with/without Silicon Nitride interlayers");


