% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


% designed to be run AFTER Create480GHzLens.m

% replace lens material with aluminum
lensmodel.TL_prototype = proto_al;
freqs = 450:5:500;
freqs = freqs * 1e9;

plotfreq = f;

[minval, freqidx] = min(abs(freqs - plotfreq));

% from https://www.vadiodes.com/images/AppNotes/VDI_Feedhorn_Summary_2020.05.04.pdf
% 400-600 horn is L = 31.0mm and aperture diameter (2a) = 3.1mm
horn = DiagonalHornClass(3.1/2, 31);
hornfunc = @(x,y,f) EField(horn, x,y,f);

sim = Simulation(lensmodel, freqs, f);

% initialize incoming plane wave
sim.initialize_E_field(hornfunc);

power_init = sim.calculate_power();

sim.propagate(144);

sim.lensCircularMaskTransform();

power_overspill = sim.calculate_power();

resim = true;
if (resim)
    lensmodel.CalcSParam(freqs);
    % transform through ideal lens pattern

end

sim.lensTransform();

power_afterlens = sim.calculate_power();

%% 

power_afterlens = reshape(power_afterlens, 1, []);

% room temp in K at 68F/20C
T_room = 273.15 + 21;
clf
plot(freqs, 1-power_afterlens)
plot(freqs, (1 - power_afterlens) * T_room)
xlabel("Frequency (Hz)")
ylabel("Effective Noise Temperature (K)");
sgtitle("Effective noise temperature of Lens (based on loss from overillumination, reflection, and conductive/dielectric loss)");




