% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


% designed to be run AFTER Create480GHzLens.m

% Temporary!!!
% this is to test the effect of replacing copper with aluminum
lensmodel.TL_prototype = proto_al;

freqs = linspace(f*0.97, f*1.03, 7);

resim = true;
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
    
    strehl_ratio = sim.StrehlRatio();
end

resim = false;

if (resim)
    zlist = linspace(0,focal_length,300);
    
    sim = Simulation(lensmodel, f, f);

    efield = @(x,y,f) (sqrt(x.^2 + y.^2) < (lensmodel.diameter/2)) * exp(0);
    % initialize incoming plane wave
    sim.initialize_E_field(efield);

    % transform through ideal lens pattern
    sim.lensTransform();

    % propagate
    sim.propagate(zlist);
end

[r,d,o] = sim.CalcLoss()

waists = sim.FitGauss();

w0 = waists(sim.designidx);
gaussfunc = @(x, y, f) exp(-(x.^2+y.^2)/w0^2);
gaussfield = gaussfunc(sim.x, sim.y, f);
coupling = sim.calculate_coupling(gaussfield);
plot = false;
if (plot)
    % plot the ideal phase of the lens
    hold off;
    S = lensmodel.S_array(:,:,sim.designidx);
    S2 = lensmodel.S_array(:,:,1);
    middle = size(sim.E_saved);
    middle = round(middle(2)/2);
    crosssection = squeeze(sim.E_saved(:,middle,2:end,sim.designidx));
    S = sim.E_current(:,:,sim.designidx);
    %S = sim.E_saved(:,:,end,sim.designidx);
    %S = S/max(S, [], [1 2]);
    %imagesc(sim.xvec, sim.yvec, 180+(angle(S) - angle(S2))*180/pi);
    %imagesc(sim.xvec, sim.yvec, abs(S).^2);
    intensity = abs(crosssection).^2;
    normalized = intensity ./ smoothdata(max(abs(crosssection).^2, [], 1));
    imagesc(zlist(2:end), sim.xvec, pow2db(intensity));
    xlabel("Z position (mm)");
    ylabel("Y position (mm)");
    colorbar;
    title("Cross-sectional Wave Propagation (normalized intensity in dB)");
    %xlim([-3 3]);
    %ylim([-3 3]);
    caxis([-70 -15]);

    plot(freqs/10^9, squeeze(strehl_ratio));
    xlabel("Frequency (GHz)");
    ylabel("Strehl Ratio (dimensionless)");
    title("Strehl Ratio as a function of frequency at z = 150mm");

end
