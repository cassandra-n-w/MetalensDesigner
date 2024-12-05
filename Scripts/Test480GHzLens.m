% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


% designed to be run AFTER Create480GHzLens.m

% Temporary!!!
% this is to test the effect of replacing copper with aluminum, and the
% dielectric with a lossier one
el_al_thin = Element('..\\..\\HFSS Models\\Polyimide 3_4 200nm Aluminum g0p2 v99 touchstone');
di_lossier = Dielectric(3.5 * ones(size(freqs)), 0.025, freqs, 0.005, 0.015); %layer thickness 5 - 15um, dielectric 3.5, loss tangent 0.015
lensmodel.TL_prototype = ProtoTL(el_al_thin, di_lossier);

%test effect of thickness being greater than expected
thickness_mult = 1.14; %11.1 um layer thickness
%thickness_mult = 1.0; %10 um layer thickness
lensmodel.SetThicknessMult(thickness_mult);

freqs = linspace(460e9, 499.9e9, 9);

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


% simulate far-field from illumination from focal point
resim = false;
if (resim)
    sim = Simulation(lensmodel, freqs, f);
    % waveguide dimensions for VDI wr-2.2
    horn = RectangularWaveguideTE10(0.570,0.285);
    hornfunc = @(x,y,f) EField(horn, x,y,f);

    sim.initialize_E_field(hornfunc);

    sim.propagate(lensmodel.focal_length);

    sim.lensTransform();

    sim.NormalizeField();

    [farfield, sintx, sinty] = sim.FarField();
    pattern_db = mag2db(abs(farfield(:,:,5)));
    % imagesc(sintx(:,5),sinty(:,5),pattern_db);
    % colorbar;
    % clim([0, max(pattern_db, [],"all")]);
    % xlim([-0.06,0.06]);
    % ylim([-0.06,0.06]);
    
    hold on;
    xcut = pattern_db(:,round(size(pattern_db, 2)/2));
    anglesx = (180/pi)*real(asin(sintx(:,5)));
    plot(anglesx,xcut-max(xcut), "--");
    ycut = pattern_db(round(size(pattern_db, 1)/2),:);
    anglesy = (180/pi)*real(asin(sinty(:,5)));
    plot(anglesy,ycut-max(ycut), "--");
    dcut = diag(pattern_db);
    angles = anglesx*sqrt(2);
    plot(angles,dcut-max(dcut), "--");
    xlim([-3,3]); 
    legend(["H-plane (meas)", ...
            "D-plane (meas)", ...
            "E-plane (meas)", ...
            "H-plane (sim)", ...
            "E-plane (sim)", ...
            "D-plane (sim)"]);

    max(xcut);

    [max_wave, ix] = max(farfield);
    [max_wave, iy] = max(max_wave);
    % ix = ix(:,iy,:);
    % 
    % radsizex = diff(asin(sintx));
    % radsizex = radsizex(ix, :);
    % radsizey = diff(asin(sinty));
    % radsizey = radsizey(iy, :);
 
    radsizex = diff(asin(sintx));
    radsizex = squeeze(radsizex(round(size(farfield,1)/2),:)); 
    radsizey = diff(asin(sinty));
    radsizey = squeeze(radsizey(round(size(farfield,2)/2),:));
    
    max_wave = squeeze(max_wave)';
    directivity = abs(((4*pi)) .* max_wave);
    directivity_db = mag2db(directivity);
    total_power = trapz(asin(sinty(:,5)),trapz(asin(sintx(:,5)), farfield.*conj(farfield) ));
end
%% 
hold on;
% data taken directly from Alain's radiometric test.
freq_meas = (460:5:500);
loss_meas_dB = [2.76, 2.77, 2.73, 2.55, 2.52, 2.42, 2.38, 2.33, 2.34];

[r,d,o] = sim.CalcLoss()

waists = sim.FitGauss();

loss = r.*d
loss_dB = -10*log10(loss)

%plot(freq_meas, loss_meas_dB, freqs/1e9, loss_dB);
plot(freqs/1e9, loss_dB);
legend(["Measured", "Simulated (nominal)", "Simulated (realistic)"]);
xlabel("Frequency (GHz)");
ylabel("Loss (dB)")
title("Lens Loss: Simulated vs Measured")
ylim([0, 3]);
hold off;

%% 

hold on;
[ks, lambdas] = Simulation.f_to_k_lambda(freqs);
diff_limit = diameter*diameter*(pi/4) * (4*pi) ./ (lambdas.*lambdas);
diff_limit_dB = 10*log10((diff_limit)');
optical_dB = -10*log10(o);
direc_sim_dB = diff_limit_dB - optical_dB;
plot(freqs/1e9, diff_limit_dB - optical_dB,"Color","#D95319"); 
title("Directivity: Simulated vs Measured");
ylabel("Directivity (dB)");
legend(["Diff. Limit", "Directivity (meas)", "", "Directivity (sim)"]);

hold off;

%% 

w0 = waists(sim.designidx);
gaussfunc = @(x, y, f) exp(-(x.^2+y.^2)/w0^2);
gaussfield = gaussfunc(sim.x, sim.y, f);
coupling = sim.calculate_coupling(gaussfield);
do_plot = false;
if (do_plot)
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
