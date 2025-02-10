
f = 480e9; % 480GHz central frequency
numfreqs=50;
freqs = linspace(f*0.2, f*2, numfreqs);

lambda = Constants.c/f;
%gridding size
g=lambda/5; % 0.11992mm ; this should be equal to lambda(500GHz, free space)/5

% total_thickness = 0.3; %16mm total thickness
% 
% % diameter of active lens area and structural lens area, respectively
% diameter = 11.25;
% struct_diam = 11.25;
% 
% focal_length = 15;

total_thickness = 2.5; %16mm total thickness

% diameter of active lens area and structural lens area, respectively
diameter = 130;
struct_diam = 130;

focal_length = 150;

% 0 is used as a placerholder for the prototype TL in this lens, since each
% TL has a different prototype
lensmodel = Lens(diameter, total_thickness,  focal_length, g, 0, struct_diam);
lensmodel.DeclareTL_1Param();

sim = Simulation(lensmodel, freqs, f);

designer = PhilDesigner(lensmodel, sim);



%%
%designer.DesignDielectricGrinDiscretizeZ(total_thickness / 11);
designer.DesignDielectricGrin();

%% 

sim = Simulation(lensmodel, f, f);
lensmodel.CalcSParam(f);
% calculate the gaussian coupling of the lens

efield = @(x,y,f) (sqrt(x.^2 + y.^2) < (lensmodel.diameter/2)) * exp(0);
% initialize incoming plane wave
sim.initialize_E_field(efield);


% transform through ideal lens pattern
sim.lensTransform();

% propagate
scan_z = focal_length*linspace(0.8,1.1,100);

sim.propagate(scan_z);
E_save = sim.E_saved;

%strehl_ratio = sim.StrehlRatio();

sim2 = Simulation(lensmodel, freqs, f);
lensmodel.CalcSParam(freqs);
[r, d, o, strehl_ratio] = sim2.CalcLoss();


%%

passivity = sqrt(abs(lensmodel.S_array(:,:,round(numfreqs/2))).^2 + abs(lensmodel.S11_array(:,:,round(numfreqs/2)).^2));

imagesc(lensmodel.xvec, lensmodel.yvec, passivity);
colorbar;

%%
phase = angle(lensmodel.S_array(:,:,40));
imagesc(lensmodel.xvec, lensmodel.yvec, phase);
colorbar;

%%

plot(freqs, strehl_ratio, freqs, r, freqs, d, freqs, o);

legend(["Strehl ratio", "antireflective efficiency", "dielectric/conductive efficiency", "optical efficiency"]);
ylim([0,1]);

%%

E_size = size(E_save);
E_plot = squeeze(E_save(round(E_size(1)/2),:,:));

imagesc(scan_z, sim2.xvec, mag2db(abs(E_plot)));
colorbar;
clim([-35, max(mag2db(abs(E_plot)), [], "all")]);