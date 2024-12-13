
f = 9e9; % 480GHz central frequency
numfreqs=201;
freqs = linspace(f*0.5, f*6, numfreqs);

lambda = Constants.c/f;
%gridding size
g=lambda/20; % 0.11992mm ; this should be equal to lambda(500GHz, free space)/5

total_thickness = 16; %16mm total thickness

% diameter of active lens area and structural lens area, respectively
diameter = 600;
struct_diam = 600;

focal_length = 800;

% 0 is used as a placerholder for the prototype TL in this lens, since each
% TL has a different prototype
lensmodel = Lens(diameter, total_thickness,  focal_length, g, 0, struct_diam);
lensmodel.DeclareTL_1Param();

sim = Simulation(lensmodel, freqs, f);

designer = PhilDesigner(lensmodel, sim);

designer.DesignDielectricGrin();

%% 

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

[r, d, o] = sim.CalcLoss();

%%

plot(freqs, strehl_ratio, freqs, r, freqs, d, freqs, o);

legend(["Strehl ratio", "antireflective efficiency", "dielectric/conductive efficiency", "optical efficiency"]);
ylim([0,1]);

%%

passivity = sqrt(abs(lensmodel.S_array(:,:,100)).^2 + abs(lensmodel.S11_array(:,:,100).^2));

imagesc(lensmodel.xvec, lensmodel.yvec, passivity);
colorbar;

%%
phase = angle(lensmodel.S_array(:,:,40));
imagesc(lensmodel.xvec, lensmodel.yvec, phase);
colorbar;
