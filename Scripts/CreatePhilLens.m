
f = 9e9; % 480GHz central frequency
numfreqs=101;
freqs = linspace(f*0.5, f*2, numfreqs);

lambda = Constants.c/f;
%gridding size
g=lambda/10; % 0.11992mm ; this should be equal to lambda(500GHz, free space)/5

total_thickness = 16; %10 um layer thickness

% diameter of active lens area and structural lens area, respectively
diameter = 600;
struct_diam = 600;

focal_length = 800;

TLblah = TL_absorber();

% 0 is used as a placerholder for the prototype TL in this lens, since each
% TL has a different prototype
lensmodel = Lens(diameter, total_thickness,  focal_length, g, 0, struct_diam);
lensmodel.DeclareTL_1Param();

sim = Simulation(lensmodel, freqs, f);

designer = PhilDesigner(lensmodel, sim);

designer.DesignDielectricGrin();
