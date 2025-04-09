
f = 480e9; % 480GHz central frequency
numfreqs=50;
freqs = linspace(f*0.4, f*2, numfreqs);

lambda = Constants.c/f;
%gridding size
g=0.11992; % 0.11992mm ; this should be equal to lambda(500GHz, free space)/5
%g=lambda*1;

total_thickness = 2.5; %16mm total thickness

% diameter of active lens area and structural lens area, respectively
diameter = 130;
struct_diam = 130;

focal_length = 150;

% total_thickness = 0.3; %16mm total thickness
% 
% % diameter of active lens area and structural lens area, respectively
% diameter = 11.25;
% struct_diam = 11.25;
% 
% focal_length = 15;

% 0 is used as a placerholder for the prototype TL in this lens, since each
% TL has a different prototype
lensmodel = Lens(diameter, total_thickness,  focal_length, g, 0, struct_diam);
lensmodel.DeclareTL_1Param();

sim = Simulation(lensmodel, freqs, f);

designer = PhilDesigner(lensmodel, sim);

%%
squaresize = linspace(0, g, 500);

epsr_polyimide = 3.4;
spacing = 0.2;
alpha = 1.15;
beta = 0.07;

index = designer.GetIndexFromSpacing(squaresize, epsr_polyimide, spacing, alpha, beta);

plot(squaresize/g, index);

%%


%%
%designer.DesignDielectricGrinDiscretizeZ(11);
disp("Designing Lens");
designer.DesignDielectricGrinDiscretizeZ(150);
%designer.DesignDielectricGrin();

%%
disp("Converting Lens to Metamaterial Design");
polyimide_dielectric = Dielectric(3.50 * ones(size(freqs)), 0.015, freqs, 0.005, 0.015);
lossless_polyimide = Dielectric(3.50 * ones(size(freqs)), 0.0, freqs, 0.005, 0.015);
aluminum_metaatom = Element('..\\..\\HFSS Models\\Polyimide Wideband 3_4 200nm Aluminum g0p2 v99 touchstone');
proto_tl = ProtoTL(aluminum_metaatom, lossless_polyimide);
tracespace = 0.003; % 3um converted to mm.
%designer.ConvertLensToMetamaterial(proto_tl, tracespace)

%%
% test path length algorithms
disp("Test Algorithms");
dr = lambda/20;
r = 0:dr:(lensmodel.diameter/2);

pathNaive = designer.FindNaivePathLength(r);
pathIdeal = designer.FindIdealPathLength(r);

Rz = designer.idealCurvature();

plot(r, pathNaive, r, pathIdeal, r, r.^2/(2*Rz));

legend(["Naive","Ideal","Curvature"]);

% plots suggest not perfectly quadratic

%% 
disp("Simulate");
sim = Simulation(lensmodel, f, f);
disp("Calc S Params");
lensmodel.CalcSParam(f);
% calculate the gaussian coupling of the lens

disp("Initialize E Field");
efield = @(x,y,f) (sqrt(x.^2 + y.^2) < (lensmodel.diameter/2)) * exp(0);
% initialize incoming plane wave
sim.initialize_E_field(efield);

disp("Lens Transform");
% transform through ideal lens pattern
sim.lensTransform();

% propagate
disp("Propagate");
scan_z = focal_length*linspace(0.0,1.3,300);

sim.propagate(scan_z);
disp("Save field");
E_save = sim.E_saved;

%strehl_ratio = sim.StrehlRatio();
disp("Calc S Params");
sim2 = Simulation(lensmodel, freqs, f);
lensmodel.CalcSParam(freqs);
disp("Calc Losses");
[ref, d, o, strehl_ratio] = sim2.CalcLoss();
disp("Done Simulating");

%%

passivity = sqrt(abs(lensmodel.S_array(:,:,round(numfreqs/2))).^2 + abs(lensmodel.S11_array(:,:,round(numfreqs/2)).^2));

imagesc(lensmodel.xvec, lensmodel.yvec, passivity);
colorbar;

%%
phase = angle(lensmodel.S_array(:,:,40));
imagesc(lensmodel.xvec, lensmodel.yvec, phase);
colorbar;

%%
gfreqs = freqs/1e9;
plot(gfreqs, strehl_ratio, gfreqs, ref, gfreqs, d, gfreqs, o);

legend(["Strehl ratio", "antireflective efficiency", "dielectric/conductive efficiency", "optical efficiency"]);
ylim([0,1]);
title("Strehl Ratio of Broadband Lens Design");
xlabel("Frequency (GHz)");
ylabel("Efficiency (unitless)");

%%

plot(freqs/1e9, strehl_ratio);

%legend(["Strehl ratio", "antireflective efficiency", "dielectric/conductive efficiency", "optical efficiency"]);
ylim([0,1]);
title("Strehl Ratio of Simulated Broadband 480 GHz Lens Design");
xlabel("Frequency (GHz)");
ylabel("Strehl Ratio (dimensionless)");

%%

E_size = size(E_save);
E_plot = squeeze(E_save(round(E_size(1)/2),:,:));

imagesc(scan_z, sim2.xvec, mag2db(abs(E_plot)));
colorbar;
max_sig = max(mag2db(abs(E_plot)), [], "all");
clim([max_sig - 18, max_sig]);

%%

imagesc(scan_z, sim2.xvec, angle(E_plot));
clim("auto");


%%
clf;
hold on;
rectangle('Position', [-diameter/2, -total_thickness/2, diameter, total_thickness], 'FaceColor',[1,1,0]);
radius = linspace(0, diameter/2, 100);
pathlength = designer.PathLengthInterp(radius);
[n_of_r, t_of_r, ar_thick] = designer.DesignIdealDielectric(pathlength); 
plot(radius, t_of_r/2, 'r--', -radius, t_of_r/2, 'r--',radius, -t_of_r/2,'r--', -radius, -t_of_r/2,'r--', "LineWidth", 3);
pattern_dims = [7, 5];
thickness = 0.005;

for i = 1:(pattern_dims(1))
    for j = 1:(pattern_dims(2))
        r = (j-0.5) * diameter/2 / (pattern_dims(2));
        
        h = (i-0.5) * total_thickness/2 / (pattern_dims(1));
        pathlength = designer.PathLengthInterp(r);
        [n_of_r, t_of_r, ar_thick] = designer.DesignIdealDielectric(pathlength);
        if h > t_of_r/2
            n = sqrt(n_of_r);
        else
            n = n_of_r;
        end
        n = n/50;
        width = n*7.5;
        rectangle('Position',[r-width/2, h-thickness/2, width, thickness], 'FaceColor',[0.5,0.5,0.5]);
        rectangle('Position',[-r-width/2, h-thickness/2, width, thickness], 'FaceColor',[0.5,0.5,0.5]);
        rectangle('Position',[r-width/2, -h-thickness/2, width, thickness], 'FaceColor',[0.5,0.5,0.5]);
        rectangle('Position',[-r-width/2, -h-thickness/2, width, thickness], 'FaceColor',[0.5,0.5,0.5]);
    end
end

axis off;
hold off;