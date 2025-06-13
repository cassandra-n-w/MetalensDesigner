

f = 54e9; % 480GHz central frequency
numfreqs=100;
freqs = linspace(f*0.8, f*1.2, numfreqs);
%freqs = linspace(451e9, 499e9, numfreqs);

lambda = Constants.c/f;
%gridding size
%g=0.11992; % 0.11992mm ; this should be equal to lambda(500GHz, free space)/5
g = 0.2 * lambda;
%g=lambda*1;

epsr_polyimide = 3.4;

total_thickness = 4.05; %16mm total thickness

% diameter of active lens area and structural lens area, respectively
diameter = 1000;
struct_diam = 1000;

focal_length = 3000;

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

designer = PhilDesigner(lensmodel, sim, sqrt(epsr_polyimide), 6.5);




%%
%% Tests on index calculation
squaresize = linspace(0.01*g, 0.985*g, 500);

spacing = 2.5/150;
alpha = 1.15;
beta = 0.035;

index_adhitya = designer.GetIndexFromSpacing(squaresize, epsr_polyimide, spacing, alpha, beta);
%index_emprical = designer.GetIndexFromPhaseSlope(squaresize/g, spacing);
%plot(squaresize/g, index_adhitya, squaresize/g, index_emprical);
plot(squaresize, index_adhitya);

%%
disp("Fresnelify")
designer.ZoneBoundsEveryTwoPi();

%%
disp("Minimize lens thickness");
[min_thickness, AR_thickness, phase_thickness] = designer.MinimizeThickness();

%designer.DesignDielectricGrinDiscretizeZ(11);
disp("Designing Lens");
%designer.DesignDielectricGrinDiscretizeZ(150);
designer.DesignDielectricGrin();

%% save design's phase output for comparison purposes.
lensmodel.CalcSParam(freqs);
phase_ideal_design = angle(lensmodel.S_array(:,:,19));

%% Convert ideal dielectric to metamaterial design
disp("Converting Lens to Metamaterial Design");
polyimide_dielectric = Dielectric(3.50 * ones(size(freqs)), 0.015, freqs, 0.005, 0.015);
lossless_polyimide = Dielectric(3.50 * ones(size(freqs)), 0.0, freqs, 0.005, 0.015);
aluminum_metaatom = Element('..\\..\\HFSS Models\\Polyimide Wideband 3_4 200nm Aluminum g0p2 v99 touchstone');
aluminum_old = Element('..\\..\\HFSS Models\\Polyimide 3_4 Aluminum g0p2 v99 touchstone');
proto_tl = ProtoTL(aluminum_metaatom, lossless_polyimide);
tracespace = 0.003; % 3um converted to mm.

%%
disp("Not actually, skipping conversion");
%designer.ConvertLensToMetamaterial(proto_tl, tracespace)

%% Tests on index calculation
squaresize = linspace(0.01*g, 0.985*g, 500);

epsr_polyimide = 3.4;
spacing = 2.5/150;
alpha = 1.15;
beta = 0.035;

index_adhitya = designer.GetIndexFromSpacing(squaresize, epsr_polyimide, spacing, alpha, beta);
%index_emprical = designer.GetIndexFromPhaseSlope(squaresize/g, spacing);
%plot(squaresize/g, index_adhitya, squaresize/g, index_emprical);
plot(squaresize/g, index_adhitya);


%%
% test path length algorithms
disp("Test Algorithms");
dr = lambda/20;
r = 0:dr:(lensmodel.diameter/2);

pathNaive = designer.FindNaivePathLength(r);
pathIdealRaw = designer.PathLengthInterpRaw(r);
pathIdeal = designer.PathLengthInterp(r);

Rz = designer.idealCurvature();

plot(r, pathNaive, r, pathIdeal, r, pathIdealRaw);

legend(["Naive","Ideal","Raw"]);

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
scan_z = focal_length*linspace(0.0,1.3,15);

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

%% Plot phase output of lens surface (
%phase_indx = 
phase = angle(lensmodel.S_array(:,:,30));
imagesc(lensmodel.xvec, lensmodel.yvec, phase);
colorbar;

%% 
phase_size = size(phase);
halfsize1 = round(phase_size(1)/2);
halfsize2 = round(phase_size(2)/2);

phase1d = (phase(halfsize1, halfsize2 : end));

phase_ideal_1d = (phase_ideal_design(halfsize1, halfsize2 : end));

phase_error = unwrap(phase_ideal_1d) - unwrap(phase1d);
%phase_error = abs(phase_ideal_1d - phase1d); 
%flip_phase_error = 2*pi - phase_error;
%phase_error(phase_error > pi) = flip_phase_error(phase_error > pi);
plotx = 1:length(phase1d);

plot(plotx,rad2deg(phase1d), plotx, rad2deg(phase_ideal_1d), plotx, rad2deg(phase_error));
xlabel("Distance from lens center (grid units)");

%% Plot central metamaterial square sizes

if (designer.is_metamaterial)
    
    tl_collection = lensmodel.TL_collection;
    
    for i = 1:length(tl_collection)
        squaresize(i) = tl_collection(i).sizes(75);
    end
    
    plot(squaresize);
end

%% Plot all loss sources
gfreqs = freqs/1e9;
plot(gfreqs, strehl_ratio, gfreqs, ref, gfreqs, d, gfreqs, o);

legend(["Strehl ratio", "antireflective efficiency", "dielectric/conductive efficiency", "optical efficiency"]);
ylim([0,1]);
title("Strehl Ratio of Broadband Lens Design");
xlabel("Frequency (GHz)");
ylabel("Efficiency (unitless)");

[bwminfreq, bwmaxfreq] = Simulation.calc_bandwidth(strehl_ratio, 0.5, gfreqs);

%% Plot strehl only

plot(freqs/1e9, strehl_ratio);

%legend(["Strehl ratio", "antireflective efficiency", "dielectric/conductive efficiency", "optical efficiency"]);
ylim([0,1]);
title("Strehl Ratio of Simulated Broadband 480 GHz Lens Design");
xlabel("Frequency (GHz)");
ylabel("Strehl Ratio (dimensionless)");

%% Plot the side view of lens beam

E_size = size(E_save);
E_plot = squeeze(E_save(round(E_size(1)/2),:,:));

imagesc(scan_z, sim2.xvec, mag2db(abs(E_plot)));
colorbar;
max_sig = max(mag2db(abs(E_plot)), [], "all");
clim([max_sig - 18, max_sig]);

%% Plot the side-view of lens phase simulation

imagesc(scan_z, sim2.xvec, angle(E_plot));
ylim([-lensmodel.diameter/2, lensmodel.diameter/2]);
clim("auto");


%% Draw a cartoon of the lens structure
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

%% print out lens specifications

disp(strcat("Diameter: ", num2str(diameter), "mm"));
disp(strcat("Focal Length: ", num2str(focal_length), "mm"));
disp(strcat("Thickness: ", num2str(designer.total_thickness), "mm : AR ", num2str(AR_thickness), "mm + phase ", num2str(phase_thickness), "mm"));
micron_tracespace = true;
if (micron_tracespace)
    tracespace = 0.003;
    disp(strcat("Min Feature Size: ", num2str(tracespace * 1000), " micron"));
else
    mil_to_mm = 0.0254;
    tracespace = 3 * mil_to_mm;
    disp(strcat("Min Feature Size: ", num2str(tracespace / mil_to_mm), " mil"));
end

layerspacing = designer.GetMinLayerSpacing(alpha, beta, tracespace);
disp(strcat("Est. Layer Spacing: ", num2str(layerspacing), "mm"));
disp(strcat("Est. Layer Count: ", num2str(designer.total_thickness/layerspacing)));

disp(strcat("Bandwidth (50% Strehl): ", num2str(bwminfreq, '%.4g'), " to ", num2str(bwmaxfreq, '%.4g'), " GHz"));