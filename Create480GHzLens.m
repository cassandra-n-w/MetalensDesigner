f = 480e9; % 480GHz central frequency

freqs = linspace(f*0.9, f*1.1, 30);

el = Element('C:\\Users\\Cassie\\My Drive\\Dropbox\\Documents\\My Work\\Research\\metamaterial_lens\\code\\HFSS Models\\Polyimide 3_4 g0p2 v99 touchstone');
di = Dielectric(3.50 * ones(size(freqs)), 0.015, freqs, 0.005, 0.015); %layer thickness 5 - 15um, dielectric 3.5, loss tangent 0.015

proto = ProtoTL(el, di);

%gridding size
g=0.11992; % 0.11992mm ; this should be equal to lambda(500GHz, free space)/5

layer_thickness = 0.010; %10 um layer thickness

diameter = 150;

focal_length = 150;

lensmodel = Lens(diameter, layer_thickness*ones(11, 1),  focal_length, g, proto);

sim = Simulation(lensmodel, freqs, f);

horn_waist = sim.calc_waist();

ideal = sim.calc_ideal_phase();

gaussapprox = sim.calc_gaussian_phase();


imagesc(sim.xvec, sim.yvec, 180+angle(ideal)*180/pi);
xlabel("X position (mm)");
ylabel("Y position (mm)");
colorbar;
title("Phase pattern over lens surface (degrees)");
xlim([-75, 75]);
ylim([-75, 75]);

