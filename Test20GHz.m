% designed to be run *after* running Load20GHzModel.m

mm_per_in = 25.4;
a1 = 1.125 * mm_per_in;
b1 = 0.800 * mm_per_in;
p = 13.1;
a = 10.67;
b = 4.32;

f0 = 20e9;

distance = 185;


myhorn = HornClass(a1, b1, p, a, b);
hornfunc = @(x,y,f) EField(myhorn, x,y,f);
%hornfunc = @myhorn.EField;%(xin,yin,fin);
hornfunc(0,0,20e9);

tempsim = Simulation(lensmodel, f0, f0);


scanpos = (-15*tempsim.dy:tempsim.dy:15*tempsim.dy);
r = sqrt((tempsim.x).^2 + (tempsim.y).^2);
mask = (r < (lensmodel.diameter/2));

% according to quasioptics by goldsmith, 7.4: Beam Radius, 
% a tapered aperture (such as our horn) has best coupling to the gaussian
% mode with w0/a = 0.35;
w0 = a1 .* 0.35;
gaussfunc = @(x,y,f) exp(-(x.^2 + y.^2) / w0^2);
tempsim.initialize_E_field(gaussfunc);
ideal_coupling = mag2db(abs(tempsim.receiver_scan(gaussfunc, 0, scanpos))) - 4.4;

gaussiantransf = tempsim.calc_gaussian_phase();
gaussianscan = sim_test(gaussiantransf, lensmodel, hornfunc, distance, scanpos, f0);

idealtransf = tempsim.calc_ideal_phase();
idealscan = sim_test(idealtransf, lensmodel, hornfunc, distance, scanpos, f0);

optimtransf_square = (lensmodel.PadSParams(zeros(tempsim.dims), f0));
optimsquarescan = sim_test(optimtransf_square, lensmodel, hornfunc, distance, scanpos-1.499, f0);

optimtransf_circle = optimtransf_square .* mask;
optimcirclescan = sim_test(optimtransf_circle, lensmodel, hornfunc, distance, scanpos-1.499, f0);

maxpower = mag2db(max(abs(optimcirclescan)));
maxpower = 0;

gaussianscan_db = mag2db(abs(gaussianscan)) - maxpower - 0.6;
idealscan_db = mag2db(abs(idealscan)) - maxpower;
optimsquarescan_db = mag2db(abs(optimsquarescan)) - maxpower;
optimcirclescan_db = mag2db(abs(optimcirclescan)) - maxpower;

VerticalScanhand1 = readtable('C:\\Users\\Cassie\\Google Drive\\Dropbox\\Documents\\My Work\\Research\\metamaterial_lens\\Measurements\\VerticalScan_hand.xlsx');

height = VerticalScanhand1.ReceiverHeight_mm__1;

powerdB = VerticalScanhand1.Power__dBm__3;

powerdB = powerdB - max(powerdB);

powerlin = 10 .^ (powerdB./10);

% calculate "center of mass" of the len's measurement
com = sum(height .* powerlin) / sum(powerlin);
height = height - com;

fitgauss = fit(height, powerlin, 'gauss1');
height = height - fitgauss.b1;
powerdB = powerdB - pow2db(fitgauss.a1) - 4.4;

p = plot(scanpos, gaussianscan_db, scanpos, idealscan_db, scanpos, optimsquarescan_db, scanpos, optimcirclescan_db, height, powerdB, scanpos, ideal_coupling);
legend("ideal gaussian lens", "experimental","optimized lens (square)", "optimized lens (circle)", "measured", "ideal coupling");

p(5).Marker = 'o';
p(5).LineStyle = 'none';
p(2).LineStyle = 'none';

function [scancoupling, sim] = sim_test(lens_transformation, lensmodel, hornfunc, distance, scanpos, f0)
    sim = Simulation(lensmodel, f0, f0);
    
    sim.initialize_E_field(hornfunc);

    sim.propagate(distance);
    
    sim.transform(lens_transformation);
    
    sim.propagate(distance);
    
    scancoupling = sim.receiver_scan(hornfunc, 0, scanpos);
end

