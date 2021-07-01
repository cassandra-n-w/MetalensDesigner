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

scanpos = (-40:tempsim.dy:40);
r = sqrt((tempsim.x).^2 + (tempsim.y).^2);
mask = (r < (lensmodel.diameter/2));


gaussiantransf = tempsim.calc_gaussian_phase();
gaussianscan = sim_test(gaussiantransf, lensmodel, hornfunc, distance, scanpos, f0);

idealtransf = tempsim.calc_ideal_phase();
idealscan = sim_test(idealtransf, lensmodel, hornfunc, distance, scanpos, f0);

optimtransf_square = (lensmodel.PadSParams(zeros(tempsim.dims), f0));
optimsquarescan = sim_test(optimtransf_square, lensmodel, hornfunc, distance, scanpos-1.499, f0);

optimtransf_circle = optimtransf_square .* mask;
optimcirclescan = sim_test(optimtransf_circle, lensmodel, hornfunc, distance, scanpos-1.499, f0);

maxpower = pow2db(max(abs(optimcirclescan)));
maxpower = 0;

gaussianscan_db = pow2db(abs(gaussianscan)) - maxpower;
idealscan_db = pow2db(abs(idealscan)) - maxpower;
optimsquarescan_db = pow2db(abs(optimsquarescan)) - maxpower;
optimcirclescan_db = pow2db(abs(optimcirclescan)) - maxpower;

plot(scanpos, gaussianscan_db, scanpos, idealscan_db, scanpos, optimsquarescan_db, scanpos, optimcirclescan_db);
legend("ideal gaussian lens", "experimental","optimized lens (square)", "optimized lens (circle)");

function [scancoupling, sim] = sim_test(lens_transformation, lensmodel, hornfunc, distance, scanpos, f0)
    sim = Simulation(lensmodel, f0, f0);
    
    sim.initialize_E_field(hornfunc);

    sim.propagate(distance);
    
    sim.transform(lens_transformation);
    
    sim.propagate(distance);
    
    scancoupling = sim.receiver_scan(hornfunc, 0, scanpos);
end

