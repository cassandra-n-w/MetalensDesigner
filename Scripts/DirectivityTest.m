lambda = 5;
frequency = Constants.c/lambda;
freqs = frequency*linspace(0.5,2,100);
lambdas = Constants.c./freqs;
diam = lambda*20;
radius = diam/2;
k = 2*pi / lambda;

lens = Lens(diam, 0, diam*3, lambda/5, 0, diam);
sim = Simulation(lens, frequency, frequency, 2);

sim.initialize_E_field(@(x,y,f) Simulation.UniformCircle(x,y,radius));
sim.calculate_power();
Simulation.calc_power2(sim.E_current(:,:,1));
gain = sim.FarFieldGainPatternFromIdx(1);
maxgain = sim.FarFieldMaxGain(1);
maxgaindb = pow2db(maxgain)
expectedgain = (2*pi*radius/lambda)^2;
expectedgaindb = pow2db(expectedgain)

maxgain./expectedgain;


%plot()


