lambda = 273;
frequency = Constants.c/lambda;
freqs = frequency*linspace(0.5,2,100);
lambdas = Constants.c./freqs;
diam = lambda*20;
radius = diam/2;
k = 2*pi / lambda;

lens = Lens(diam, 0, diam*3, lambda/5, 0, diam);
sim = Simulation(lens, frequency, frequency, 10);

sim.initialize_E_field(@(x,y,f) Simulation.UniformCircle(x,y,radius));
sim.calculate_power();

Simulation.calc_power2(sim.E_current(:,:,1));
[gain, theta, phi] = sim.FarFieldGainPatternFromIdx(1);
maxgain = sim.FarFieldMaxGain(1);
maxgaindb = pow2db(maxgain)
expectedgain = (2*pi*radius/lambda)^2;
expectedgaindb = pow2db(expectedgain)

maxgain./expectedgain;

gainsize = size(gain);


[mintheta, thetaidx] = min(theta, [], 1);
[mintheta, minrow] = min(mintheta);
selectPos = ((phi(:,minrow) > 0) - 0.5) * 2;

assert(all(abs(selectPos) == ones(size(selectPos))));

plot(rad2deg(theta(:, minrow)).*selectPos, pow2db(gain(:,minrow)));
ylim([0,maxgaindb]);

%%

efield = sim.E_current(:,:,1);
Simulation.calc_power(efield)
Nx = gainsize(1);
Ny = gainsize(2);
testfft = fftshift(sqrt(1/(Nx*Ny)) * fft2((efield)));
Simulation.calc_power(testfft)


