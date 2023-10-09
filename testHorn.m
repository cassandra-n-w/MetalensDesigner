% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


lambda = 1.499; % cm
radius = 25.4/2;
dx = 0.15;
dy = 0.15;
k = 2*pi/lambda;

% simulation granularity settings
dims = [1000,1000];
xbounds = [-(dims(1)-1)/2, (dims(1)-1)/2];
ybounds = [-(dims(2)-1)/2, (dims(2)-1)/2];

[x, y] = meshgrid(linspace(xbounds(1), xbounds(2), dims(1)), linspace(ybounds(1), ybounds(2), dims(2)));

x = x * dx;
y = y * dy;

sigma = 1;

% lens model, sort of
r = sqrt(x.^2 + y.^2);
mag = r < radius;
Rz = 10.67; % cm: radius of curvature
inctheta = 0 * pi/180; % angle of incidence in radians
incphase = exp(1.0i * 2*pi * sin(inctheta) * x / lambda);
phasetrans = exp(-1.0i * pi * r.^2 / (lambda * Rz));
E_lens = mag .* incphase .* phasetrans;

% horn model
cm_per_in = 2.54;
a1 = 1.125 * cm_per_in;
b1 = 0.800 * cm_per_in;
p = 1.31;
a = 1.067;
b = 0.432;

E_horn = HornField(a1, b1, p, a, b, k, x, y);

E_gaussian = exp( -(x.^2 + y.^2)/(2 * sigma^2) ); % gaussian beam


E = E_horn;
[A, kx, ky] = FourierEtoA(E, dx, dy);

zmax = 18.5;
z = 0:0.2:zmax;

Ez1 = FourierAtoE(A, z, kx, ky, k);

E2 = Ez1(:,:,end) .* E_lens;
[A, kx, ky] = FourierEtoA(E2, dx, dy);
zmax = 18.5;
z_tot = -z;
z = 0:0.2:zmax;
z_tot = [z_tot, z];
Ez2 = FourierAtoE(A, z, kx, ky, k);

Ez = cat(3, Ez1, Ez2);

xsection = Ez(ceil(dims(2)/2), :, :);
xsectionreshape = [];
xsectionreshape(:, :) = xsection(1, :, :);

ysection = Ez(:, ceil(dims(1)/2), :);
ysectionreshape = [];
ysectionreshape(:, :) = ysection(:, 1, :);

powerdb = 20*log10(abs(ysectionreshape));
phase = angle(ysectionreshape);

maxpowerz = max(powerdb, [], 1);
relpowerdb = powerdb - repmat(maxpowerz, dims(1), 1);

subplot(1,2,1);
imagesc([min(z_tot) max(z_tot)], [min(x) max(x)], relpowerdb, [-30,0]);
ylim([-20, 20]);
title('Power of focused beam');
%subtitle('');
xlabel('Distance from lens (cm)');
ylabel('y (cm)');
c = colorbar;
c.Label.String = 'Power (dB relative beam center)';

subplot(1,2,2);
imagesc([min(z_tot) max(z_tot)], [min(x) max(x)], phase);
ylim([-20,20]);
title('Phase of focused beam');
xlabel('Distance from lens (cm)');
ylabel('y (cm)');
c = colorbar;
c.Label.String = 'Phase (radians)';
