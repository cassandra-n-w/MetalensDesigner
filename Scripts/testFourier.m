% This file is part of MetalensDesigner.
% 
% MetalensDesigner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% 
% MetalensDesigner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with MetalensDesigner. If not, see <https://www.gnu.org/licenses/>. 


lambda = 1.5; % cm
radius = 25.4/2;
dx = 0.15;
dy = 0.15;
k = 2*pi/lambda;

dims = [1000,1000];
xbounds = [-(dims(1)-1)/2, (dims(1)-1)/2];
ybounds = [-(dims(2)-1)/2, (dims(2)-1)/2];

[x, y] = meshgrid(linspace(xbounds(1), xbounds(2), dims(1)), linspace(ybounds(1), ybounds(2), dims(2)));

x = x * dx;
y = y * dy;

sigma = 1;

r = sqrt(x.^2 + y.^2);
mag = r < radius;
Rz = 10.67; % cm: radius of curvature
inctheta = 0 * pi/180; % angle of incidence in radians
incphase = exp(1.0i * 2*pi * sin(inctheta) * x / lambda);
phasetrans = exp(-1.0i * pi * r.^2 / (lambda * Rz));

E = mag .* incphase .* phasetrans;
% exp( -(x.^2 + y.^2)/(2 * sigma^2) ); % gaussian beam

[A, kx, ky] = FourierEtoA(E, dx, dy);

zmax = 13;
z = 0:0.1:zmax;

Ez = FourierAtoE(A, z, kx, ky, k);

section = Ez(ceil(dims(2)/2), :, :);

sectionreshape = [];
sectionreshape(:, :) = section(1, :, :);

powerdb = 20*log10(abs(sectionreshape));

maxpowerz = max(powerdb, [], 1);

relpowerdb = powerdb - repmat(maxpowerz, dims(1), 1);
phase = angle(sectionreshape);

subplot(1,2,1);
imagesc([min(z) max(z)], [min(x) max(x)], powerdb, [-30, 25]);
title('Power of focused beam');
%subtitle('');
xlabel('Distance from lens (cm)');
ylabel('y (cm)');
c = colorbar;
c.Label.String = 'Power (dB relative to inc. plane wave)';

subplot(1,2,2);
imagesc([min(z) max(z)], [min(x) max(x)], phase);
title('Phase of focused beam');
xlabel('Distance from lens (cm)');
ylabel('y (cm)');
c = colorbar;
c.Label.String = 'Phase (radians)';
