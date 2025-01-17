freqs = 1e9 * logspace(0,3,500);
freqs2 = 1e9 * logspace(0,3,500);
lengths = logspace(0,3,50000);

dielectric = Dielectric(35, 0, freqs, 1, 1000);

abcd = dielectric.abcd_interp(freqs2, 100);

s = abcd2s(abcd, Constants.Z0);

total = abs(s(1,1,:)).^2 + abs(s(2,1,:)).^2;
total = squeeze(total);
hold off;
plot(freqs2, total);
hold on;
%% 


abcd = dielectric.abcd_interp(freqs(250), lengths);
abcd = reshape(abcd, 2,2,[]);
s = abcd2s(abcd, Constants.Z0);

total = abs(s(1,1,:)).^2 + abs(s(2,1,:)).^2;
total = squeeze(total);

hold off;
plot(lengths, total);
hold on;
plot(lengths, squeeze(abs(s(1,1,:)).^2));
plot(lengths, squeeze(abs(s(2,1,:)).^2));
