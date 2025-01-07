freqs = 1e9 * logspace(0,3,500);
freqs2 = 1e9 * logspace(3,3,5000);
lengths = logspace(2.1,2.9,50000);

dielectric = Dielectric(35, 0, freqs, 100, 1000);

abcd = dielectric.abcd_interp(freqs2, 100);

s = abcd2s(abcd, Constants.Z0);

total = abs(s(1,1,:)).^2 + abs(s(2,1,:)).^2;
total = squeeze(total);

plot(freqs2, total);

abcd = dielectric.abcd_interp(freqs(10), lengths);
abcd = reshape(abcd, 2,2,[]);
s = abcd2s(abcd, Constants.Z0);

total = abs(s(1,1,:)).^2 + abs(s(2,1,:)).^2;
total = squeeze(total);

hold off;
plot(lengths, total);
hold on;
plot(lengths, squeeze(abs(s(1,1,:)).^2));
plot(lengths, squeeze(abs(s(2,1,:)).^2));
