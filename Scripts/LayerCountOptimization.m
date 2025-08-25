f = 480e9;

lambda = Constants.c / f;

n_poly = 1.87;
n_poly = sqrt(2.2);

n_max_hypothetical = n_poly^2 + sqrt(n_poly^4 - n_poly^2);

n_max = linspace(2, 15, 200);

g = lambda/5;
a = 0.003;

alpha = 0.95;
beta = 0.07;

diameter = 130;

focal_length = 150;

OPDR = -diameter^2/(8 * focal_length);

thickness = (lambda / 2) ./ (n_poly + sqrt(n_max)) + OPDR ./ (n_poly^2 - n_max);

spacing = alpha * (2 * g / pi) * log(csc(pi * a / (2 * g))) ./ ((n_max.^2 / n_poly^2) - 1) - beta * lambda * a / g;

plot(n_max, thickness ./ spacing, n_max, thickness, n_max, spacing);
ylim([0, 500])

