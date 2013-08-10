% Plot oscillators used in test.

% parameters
K_k = [1., 10., 20., 40., 80., 160.]; % spring constants for each state
O_k = [0., 0.5, 1.0, 1.5, 2.0, 2.5]; % offsets for spring constants
beta = 1.0 % inverse temperature

% determine number of oscillators
K = length(K_k);

% compute standard deviations for oscillators
sigma_k = (beta * K_k).^-0.5

% compute min and max extents
nsigma = 4.0 % number of standard deviations to go out.
xmin = O_k(1)
xmax = O_k(1);
for k = 1:K
  xmin = min(xmin, O_k(k) - nsigma * sigma_k(k));
  xmax = max(xmax, O_k(k) + nsigma * sigma_k(k));
end

% generate plot
npoints = 1000;
x = linspace(xmin, xmax, npoints);
y_kn = zeros(K, npoints);
for k = 1:K
  y_kn(k,:) = 1/(sqrt(2.0 * pi) * sigma_k(k)) * exp(-(x - O_k(k)).^2 / (2 * sigma_k(k)));
end 
plot(x, y_kn);
xlabel('x');
ylabel('p(x)');

% Write file
filename = 'oscillators.eps';
print('-depsc', filename);
unix(sprintf('epstopdf %s', filename));

