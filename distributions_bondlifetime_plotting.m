
% bond life time generation

% Parameters
alpha = 0.5;
lambda = 0.1;
num_samples = 1000; % Number of lifetimes to generate

% Define the equation to solve for t for the truncated power law
inverse_cdf_power_law = @(U, alpha, lambda) fsolve(@(t) t^(-alpha) * exp(-lambda * t) - (1 - U), 1.0);

% Generate random samples for the truncated power law
U = rand(num_samples, 1);
lifetimes_power_law = arrayfun(@(u) inverse_cdf_power_law(u, alpha, lambda), U);

% Generate random samples for the exponential distribution
lifetimes_exponential = -log(rand(num_samples, 1)) / lambda;


% Theoretical PDF for truncated power law
t_values = linspace(min(lifetimes_power_law), max(lifetimes_power_law), 1000);
pdf_values_power_law = (alpha * t_values.^(-alpha - 1) + lambda * t_values.^(-alpha)) .* exp(-lambda * t_values);

% Theoretical PDF for exponential distribution
pdf_values_exponential = lambda * exp(-lambda * t_values);

figure;
pbaspect([1 1 1]);
histogram(lifetimes_power_law, 50, 'Normalization', 'pdf','LineWidth', 1.5);
hold on;
plot(t_values, pdf_values_power_law, 'r', 'LineWidth', 2);
%title('Histogram of Generated Lifetimes - Truncated Power Law with Theoretical PDF');
xlabel('Bond lifetime (seconds)');
ylabel('PDF');
%legend('Generated Data', 'Theoretical PDF');
hold off;

figure;
pbaspect([1 1 1]);
histogram(lifetimes_exponential, 50, 'Normalization', 'pdf','LineWidth', 1.5);
hold on;
plot(t_values, pdf_values_exponential, 'r', 'LineWidth', 2);
%title('Histogram of Generated Lifetimes - Exponential Distribution with Theoretical PDF');
xlabel('Bond lifetime (seconds)');
ylabel('PDF');
%legend('Generated Data', 'Theoretical PDF');
hold off;