
% bond life time generation

function [lifetimes_power, lifetimes_exponential] = distribution_bond_lifetime(alpha, lambda, num_samples)

% Define the equation to solve for t for the truncated power law
inverse_cdf_power_law = @(U, alpha, lambda) fsolve(@(t) t^(-alpha) * exp(-lambda * t) - (1 - U), 1.0);

% Generate random samples for the truncated power law
U = rand(num_samples, 1);
lifetimes_power = arrayfun(@(u) inverse_cdf_power_law(u, alpha, lambda), U);

% Generate random samples for the exponential distribution
lifetimes_exponential = -log(rand(num_samples, 1)) / lambda;

end