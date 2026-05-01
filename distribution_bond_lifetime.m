
% bond life time generation

function [bond_lifetimes] = distribution_bond_lifetime(Num_conformation, conformation_diffusivity, num_samples);

N = Num_conformation;     
C = conformation_diffusivity; 

l = 1:1:100;
t = 1:0.1:850;

sum = 0;
sumP = 0;

for i = 1:length(l)
    % S(t) survival probability calculation
    n0 = 1;
    a0 = 1/(2 * i -1);
    a1 = sin(((2 * i - 1) * pi / (2 * N)) * n0);
    b = C * ((2 * i - 1) * pi/(2*N))^2;
    a2 = exp(- b .* t);
    d = a0 * a1 * a2;
    sum = sum + d;

    % bond lifetime distribribution function (-ds/dt)
    pl1 = (2 * i - 1) * sin(((2 * i - 1) * pi / (2 * N)));
    pl2 = C *  ((2 * i - 1) * pi/(2*N))^2;
    pl3 = pl1 * exp(- pl2 .* t);
    sumP = sumP + pl3;
end

Sur_Probability = (4/pi) * sum;
lifetime_dstr_function = (C * pi)/N^2 * sumP;

[xData, yData] = prepareCurveData( t, Sur_Probability );
ft = fittype( 'x^a*exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.8 0.08];
[fitresult, gof] = fit( xData, yData, ft, opts );
a = fitresult.a;
b = fitresult.b;

inverse_cdf_power_law = @(U, a, b) fsolve(@(t) t^(a) * exp(-b * t) - (1 - U), 1.0);
U = rand(num_samples, 1);
bond_lifetimes = arrayfun(@(u) inverse_cdf_power_law(u, a, b), U);



%%% old set

% % Define the equation to solve for t for the truncated power law
% inverse_cdf_power_law = @(U, alpha, lambda) fsolve(@(t) t^(-alpha) * exp(-lambda * t) - (1 - U), 1.0);
% 
% % Generate random samples for the truncated power law
% U = rand(num_samples, 1);
% lifetimes_power = arrayfun(@(u) inverse_cdf_power_law(u, alpha, lambda), U);
% 
% % Generate random samples for the exponential distribution
% lifetimes_exponential = -log(rand(num_samples, 1)) / lambda;

end