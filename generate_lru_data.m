function [lambda_vec, T_vec, c_vec] = generate_lru_data(seed)
%GENERATE_LRU_DATA  Generate LRU parameters for the spare parts problem.
%
%  [lambda_vec, T_vec, c_vec] = generate_lru_data(seed)
%
%  Returns randomised LRU parameters centred on representative base values
%  for a 9-component system.  Randomisation ensures variability across
%  experiments while remaining physically plausible.
%
%  Inputs
%    seed       — integer seed for reproducibility (default: 42)
%
%  Outputs
%    lambda_vec — failure rates  [failures / time unit],  1×9
%    T_vec      — mean repair / resupply times [time units], 1×9  (integers)
%    c_vec      — unit spare part costs [currency units],   1×9  (integers)
%
%  Base values represent a plausible industrial system.  Each parameter is
%  perturbed by ±20 % uniform noise around the base value.

if nargin < 1
    seed = 42;
end

% --- Base (nominal) values ---
lambda_base = 1/1000 * [50, 40, 45, 50, 25, 48, 60, 35, 15];
T_base      =          [ 6,  8, 14, 25, 12, 18, 33,  8, 12];
c_base      =          [12, 14, 21, 20, 11, 45, 75, 30, 22];

n = length(lambda_base);

% --- Apply ±20 % uniform perturbation ---
rng(seed);
lambda_vec = lambda_base .* (1 + 0.2 * rand(1, n));
T_vec      = round(T_base  .* (1 + 0.2 * rand(1, n)));
c_vec      = round(c_base  .* (1 + 0.2 * rand(1, n)));

end
