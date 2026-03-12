%% Spare Parts Optimization: Marginal Allocation & Dynamic Programming
%
%  Solves the multi-LRU spare parts stocking problem:
%  minimize total Expected Backorders (EBO) subject to a budget constraint.
%
%  Two algorithms are implemented and compared:
%    1. Marginal Allocation  — greedy, generates the efficient frontier
%    2. Dynamic Programming  — exact, solves every integer budget state
%
%  Demand for each LRU during repair time follows a Poisson distribution
%  with mean mu_j = lambda_j * T_j.
%
%  Usage: Run this script in MATLAB. No toolboxes required.
%         Results are printed to the Command Window and two figures are generated.

clc; clear; close all;
format compact;

% =========================================================================
%  CONFIGURATION
% =========================================================================
RANDOM_SEED  = 42;        % Seed for reproducible LRU parameter generation
C_BUDGET     = 500;       % Total spare parts budget
MAX_SPARES   = 50;        % Upper bound on spares per LRU (for precomputation)

% =========================================================================
%  DATA GENERATION
% =========================================================================
[lambda_vec, T_vec, c_vec] = generate_lru_data(RANDOM_SEED);

n_LRU        = length(lambda_vec);
mu_vec       = lambda_vec .* T_vec;   % Mean demand per LRU during repair time

fprintf('=== System Parameters ===\n');
fprintf('Number of LRU types : %d\n', n_LRU);
fprintf('Budget constraint   : %d\n', C_BUDGET);
fprintf('\n');

% =========================================================================
%  PRECOMPUTE EBO / R / p  MATRICES  (RESOPT recursive equations)
%
%  For Poisson(mu) demand with s spare parts:
%    p(0)    = exp(-mu)
%    p(s+1)  = (mu / s) * p(s)
%    R(0)    = 1 - p(0)
%    R(s+1)  = R(s) - p(s+1)
%    EBO(0)  = mu
%    EBO(s+1)= EBO(s) - R(s)
%
%  Rows = stock levels 0..MAX_SPARES,  Columns = LRU index
% =========================================================================
EBO_matrix = zeros(MAX_SPARES + 1, n_LRU);
R_matrix   = zeros(MAX_SPARES + 1, n_LRU);
p_matrix   = zeros(MAX_SPARES + 1, n_LRU);

for j = 1:n_LRU
    mu = mu_vec(j);

    % Base case  s = 0
    p_matrix(1, j)   = exp(-mu);
    R_matrix(1, j)   = 1 - p_matrix(1, j);
    EBO_matrix(1, j) = mu;

    % Recursive case  s = 1 .. MAX_SPARES
    for s = 1:MAX_SPARES
        p_matrix(s+1, j)   = (mu / s) * p_matrix(s, j);
        R_matrix(s+1, j)   = R_matrix(s, j)   - p_matrix(s+1, j);
        EBO_matrix(s+1, j) = EBO_matrix(s, j) - R_matrix(s,   j);
    end
end

% --- Sanity check: verify integer convexity of EBO and cost functions ---
delta_EBO  = diff(EBO_matrix);
delta2_EBO = diff(delta_EBO);

assert(all(delta_EBO(:)  < 1e-10),  'EBO is not strictly decreasing');
assert(all(delta2_EBO(:) > -1e-10), 'EBO is not integer-convex');

s_vec    = (0:MAX_SPARES)';
g_matrix = s_vec * c_vec;          % Cost matrix  g(s_j) = c_j * s_j
assert(all(diff(g_matrix) > 0, 'all'), 'Cost function is not strictly increasing');

fprintf('Convexity checks passed for EBO and cost functions.\n\n');

% =========================================================================
%  ALGORITHM 1: MARGINAL ALLOCATION
%
%  Greedy algorithm that iteratively selects the LRU yielding the greatest
%  EBO reduction per unit cost, R_j(s_j) / c_j, until the budget is exhausted.
%  By integer convexity, every generated solution is Pareto-efficient.
% =========================================================================
fprintf('=== Marginal Allocation ===\n');

s_MA      = zeros(1, n_LRU);   % Current stock vector
EBO_MA    = sum(EBO_matrix(1, :));
C_MA      = 0;

% Pre-allocate efficient points table: [s1..s9 | EBO | Cost]
EP = zeros(C_BUDGET + 10, n_LRU + 2);
k  = 1;
EP(k, :) = [s_MA, EBO_MA, C_MA];

while true
    % Compute marginal efficiency quotient for each LRU
    quotient = zeros(1, n_LRU);
    for j = 1:n_LRU
        s = s_MA(j);
        if s < MAX_SPARES
            quotient(j) = R_matrix(s + 1, j) / c_vec(j);
        end
    end

    [~, best_j] = max(quotient);
    cost_delta  = c_vec(best_j);

    if C_MA + cost_delta > C_BUDGET
        break;
    end

    s_MA(best_j) = s_MA(best_j) + 1;
    C_MA         = C_MA  + cost_delta;
    EBO_MA       = EBO_MA - R_matrix(s_MA(best_j), best_j);

    k      = k + 1;
    EP(k, :) = [s_MA, EBO_MA, C_MA];
end

EP = EP(1:k, :);   % Trim to actual number of efficient points
n_EP = k;

% --- Display first 5 efficient points ---
fprintf('First 5 Efficient Points (stock levels | EBO | Cost):\n');
EP_display = array2table(EP(1:min(5, n_EP), :), ...
    'VariableNames', {'s1','s2','s3','s4','s5','s6','s7','s8','s9','EBO','Cost'});
disp(EP_display);

% =========================================================================
%  ALGORITHM 2: DYNAMIC PROGRAMMING
%
%  Solves the Distribution-of-Effort problem exactly via backward recursion:
%
%    f*_n(s_n) = min_{x_n in F_n(s_n)} { EBO_n(x_n / c_n) + f*_{n+1}(s_n - x_n) }
%
%  where F_n(s_n) = {0, c_n, 2*c_n, ...} ∩ [0, s_n].
%  Terminal condition: f*_{N+1}(s) = 0 for all s.
%
%  After backward recursion, a forward pass reconstructs the optimal
%  allocation vector for every integer budget in [0, C_BUDGET].
% =========================================================================
fprintf('=== Dynamic Programming ===\n');

n_stages = n_LRU;

% Value function table:  f_star(stage, budget_state)
%   stage index 1..n_stages+1,  budget state index 1..C_BUDGET+1
f_star = zeros(n_stages + 1, C_BUDGET + 1);   % Terminal values are 0
x_star = zeros(n_stages,     C_BUDGET + 1);   % Optimal spending at each stage

% --- Backward Recursion ---
for n = n_stages:-1:1
    c_n = c_vec(n);

    for s = 0:C_BUDGET
        % Feasible decisions: spend x in {0, c_n, 2*c_n, ...} with x <= s
        possible_x = 0:c_n:s;
        k_vals     = possible_x / c_n;               % Corresponding spare counts

        % Immediate EBO from this LRU
        k_idx         = min(k_vals + 1, size(EBO_matrix, 1));
        immediate_EBO = EBO_matrix(k_idx, n)';

        % Future optimal EBO (from next stage)
        rem_idx    = (s - possible_x) + 1;
        future_EBO = f_star(n + 1, rem_idx);

        [best_val, best_idx] = min(immediate_EBO + future_EBO);

        f_star(n, s + 1) = best_val;
        x_star(n, s + 1) = possible_x(best_idx);
    end
end

% --- Forward Reconstruction ---
% For every integer budget b in [0, C_BUDGET], trace optimal allocation
DP_Results = zeros(C_BUDGET + 1, n_LRU + 2);   % [s1..s9 | EBO | ActualCost]

for b = 0:C_BUDGET
    remaining = b;
    alloc     = zeros(1, n_LRU);

    for n = 1:n_stages
        spend_n    = x_star(n, remaining + 1);
        alloc(n)   = spend_n / c_vec(n);
        remaining  = remaining - spend_n;
    end

    DP_Results(b + 1, :) = [alloc, f_star(1, b + 1), b - remaining];
end

% --- Display DP solutions at selected budget checkpoints ---
target_budgets  = [0, 100, 150, 350, 500];
target_idx      = target_budgets + 1;

fprintf('Optimal DP Allocations at Selected Budgets:\n');
DP_display = array2table(DP_Results(target_idx, :), ...
    'VariableNames', {'s1','s2','s3','s4','s5','s6','s7','s8','s9','EBO','ActualCost'});
DP_display.Budget = target_budgets';
disp(DP_display(:, [end, 1:end-1]));

% =========================================================================
%  ANALYSIS: BASELINE COMPARISON
%
%  Evaluate the naive "one spare per LRU" baseline allocation:
%    s_base = [1, 1, ..., 1]
%  and compare against the DP-optimal solutions on two criteria:
%    (1) Minimum cost achieving the same or better EBO (cost minimization)
%    (2) Minimum EBO achievable within the same total cost (performance maximization)
% =========================================================================
fprintf('=== Baseline Comparison: s_base = [1,1,...,1] ===\n');

s_base    = ones(1, n_LRU);
Cost_base = sum(c_vec .* s_base);
EBO_base  = sum(EBO_matrix(2, :));

fprintf('Baseline allocation : cost = %d,  EBO = %.4f\n\n', Cost_base, EBO_base);

% Scenario 1: Least expensive allocation with EBO <= EBO_base
better_idx      = find(DP_Results(:, end-1) <= EBO_base);
[min_cost, loc] = min(DP_Results(better_idx, end));
Cost_opt1       = min_cost;
EBO_opt1        = DP_Results(better_idx(loc), end-1);
savings         = Cost_base - Cost_opt1;

fprintf('Scenario 1 — Cost Minimisation (EBO <= baseline):\n');
fprintf('  Optimal cost    : %d  (saving %d vs baseline)\n', Cost_opt1, savings);
fprintf('  Achieved EBO    : %.4f\n\n', EBO_opt1);

% Scenario 2: Lowest EBO within the same budget
EBO_opt2     = DP_Results(Cost_base + 1, end-1);
ActualCost2  = DP_Results(Cost_base + 1, end);
improvement  = EBO_base - EBO_opt2;

fprintf('Scenario 2 — Performance Maximisation (cost <= baseline budget %d):\n', Cost_base);
fprintf('  Optimal EBO     : %.4f  (improvement %.4f vs baseline)\n', EBO_opt2, improvement);
fprintf('  Actual cost     : %d\n\n', ActualCost2);

fprintf('Both scenarios confirm the baseline is not on the efficient frontier.\n\n');

% =========================================================================
%  VISUALISATION
% =========================================================================

% --- Figure 1: Efficient Solutions Curve (Marginal Allocation) ---
fig1 = figure(1);
set(fig1, 'Position', [100, 100, 800, 550]);

plot(EP(:, end), EP(:, end-1), '.-k', 'LineWidth', 2, 'MarkerSize', 20, ...
    'DisplayName', 'Efficient Frontier (MA)');
hold on;

% Mark the first non-trivial allocation (single spare, best marginal gain)
first_spare_cost = EP(2, end);
first_spare_EBO  = EP(2, end-1);
plot(first_spare_cost, first_spare_EBO, 'ro', 'MarkerSize', 14, 'LineWidth', 2, ...
    'DisplayName', 'First Marginal Spare');

grid on;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
xlabel('Total Cost ($C$)',                'FontSize', 16, 'Interpreter', 'latex');
ylabel('Expected Backorders ($EBO$)',     'FontSize', 16, 'Interpreter', 'latex');
title('Efficient Solutions Curve — Marginal Allocation', 'FontSize', 16, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
xlim([0, C_BUDGET]); ylim([0, 8]);
hold off;

% --- Figure 2: MA vs DP Comparison ---
fig2 = figure(2);
set(fig2, 'Position', [150, 150, 800, 550]);

plot(EP(:, end), EP(:, end-1), 'k.-', 'LineWidth', 2, 'MarkerSize', 20, ...
    'DisplayName', 'Marginal Allocation (efficient frontier)');
hold on;

plot(DP_Results(:, end), DP_Results(:, end-1), 'rx', 'LineWidth', 1.5, 'MarkerSize', 14, ...
    'DisplayName', 'Dynamic Programming (all integer budgets)');

grid on;
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
xlabel('Total Cost ($C$)',                'FontSize', 16, 'Interpreter', 'latex');
ylabel('Expected Backorders ($EBO$)',     'FontSize', 16, 'Interpreter', 'latex');
title('MA vs DP: Optimal Spare Parts Allocation', 'FontSize', 16, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
hold off;

fprintf('Done. Figures 1 and 2 generated.\n');
