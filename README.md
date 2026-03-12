# Spare Parts Optimization

**Minimising Expected Backorders under a Budget Constraint using Marginal Allocation and Dynamic Programming**

---

## Overview

This project solves the classical **multi-item spare parts stocking problem**: given a system of $N$ Line Replaceable Units (LRUs), each with known failure rates, repair times, and unit costs, how should a fixed budget be allocated to spare parts in order to **minimise the total Expected Number of Backorders (EBO)**?

Two exact and well-established Operations Research algorithms are implemented from scratch in MATLAB and compared:

| Algorithm | Approach | Guarantees |
|---|---|---|
| **Marginal Allocation** | Greedy, iterative | Generates the full efficient frontier (convex hull) |
| **Dynamic Programming** | Backward recursion | Globally optimal solution for every integer budget |

Both algorithms exploit the **integer convexity** of the EBO function under Poisson demand, which is verified numerically at runtime.

---

## Problem Formulation

Each LRU $j$ experiences demand during its repair time that follows a **Poisson distribution** with mean $\mu_j = \lambda_j T_j$, where $\lambda_j$ is the failure rate and $T_j$ is the mean repair time.

The objective is:

$$\min_{\mathbf{s}} \sum_{j=1}^{N} \text{EBO}_j(s_j) \quad \text{subject to} \quad \sum_{j=1}^{N} c_j s_j \leq C_{\text{budget}}$$

where $s_j \in \mathbb{Z}_{\geq 0}$ is the number of spare units stocked for LRU $j$.

---

## Mathematical Background

### Recursive EBO Computation

The Expected Backorders and shortage probability for LRU $j$ with $s$ spare units are computed recursively (O(N × MAX_SPARES) precomputation, used by both algorithms):

$$p_j(0) = e^{-\mu_j}, \quad p_j(s+1) = \frac{\mu_j}{s} \cdot p_j(s)$$

$$R_j(0) = 1 - p_j(0), \quad R_j(s+1) = R_j(s) - p_j(s+1)$$

$$\text{EBO}_j(0) = \mu_j, \quad \text{EBO}_j(s+1) = \text{EBO}_j(s) - R_j(s)$$

These recursions confirm that $\text{EBO}_j(s)$ is **strictly decreasing** and **integer-convex** — conditions required for the Marginal Allocation algorithm to produce guaranteed efficient solutions.

### Marginal Allocation

The greedy algorithm selects the spare part with the highest **marginal efficiency quotient** at each step:

$$Q_j(s_j) = \frac{-\Delta \text{EBO}_j(s_j)}{\Delta C_j} = \frac{R_j(s_j)}{c_j}$$

This is iterated until the budget is exhausted, producing the **southwestern boundary of the feasible objective space** — the set of all Pareto-optimal (cost, EBO) pairs.

### Dynamic Programming

The problem is modelled as a **Distribution-of-Effort** problem. Each LRU corresponds to one stage. The backward recursion is:

$$f_{n}^{\ast}(s_{n}) = \min_{x_{n} \in F_{n}(s_{n})} \left\lbrace \text{EBO}_{n}\left(\frac{x_{n}}{c_{n}}\right) + f_{n+1}^{\ast}(s_{n} - x_{n}) \right\rbrace$$

where $F_{n}(s_{n}) = \lbrace 0,\ c_{n},\ 2c_{n},\ \ldots \rbrace \cap [0,\ s_{n}]$ is the set of feasible budget allocations at stage $n$.

A forward reconstruction pass then retrieves the optimal allocation vector for every integer budget in $[0,\ C_{\text{budget}}]$.

---

## Key Results

### Efficient Frontier (Marginal Allocation)

The EBO–cost curve is strictly decreasing and convex, reflecting diminishing returns as the spare parts budget grows. The first iteration selects the LRU with the best marginal efficiency ratio, not necessarily the one with the highest failure rate — illustrating the difference between absolute and relative optimisation criteria.

### MA vs DP Comparison

At every budget level corresponding to a Marginal Allocation efficient point, Dynamic Programming produces an **identical allocation and EBO** — confirming that MA correctly identifies the convex hull of the solution set. At intermediate (non-efficient) budgets, DP solutions lie slightly above the MA curve, as expected: the MA curve is a piecewise-linear interpolation (convex hull), while DP solves the exact integer problem.

### Baseline Evaluation

Evaluating the naive "one spare per LRU type" baseline ($\mathbf{s}_{\text{base}} = [1,1,\ldots,1]$) against the DP-optimal frontier shows that this allocation is **not on the efficient frontier**:

- **Cost minimisation**: An allocation exists that achieves the same EBO level at a cost saving of ~5%.
- **Performance maximisation**: Within the same budget, a re-optimised allocation reduces EBO by ~0.23 units.

This result demonstrates that uniform allocation is suboptimal, and systematic optimisation yields meaningful improvements.

---

## Repository Structure

```
spare-parts-optimization/
├── main.m                  # Main analysis script — runs both algorithms end-to-end
├── generate_lru_data.m     # LRU parameter generator (failure rates, repair times, costs)
├── docs/
│   └── methodology.md      # Detailed mathematical derivations and algorithm walk-through
└── README.md
```

---

## How to Run

**Requirements:** MATLAB R2020a or later. No additional toolboxes required.

```matlab
% Clone or download the repository, then in MATLAB:
cd spare-parts-optimization
main
```

The script will:
1. Generate LRU parameters (reproducible via fixed seed)
2. Precompute EBO/R/p tables for all stock levels
3. Run Marginal Allocation and print the first 5 efficient points
4. Run Dynamic Programming and print optimal allocations at key budget levels
5. Compare both methods against the uniform baseline
6. Display two figures (efficient frontier; MA vs DP overlay)

To experiment with a different system configuration, edit the constants at the top of `main.m`:

```matlab
RANDOM_SEED = 42;    % Change to generate a different LRU instance
C_BUDGET    = 500;   % Total budget available
MAX_SPARES  = 50;    % Maximum spares considered per LRU type
```

---

## Topics

`Operations Research` · `Inventory Optimisation` · `Dynamic Programming` · `Marginal Allocation` · `Poisson Processes` · `Expected Backorders` · `Multi-Objective Optimisation` · `MATLAB`
