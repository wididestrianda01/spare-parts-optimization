# Methodology: Spare Parts Optimisation

## 1. Problem Setting

A system consists of $N = 9$ **Line Replaceable Units (LRUs)** — independent components that can fail and be replaced. When an LRU fails, it is sent for repair. During the repair period (duration $T_j$), if no spare unit is on the shelf, the system experiences a **backorder**: it cannot be restored to full operational status.

Each LRU $j$ has the following parameters:

| Parameter | Symbol | Description |
|---|---|---|
| Failure rate | $\lambda_j$ | Failures per unit time |
| Repair time | $T_j$ | Mean time until repaired unit is returned |
| Unit cost | $c_j$ | Cost of purchasing one spare unit |
| Mean demand | $\mu_j = \lambda_j T_j$ | Expected demand during one repair cycle |

The number of demands during a repair cycle follows a **Poisson distribution** with mean $\mu_j$.

**Objective:** Choose stock levels $s_j \geq 0$ (integers) to minimise total Expected Backorders subject to a total budget $C_{\text{budget}}$:

$$\min_{\mathbf{s} \in \mathbb{Z}_{\geq 0}^N} \; \sum_{j=1}^{N} \text{EBO}_j(s_j) \quad \text{s.t.} \quad \sum_{j=1}^{N} c_j s_j \leq C_{\text{budget}}$$

---

## 2. Recursive EBO Formulae

For a Poisson($\mu$) demand and $s$ spare units on the shelf, the key quantities are derived from the **METRIC / RESOPT** framework.

### Probability Mass Function

$$p(0) = e^{-\mu}, \qquad p(s+1) = \frac{\mu}{s} \cdot p(s), \quad s \geq 1$$

### Shortage Probability

$R(s)$ is the probability that demand exceeds stock (a shortage occurs):

$$R(0) = 1 - p(0), \qquad R(s+1) = R(s) - p(s+1)$$

### Expected Backorders

$\text{EBO}(s)$ is the expected number of unfilled demands:

$$\text{EBO}(0) = \mu, \qquad \text{EBO}(s+1) = \text{EBO}(s) - R(s)$$

These recursions are $O(N \cdot S_{\max})$ to compute and are precomputed once for stock levels $s = 0, 1, \ldots, S_{\max}$.

### Integer Convexity

Two properties are critical for the Marginal Allocation algorithm to be valid:

**Property 1 (strictly decreasing):**
$$\Delta \text{EBO}(s) = \text{EBO}(s+1) - \text{EBO}(s) = -R(s) < 0$$

Since $R(s) > 0$ for all finite $s$, the EBO is strictly decreasing in $s$.

**Property 2 (integer-convex):**
$$\Delta^2 \text{EBO}(s) = \text{EBO}(s+2) - 2\,\text{EBO}(s+1) + \text{EBO}(s) = p(s+1) \geq 0$$

Since $p(s+1) \geq 0$, the EBO function satisfies the discrete convexity condition.

The cost function $g_j(s_j) = c_j \cdot s_j$ is linear, hence trivially strictly increasing and integer-convex.

Both conditions are verified numerically in `main.m` before the algorithms run.

---

## 3. Marginal Allocation Algorithm

### Motivation

The multi-item stocking problem is a **multi-objective optimisation** problem: minimise EBO and minimise cost. The set of Pareto-optimal solutions forms the *efficient frontier*. Because the objective and constraint functions are integer-convex, the Marginal Allocation algorithm traces this frontier exactly via a greedy strategy.

### Algorithm

**Initialise:** $\mathbf{s}^{(0)} = \mathbf{0}$, $C^{(0)} = 0$, $\text{EBO}^{(0)} = \sum_j \mu_j$.

**Iterate (k = 1, 2, …):**

1. Compute the marginal efficiency quotient for each LRU at its current stock level:
$$Q_j(s_j) = \frac{R_j(s_j)}{c_j}$$

2. Select the LRU with the highest quotient:
$$\ell^* = \arg\max_j \; Q_j(s_j)$$

3. If $C^{(k-1)} + c_{\ell^*} > C_{\text{budget}}$, **stop**.

4. Update: $s_{\ell} \leftarrow s_{\ell} + 1$, $\; C \leftarrow C + c_{\ell}$, $\; \text{EBO} \leftarrow \text{EBO} - R_{\ell}(s_{\ell})$.

5. Record the new state as an **efficient point**.

### Correctness

Because the EBO function is integer-convex, the efficiency quotient $Q_j(s_j)$ is non-increasing in $s_j$. This ensures the greedy selection always produces a point on the southwestern boundary (convex hull) of the feasible objective space.

### Complexity

$O(K \cdot N)$ where $K$ is the total number of spares purchased (at most $C_{\text{budget}} / \min_j c_j$).

---

## 4. Dynamic Programming Algorithm

### Formulation

The problem is cast as a **Distribution-of-Effort** problem. Each of the $N$ LRUs corresponds to one *stage*. The state variable $s_n$ is the remaining budget available to stages $n, n+1, \ldots, N$.

**Decision:** $x_{n} \in F_{n}(s_{n}) = \lbrace 0,\ c_{n},\ 2c_{n},\ \ldots \rbrace \cap [0,\ s_{n}]$ — the budget allocated to LRU $n$.

**State transition:** $s_{n+1} = s_{n} - x_{n}$ (linear, so the Principle of Optimality holds).

**Stage contribution:** $P_{n}(x_{n}) = \text{EBO}_{n}\!\left(\dfrac{x_{n}}{c_{n}}\right)$.

**Bellman recursion (backward, from stage $N$ to stage $1$):**

$$f_{n}^{*}(s_{n}) = \min_{x_{n} \in F_{n}(s_{n})} \left\{ P_{n}(x_{n}) + f_{n+1}^{*}(s_{n} - x_{n}) \right\}$$

**Terminal condition:** $f_{N+1}^{*}(s) = 0 \;\forall s$ (no future stages).

### Implementation

The backward recursion fills a table $f^{*}$ of size $(N+1) \times (C_{\text{budget}}+1)$ and records the optimal policy $x_{n}^{*}(s_{n})$. A forward reconstruction pass then traces the optimal allocation vector for any desired budget level.

### Complexity

$O(N \cdot C_{\text{budget}}^2 / \min_j c_j)$ in the worst case; in practice, the inner loop over feasible decisions is $O(s / c_n)$ per state, making the total complexity manageable for $C_{\text{budget}} = 500$.

### Monotonicity

Unlike Marginal Allocation, Dynamic Programming does **not** guarantee that allocations are non-decreasing in the budget. As the budget increases, it may become worthwhile to reallocate funds away from a cheaper-but-less-effective LRU toward a more expensive-but-highly-effective one. This was confirmed numerically: 126 monotonicity violations were identified across all integer budgets in $[0, 500]$, consistent with the theory.

---

## 5. Comparison: MA vs DP

| Property | Marginal Allocation | Dynamic Programming |
|---|---|---|
| Solution quality | Efficient points only | Globally optimal for every budget |
| Agreement | Identical at efficient budget levels | Identical at MA efficient points |
| Intermediate budgets | Linear interpolation (convex hull) | Exact integer solution |
| Allocation monotonicity | Guaranteed (never removes a spare) | Not guaranteed |
| Computational cost | Low — $O(KN)$ | Higher — $O(N \cdot C^2 / c_{\min})$ |

The two methods are complementary: MA quickly produces the efficient frontier (useful for budget planning and trade-off analysis), while DP finds exact solutions for any specific budget target.

---

## 6. Baseline Suboptimality Analysis

A uniform baseline allocation $\mathbf{s}_{\text{base}} = [1, 1, \ldots, 1]$ (one spare per LRU type) is a common naive policy. Evaluating it against the DP-optimal frontier reveals:

- **It does not lie on the efficient frontier**, meaning it is simultaneously more expensive than necessary for its EBO level and worse in EBO than achievable within its budget.
- Re-optimising the allocation within the same budget reduces EBO by approximately 7%.
- Alternatively, the same EBO level can be achieved at approximately 5% lower cost.

This quantifies the value of systematic optimisation over uniform stocking policies.
