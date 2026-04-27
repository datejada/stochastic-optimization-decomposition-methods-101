# Short Introduction to Stochastic Optimization

## Deterministic vs Stochastic Fixed-Charge Transportation (FCTP/SFCTP)

---

> **Audience**  
> This note is for readers who already know LP/MIP modeling and want a first, practical introduction to stochastic optimization through the fixed-charge transportation example in this project.

---

## 1. Why Stochastic Optimization?

In deterministic optimization, uncertain inputs (like demand) are replaced by a single value (one scenario, or an average). This is often easy to model, but it can produce decisions that are fragile when reality differs from that single value.

Stochastic optimization explicitly models uncertainty through a set of scenarios with probabilities and optimizes expected performance, often with recourse decisions after uncertainty is revealed.

For two-stage models, the idea is:

1. **First stage (here-and-now):** decisions made before demand is known.
2. **Second stage (wait-and-see / recourse):** operational adjustments after a scenario is realized.

In the SFCTP code, the first-stage decision is arc investment $Y_{ij}$; second-stage decisions are flows $X_{ij}^s$ and unmet-demand slack $Z_j^s$.

---

## 2. Deterministic FCTP in [FCTP.jl](FCTP.jl)

The deterministic fixed-charge transportation model solves:

$$\min \sum_{i,j} F_{ij}Y_{ij} + \sum_{i,j} C_{ij}X_{ij}$$

subject to supply, demand, and linking constraints:

$$\sum_j X_{ij} \le A_i, \quad \sum_i X_{ij} \ge B_j, \quad X_{ij} \le \min(A_i,B_j)Y_{ij}, \quad Y_{ij}\in\{0,1\}$$

### What this means

- You choose which arcs to invest in via binary $Y_{ij}$ and pay fixed costs $F_{ij}$.
- Then you ship flow $X_{ij}$ and pay variable costs $C_{ij}$.
- Demand vector $B$ is assumed known and fixed.

### In the file

[FCTP.jl](FCTP.jl) solves multiple deterministic instances:

1. Scenario 1 demand: $[21,51,31]$
2. Scenario 2 demand: $[32,22,52]$
3. Scenario 3 demand: $[53,33,23]$
4. Weighted average-demand case

This is useful for comparison, but each solve assumes one demand realization is "the truth".

---

## 3. Stochastic SFCTP in [SFCTP.jl](SFCTP.jl)

The stochastic model introduces scenarios $s=1,\dots,S$ with probabilities $p_s$ and scenario-dependent demand $B_{sj}$.

In the example data:

- $S=3$
- $p = [0.5,\;0.3,\;0.2]$
- scenario demands are rows of matrix $B$

The SFCTP objective is:

$$\min \sum_{i,j} F_{ij}Y_{ij} + \sum_s p_s\left(\sum_{i,j} C_{ij}X_{ij}^s + \kappa\sum_j Z_j^s\right)$$

with $\kappa=100$ as non-supplied demand penalty (`NSD_COST`).

### Key difference from deterministic model

- **Deterministic:** one $B_j$, one operational plan.
- **Stochastic:** one shared investment plan $Y$ across all scenarios, plus scenario-specific recourse $(X^s,Z^s)$.

This non-anticipativity (same $Y$ for all scenarios) is exactly what makes the model "stochastic programming" instead of solving each scenario independently.

---

## 4. Side-by-Side Comparison

| Feature | Deterministic FCTP | Stochastic SFCTP |
|---|---|---|
| Demand input | Single vector $B$ | Scenario matrix $B_{sj}$ + probabilities $p_s$ |
| Investment decision $Y$ | Optimized for one demand case | Single robust plan across all scenarios |
| Flow decision | One flow matrix $X$ | One flow matrix per scenario $X^s$ |
| Shortage handling | Typically infeasible if demand exceeds supply (in [FCTP.jl](FCTP.jl), infeasible by assertion) | Explicit shortage slack $Z^s$ with penalty |
| Objective meaning | Cost for one assumed future | Expected cost over uncertain futures |

---

## 5. Stochastic Performance Measures in [SFCTP.jl](SFCTP.jl)

The file computes the standard quality metrics used in stochastic programming.

## 5.1 Recourse Problem Value (RP)

This is the optimal objective of the stochastic model (the true two-stage model):

$$RP = z_{SP}$$

In code: `recourse_problem_solution = objective_value(stochastic)`.

Interpretation: best expected cost when uncertainty is modeled directly.

## 5.2 Wait-and-See (WS)

Solve each scenario independently as if you knew it in advance, then take expected value:

$$WS = \sum_s p_s z_s^{\text{det-opt}}$$

In code: `wait_and_see_solution = sum(P[s] * objective_value(determ_sc[s]) for s = 1:S)`.

Interpretation: lower bound for minimization (perfect foresight benchmark). Notice that this lower bound is unachievable in practice because you can't know the future demand with certainty.

## 5.3 EVPI (Expected Value of Perfect Information)

$$EVPI = RP - WS$$

In code: `expected_value_perfect_information = recourse_problem_solution - wait_and_see_solution`.

Interpretation: maximum amount worth paying for perfect demand information. For minimization, typically nonnegative.

## 5.4 Expected Value (EV) solution and EEV

1. Build an expected-demand instance with
$$\bar{B}_j = \sum_s p_s B_{sj}$$
2. Solve it to get investment decision $Y^{EV}$.
3. Fix $Y=Y^{EV}$ in each real scenario model and re-optimize recourse.
4. Compute expected resulting cost:
$$EEV = \sum_s p_s z_s(Y^{EV})$$

In the example code this is exactly steps 4.1 to 4.4.

Interpretation: performance of the "optimize the average" policy under true uncertainty.

## 5.5 VSS (Value of the Stochastic Solution)

$$VSS = EEV - RP$$

In code: `value_stochastic_solution = expectation_expected_value_solution - recourse_problem_solution`.

Interpretation: expected benefit of solving the stochastic model instead of the average-value deterministic model. For minimization, typically nonnegative.

---

## 7. Out-of-Sample Check in [SFCTP.jl](SFCTP.jl)

The script also performs an out-of-sample test with demand $[20,20,70]$:

1. Compute the out-of-sample optimum (benchmark for that specific unseen demand).
2. Fix $Y$ from stochastic solution and evaluate cost + unmet demand.
3. Fix $Y$ from deterministic scenario-1 solution and evaluate again.

This is excellent practice because it tests **policy robustness**, not only in-sample objective value.

---

## 8. Practical Guidelines

1. Use deterministic FCTP when uncertainty is genuinely negligible or when you only need a quick baseline.
2. Use SFCTP when first-stage decisions are expensive/hard to change and demand uncertainty is material.
3. Report `RP`, `WS`, `EVPI`, `EEV`, and `VSS` together; each one answers a different managerial question.
4. Add out-of-sample stress tests (as in the example code) to avoid overconfidence in in-sample performance.

---

## 9. Minimal Formula Summary

$$RP = z_{SP}, \qquad WS = \sum_s p_s z_s^{\text{det-opt}}, \qquad EVPI = RP - WS$$

$$\bar{B} = \sum_s p_s B_s, \qquad EEV = \sum_s p_s z_s(Y^{EV}), \qquad VSS = EEV - RP$$

For minimization problems under standard assumptions:

$$WS \le RP \le EEV \quad \Rightarrow \quad EVPI \ge 0, \; VSS \ge 0$$

---

## 10. References

- Birge, J. R., & Louveaux, F. (2011). *Introduction to Stochastic Programming* (2nd ed.). Springer.
- Shapiro, A., Dentcheva, D., & Ruszczynski, A. (2009). *Lectures on Stochastic Programming: Modeling and Theory* (2nd ed.). SIAM.
- Kall, P., & Wallace, S. W. (1994). *Stochastic Programming*. Wiley.
- Rockafellar, R. T., & Wets, R. J.-B. (1991). Scenarios and policy aggregation in optimization under uncertainty. *Mathematics of Operations Research*, 16(1), 119-147.

---

*Companion files in this project:*  

- Deterministic model: [FCTP.jl](FCTP.jl)
- Stochastic model + metrics: [SFCTP.jl](SFCTP.jl)  
- Advanced decomposition methods: [SFCTP-CG.jl](SFCTP-CG.jl), [SFCTP-CG-B&P.jl](SFCTP-CG-B&P.jl)
