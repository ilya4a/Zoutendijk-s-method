# ZoitendijkMethod

A small C++ library for solving constrained optimization problems with the Zoutendijk feasible directions method.

## Features

- Automatic differentiation via `autodiff`
- Auxiliary simplex subproblems solved with `CLP`
- Linear algebra based on `Eigen`
- Public API wrapped through a function interface (`FuncWrap`)
- Example usage provided in `examples/`


## Dependencies

- CMake 3.16+
- Eigen3
- autodiff
- CLP

## Build

```bash
cmake -S . -B build
cmake --build build
```

## 2.1.2 Problem Statement

We consider the constrained optimization problem

$$
\min \varphi_0(x), \qquad x \in S,
$$

where

$$
S = \{x \in \mathbb{R}^n \mid \varphi_i(x) \le 0,\ i=1,\dots,m,\ Ax=b\},
$$

$$
A \in \mathbb{R}^{l \times n}, l < n, and b \in \mathbb{R}^l
$$

The method is applied under the standard assumptions of compact feasibility, Slater’s condition, convexity and smoothness of the objective and constraints, and bounded/Lipschitz gradients.



## Algorithm Description

At each iteration, the method works with the set of nearly active constraints

$$
J_\delta(x) = \{\, i \in \{1,\dots,m\} \mid -\delta \le \varphi_i(x) \le 0 \,\}.
$$

Only the constraints from this set are used in the linearized auxiliary problem.

To find a feasible descent direction \(s_k\), the following auxiliary linear programming problem is solved:

$$
\min \eta
$$

subject to

$$
\nabla \varphi_0(x_k)^T s \le \xi_0 \eta,
$$

$$
\nabla \varphi_i(x_k)^T s \le \xi_i \eta, \qquad i \in J_{\delta_k}(x_k),
$$

$$
A s = 0,
$$

$$
-1 \le s^{(j)} \le 1, \qquad j = 1,\dots,n.
$$

Here $\xi_0,\xi_1,\dots,\xi_m > 0$ are weighting coefficients (they are needed for special problems and can be taken as 1), and the bound on \(s\) makes the subproblem a standard linear program. In this project, it is solved using an external simplex solver.

The main update step accepts the direction only if it provides sufficient decrease and preserves feasibility; otherwise, the parameter $\delta_k$ is reduced.

A key feature of the method is that it uses only the nearly active constraints $J_\delta(x)$, which makes it effective even for problems with a large number of nonlinear constraints.
### Initial approximation

The starting point $x_0$ is built from a separate feasibility-type auxiliary problem. If the initial point is not feasible, the algorithm first drives it into the feasible region and only then proceeds with the main optimization loop.

The initial feasible point is obtained from the auxiliary problem

$$
\min \eta
$$

subject to

$$
\varphi_i(x)\le \eta,\qquad i=1,\dots,m,
$$

$$
Ax=b.
$$

### Step update

If the computed value $\eta_k$ is sufficiently negative, the method moves in the direction \(s_k\) with a step length chosen by backtracking. The step is accepted only if the objective decreases and all constraints remain satisfied. Otherwise, the threshold parameter $\delta_k$ is reduced.
