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

\[
\min \varphi_0(x), \qquad x \in S,
\]

where

\[
S = \{x \in \mathbb{R}^n \mid \varphi_i(x) \le 0,\ i=1,\dots,m,\ Ax=b\},
\]

\(A \in \mathbb{R}^{l \times n}\), \(l < n\), and \(b \in \mathbb{R}^l\).

The method is applied under the standard assumptions of compact feasibility, Slater’s condition, convexity and smoothness of the objective and constraints, and bounded/Lipschitz gradients.

## 2.1.4 Algorithm Description

The Zoutendijk feasible directions method solves the original nonlinear constrained problem iteratively by constructing a descent direction at each step.

At iteration \(k\), an auxiliary linear programming problem is solved to find a feasible direction \(s_k\) and a value \(\eta_k\):

- the gradient of the objective and the gradients of the nearly active constraints are included;
- equality constraints are enforced by \(As = 0\);
- the direction is bounded component-wise, which makes the auxiliary problem linear.

This auxiliary problem is solved with the simplex method from an external library. In this project, the simplex solver is used only for the inner linear subproblem, while the nonlinear gradients are obtained automatically.

### Initial approximation

The starting point \(x_0\) is built from a separate feasibility-type auxiliary problem. If the initial point is not feasible, the algorithm first drives it into the feasible region and only then proceeds with the main optimization loop.

### Step update

If the computed value \(\eta_k\) is sufficiently negative, the method moves in the direction \(s_k\) with a step length chosen by backtracking. The step is accepted only if the objective decreases and all constraints remain satisfied. Otherwise, the threshold parameter \(\delta_k\) is reduced.

### Main feature

A key advantage of the method is that it works with the set of nearly active constraints \(J_\delta(x)\). Because of this, it remains effective even when the problem contains many nonlinear constraints: at each iteration only the constraints close to being active are used in the auxiliary subproblem, which keeps the linearized problem compact and tractable.