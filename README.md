# RIPM

**Riemannian Interior Point Methods (RIPM)** for Constrained Optimization on Manifolds

This is a *primal-dual interior point methods* solver for nonlinear
optimization problems on Riemannian manifolds, which aims to minimize the
cost function in the given problem structure with (in)equality
constraints.

For a description of the algorithm and theorems offering convergence
guarantees, see the references below:\
Z. Lai and A. Yoshise.
*Riemannian Interior Point Methods for Constrained Optimization on Manifolds.* [[arXiv]](https://arxiv.org/abs/2203.09762)

An earlier version of this article has been circulated under the title
"Superlinear and Quadratic Convergence of Riemannian Interior Point Methods".

In addition, we have a separate document with detailed instructions for
the implementation. The documentation is available on github.

The code is based on matalb solver 'Manopt'.\
Before running the codes, you must install the solver ['Manopt'](https://www.manopt.org/tutorial.html)!
