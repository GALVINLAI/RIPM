# Riemannian Interior Point Methods (RIPM)

Riemannian Interior Point Methods (RIPM) for Constrained Optimization on Manifolds

This is a *primal-dual interior point methods* solver for nonlinear
optimization problems on Riemannian manifolds, which aims to minimize the
cost function in the given problem structure with (in)equality constraints：

$$
\begin{array}{cl}
\min _{x \in \mathbb{M}} & f(x) \\
\text { s.t. } & h(x)=0, \text { and } g(x) \leq 0.
\end{array}
$$

For a description of the algorithm and theorems offering convergence guarantees, see the references:\
Z. Lai and A. Yoshise. *Riemannian Interior Point Methods for Constrained Optimization on Manifolds.* [[arXiv]](https://arxiv.org/abs/2203.09762)

**We prepared a additional document [Implement Note of Global RIPM](Implement_Note_of_Global_RIPM.pdf) with detailed instructions for the implementation.**

**[New!] We upload a [new toturial](NewNoteOfRIPM.pdf)**


The code is based on matalb solver 'Manopt'.\
Before running the codes, you must install the solver ['Manopt'](https://www.manopt.org/tutorial.html)!
