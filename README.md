# Riemannian Interior Point Methods (RIPM)

## 1. overview

Riemannian Interior Point Methods (RIPM) for Constrained Optimization on Manifolds

This is a *primal-dual interior point methods* solver for nonlinear optimization problems on Riemannian manifolds, which aims to minimize the cost function in the given problem structure with (in)equality constraints

$$
\begin{array}{cl}
\min _{x \in \mathbb{M}} & f(x) \\
\text { s.t. } & h(x)=0, \text { and } g(x) \leq 0.
\end{array}
$$

For a description of the algorithm and theorems offering convergence guarantees, see the references:\
Z. Lai and A. Yoshise. *Riemannian Interior Point Methods for Constrained Optimization on Manifolds.* [[arXiv]](https://arxiv.org/abs/2203.09762)

**We prepared a additional document [Implement Note of Global RIPM](Implement_Note_of_Global_RIPM.pdf) with detailed instructions for the implementation.**

**[New!] We upload a new [tutorial](tutorial.pdf) for quick start.**



## 2. Numerical Experiments and Codes

### Part I

The codes below correspond to the **Numerical Experiments Section** in the paper. 

By running the `Boss_*.m` files, the numerical experiments are executed with the same settings as in the paper. The results will be output in `csv` format in the `.\numerical results` folder. Note: It may take 1-2 days to run the full experiment through the `.\bosses\cmd.m` file.

| MATLAB code (.\bosses\\)                   | Corresponding questions in our paper |
| ------------------------------------------ | ------------------------------------ |
| Boss_1_fixedrank_NLRM.m                    | Problem I (NLRM)                     |
| Boss_2_Stiefel_NonnegativeProjection.m     | Problem II (Model_St)                |
| Boss_3_Ob_equality_NonnegativeProjection.m | Problem II (Model_Ob)                |

The corresponding `client_*.m` files combines specific problems with various solvers. They are used in the `Boss_*.m` files.


| MATLAB code (.\clients\\)                  | Corresponding questions in our paper |
| ------------------------------------------ | ------------------------------------ |
| client_fixedrank_NLRM.m                    | Problem I (NLRM)                     |
| client_Stiefel_NonnegativeProjection.m     | Problem II (Model_St)                |
| client_Ob_equality_NonnegativeProjection.m | Problem II (Model_Ob)                |

### Part II

The codes below do not appear in the paper and are used to test various example problems for our RIPM.

| MATLAB code (.\examples\\)                    | M         | f                 | g           | h                      |
| --------------------------------------------- | --------- | ----------------- | ----------- | ---------------------- |
| Euc_linear_nonnega_sphereEq                   | Euclidean | linear            | nonnegative | sphere x'*x=1          |
| Euc_linear_or_projection_nonnega_orthogonalEq | Euclidean | linear/projection | nonnegative | orthogonality X'*X-I=0 |
| Euc_projection_nonnega                        | Euclidean | projection        | nonnegative | -                      |
| Euc_projection_nonnega_symmetricEq            | Euclidean | projection        | nonnegative | symmetry X-X'=0        |
| Fixedrank_MatrixCompletion_nonnega            | fixedrank | matrix completion | nonnegative | -                      |
| Fixedrank_MatrixCompletion_nonnega_reliableEq | fixedrank | matrix completion | nonnegative | reliable sampled data  |
| Fixedrank_projection_nonnega                  | fixedrank | projection        | nonnegative | -                      |
| Ob_ONMF_StiefelEq                             | oblique   | -trace(X'*AAt*X)  | nonnegative | norm(X*V,'fro')^2-1    |
| Ob_linear_or_projection_nonnega_StiefelEq     | oblique   | linear/projection | nonnegative | norm(X*V,'fro')^2-1    |
| Sp_linear_nonnega                             | sphere    | linear            | nonnegative | -                      |
| Sp_linear_nonnega_linearEq                    | sphere    | linear            | nonnegative | linear                 |
| Sp_quadratic_nonnega                          | sphere    | quadratic         | nonnegative | -                      |
| Stiefel_linear_or_projection_nonnega          | stiefel   | linear/projection | nonnegative | -                      |
| Sym_projection_nonnega                        | symmetric | projection        | nonnegative | -                      |
| run_all_examples                              |           |                   |             |                        |

