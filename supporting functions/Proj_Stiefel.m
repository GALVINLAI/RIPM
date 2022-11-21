function Z= Proj_Stiefel(X)
% Projection onto Stiefel manifold whose size is euqal to X.

[U,~,V]=svd(X,'econ');
Z=U*V';

% M \in \mathbb{R}^{n \times k}
% [U,S,V] = svd(M, >econ=); Gives a compact form of SVD for both n < k and n >= k.

% Finally, $\mathcal{P}_{\mathcal{S}_{n, p}}(X)=U V^{\top}$ denotes the
% projection to Stiefel manifold, where $X=U \Sigma V^{\top}$ is
% the economic SVD of $X$ with $U \in \mathcal{S}_{n, p}, V \in
% \mathcal{S}_{p, p}$ and $\Sigma$ is $p \times p$ diagonal matrix with the
% singular values of $X$ on its diagonal.

% SEE
% http://www.optimization-online.org/DB_FILE/2021/10/8628.pdf

end

