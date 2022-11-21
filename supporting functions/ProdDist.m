function val = ProdDist(M, x, s1, s2)

% Distance for product tangent space of
% M x Euc_eq x Euc_ineq x Euc_ineq 
% appeared in RIPM.

% function val = ProdDist(M, x, s1, s2)

% only for 
% s.dx, s.dy,
% and
% s.dx, s.dy, s.dz, s.ds.

val = ProdNorm(M, x, ProdLincomb(M, x, 1, s1, -1, s2));

end