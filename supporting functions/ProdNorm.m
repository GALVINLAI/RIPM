function val = ProdNorm(M, x, s)

% Norm for product tangent space of
% M x Euc_eq x Euc_ineq x Euc_ineq 
% appeared in RIPM.

% function val = ProdNorm(M, x, s)

% only for 
% s.dx, s.dy,
% and
% s.dx, s.dy, s.dz, s.ds.

val =  sqrt(ProdInner(M, x, s, s));

end