function val = ProdInner(M, x, s1, s2)
% 特别注意: 默认Euc_eq x Euc_ineq x Euc_ineq 装配的是普通内积,
% 但是如果Euc_eq = Skew space, 普通内积或可能影响对称性.

% Inner product for product tangent space of
% M x Euc_eq x Euc_ineq x Euc_ineq 
% appeared in RIPM.

% function val = ProdInner(M, x, s1, s2)

% only for 
% s.dx, s.dy,
% and
% s.dx, s.dy, s.dz, s.ds.

if isempty(s1.dy) 
    % note that 1 + [] = [] in Matlab!
    val = M.inner(x, s1.dx, s2.dx);
else
    val = M.inner(x, s1.dx, s2.dx) + s1.dy(:)'*s2.dy(:);
end

if isfield(s1, 'dz') && isfield(s2, 'ds')
    %MODIFIED: For RIPM on ProductManifold
    if isstruct(s1.dz)
        val = val + innerProductStructs(s1.dz,s2.dz) + innerProductStructs(s1.ds,s2.ds);
    else
        val = val + s1.dz(:)'*s2.dz(:) + s1.ds(:)'*s2.ds(:);
    end
    %val = val + s1.dz(:)'*s2.dz(:) + s1.ds(:)'*s2.ds(:);
end

end