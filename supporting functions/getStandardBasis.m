function basis = getStandardBasis(M)
% only for problem.Euc_eq is one of the following:
%
% euclideanfactory; skewsymmetricfactory; symmetricfactory; euclideansparsefactory
%

name = M.name();
dim = M.dim();

%% euclideanfactory(m, n)
if contains(name,'Euclidean space R^')
    basis = cell(dim, 1);
    ZERO = M.zerovec();
    for k = 1 : dim
        basis{k} = ZERO;
        basis{k}(k) = 1;
    end
end

%% skewsymmetricfactory(n)
if contains(name,'Skew-symmetric matrices of size')
    basis = getEuclideanVecBasis(dim);
    for k = 1 : dim
        basis{k} = Vec2Skew(basis{k});
    end
end

%% symmetricfactory(n)
if contains(name,'Symmetric matrices of size')
    basis = getEuclideanVecBasis(dim);
    for k = 1 : dim
        basis{k} = Vec2Sym(basis{k});
    end
end
% e.g. for M=symmetricfactory(3): standard basis are
% 1	0 0 % 0	0 0 % 0	0 0 
% 0	0 0 % 0	1 0 % 0	0 0
% 0	0 0 % 0	0 0 % 0	0 1
% and
% 0	1 0 % 0	0 1 % 0	0 0
% 1	0 0 % 0	0 0 % 0	0 1
% 0	0 0 % 1	0 0 % 0	1 0
% But the latter three are not unit normally under unsual inner product.
% All of them have norm sqrt(2)!!
% This point is the same to skewsymmetricfactory(n).

%% euclideansparsefactory(A)
if contains(name,'Euclidean space R^') && contains(name,'with fixed sparsity pattern containg')
    basis = cell(dim, 1);
    Pattern = M.rand(); % 
    [m, n] = size(Pattern);
    nonzero_indix = find(Pattern); % linear index
    [row, col] = ind2sub([m, n], nonzero_indix); % linear index 2 subscript index
    S = sparse(m, n);
    for k = 1 : dim
        basis{k} = S; 
        basis{k}(row(k),col(k)) = 1;
    end
end

end


%% local funs
function basis = getEuclideanVecBasis(n)
    basis = cell(n, 1);
    ZERO = zeros(n,1);
    for k = 1 : n
        basis{k} = ZERO;
        basis{k}(k) = 1;
    end
end




