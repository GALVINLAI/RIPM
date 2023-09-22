function prod = ProdLincomb(M, x, a1, s1, a2, s2)
% 线性结构. Euc_eq是具体那种子空间没有影响.

% Linear combination function for product tangent space of
% M x Euc_eq x Euc_ineq x Euc_ineq 
% appeared in RIPM.

% function prod = ProdLincomb(M, x, a1, s1)
% function prod = ProdLincomb(M, x, a1, s1, a2, s2)

% only for 
% s.dx, s.dy,
% and
% s.dx, s.dy, s.dz, s.ds.

if nargin == 4 % for ProdLincomb(M, x, a1, s1)
    prod.dx = M.lincomb(x, a1, s1.dx);% This is valid since fun matrixlincomb() used in Manopt.
    % 下段代码块足以应对结构体s的两种类型.
    prod.dy = a1*s1.dy; % note that 1 * [] = [] in Matlab!
    if isfield(s1, 'dz')
        %MODIFIED: For RIPM on ProductManifold
        if isstruct(s1.dz)  
            prod.dz = multiplyStructWithScalar(s1.dz, a1);
            prod.ds = multiplyStructWithScalar(s1.ds, a1);
        else
            prod.dz = a1*s1.dz;
            prod.ds = a1*s1.ds;
        end
        % prod.dz = a1*s1.dz;
        % prod.ds = a1*s1.ds;
    end
elseif nargin == 6 % for ProdLincomb(M, x, a1, s1, a2, s2)
    prod.dx = M.lincomb(x, a1, s1.dx, a2, s2.dx);
    % 下段代码块足以应对结构体s的两种类型.
    prod.dy = a1*s1.dy + a2*s2.dy; % note that  [] + [] = [] in Matlab!
    if isfield(s1, 'dz') && isfield(s2, 'ds')
        %MODIFIED: For RIPM on ProductManifold
        if isstruct(s1.dz)  
            prod.dz = addStructs(multiplyStructWithScalar(s1.dz, a1), multiplyStructWithScalar(s2.dz, a2));
            prod.ds = addStructs(multiplyStructWithScalar(s1.ds, a1), multiplyStructWithScalar(s2.ds, a2));
        else
            prod.dz = a1*s1.dz + a2*s2.dz;
            prod.ds = a1*s1.ds + a2*s2.ds;
        end
        % prod.dz = a1*s1.dz + a2*s2.dz;
        % prod.ds = a1*s1.ds + a2*s2.ds;
    end
else
    error('ProdLincomb takes either 4 or 6 inputs.');
end

end