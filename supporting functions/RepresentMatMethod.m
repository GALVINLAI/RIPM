function [NTdir, RepresentMat, RepresentMatOrder] = RepresentMatMethod(Aw, Hxaj, b, M, x, y, condet, problem)

warning('off','MATLAB:nearlySingularMatrix');

%% Case 1. has_eq_cost = false.
% In this case, T = Aw is an operator to and from TxM.
if ~condet.has_eq_cost
    try
        Bx = getStandardBasis(M); % when M is trivial (linear) manifold.
    catch
        Bx = tangentorthobasis(M, x); % an orthonormal basis of TxM
    end

    % a matrix with respect to Bx; Aw_mat is symmetric.
    % the next code is equal to run:
    % Aw_mat = operator2matrix(M, x, x, Aw, Bx, Bx);
    % operator2matrix is given by Manopot for general operator, without using symmetry of Aw.
    Aw_mat = SelfAdj_operator2matrix(M, x, Aw, Bx);

    c_vec = tangent2vec(M, x, Bx, b.dx); % Expands tangent vector c into an orthonormal basis Bx

    % direct method
    sol_vec = linsolve(Aw_mat, c_vec, struct('SYM', true));
    NTdir.dx = lincomb(M, x, Bx, sol_vec);% restore solution; note that it is not M.lincomb()
    NTdir.dy = [];
    RepresentMat = Aw_mat;
    RepresentMatOrder = M.dim();
end


%% Case 2. has_eq_cost = ture.
% In this case, T is an operator to and from TxM x Euc_eq.
if condet.has_eq_cost
    % get product manifold M_Euc_eq := M x Euc_eq
    M_Euc_eq = productmanifold(struct('dx', M, 'dy', problem.Euc_eq));
    xy.dx = x; xy.dy = y; % xy = (x,y) in M x Euc_eq

    % get basis Bxy of M_Euc_eq
    try
        Bx = getStandardBasis(M); % when M is trivial (linear) manifold.
    catch
        Bx = tangentorthobasis(M, x); % an orthonormal basis of TxM
    end

    By = getStandardBasis(problem.Euc_eq);

    dimM = problem.dimM;
    dimEuc_eq = problem.dimEuc_eq;
    dimProd = dimM + dimEuc_eq;
    Bxy = cell(1, dimProd);
    for i = 1: dimM
        Bxy{i}.dx = Bx{i};
        Bxy{i}.dy = problem.ZeroTyEuc;
    end
    for i = dimM+1: dimProd
        Bxy{i}.dx = problem.ZeroTxM;
        Bxy{i}.dy = By{i-dimM};
    end

    % Under this basis Bxy, the following codes return a saddle-point
    % linear system whose matrix has the form
    %
    %
    %           [ HessLag_mat + THETA_mat | Hx_mat]
    %   T_mat = -----------------------------------
    %           [ Hx_mat'                 | O     ]
    %
    % where
    % - n:= dimM, l:=dimEuc_eq,
    % - HessLag_mat and THETA_mat are symmetric and n x n,
    % - Hx_mat is l x n, with l <= n,
    % - O is l x l zero matrix.

    % the next code is equal to run:
    % HessLag_mat = SelfAdj_operator2matrix(M, x, x, HessLag, Bx, Bx);
    % THETA_mat = SelfAdj_operator2matrix(M, x, x, THETA, Bx, Bx);
    % Aw_mat = HessLag_mat + THETA_mat;
    Aw_mat = SelfAdj_operator2matrix(M, x, Aw, Bx);

    % the next code is equal to run:
    % Hx_mat = operator2matrix(problem.Euc_eq, y, x, Hx, By, Bx, M);
    % T_mat =[Aw_mat,  Hx_mat; ...
    %        Hx_mat', zeros(dimEuc_eq, dimEuc_eq)];
    % But, Hxaj is cheaper, because Hx needs to call
    % orthogonal projection onto tangent space.
    Hxaj_mat = operator2matrix(M, x, y, Hxaj, Bx, By, problem.Euc_eq); % Hx_mat'=Hxaj_mat.
    T_mat =[Aw_mat,  Hxaj_mat'; ...
        Hxaj_mat, zeros(dimEuc_eq, dimEuc_eq)];

    cq_vec = tangent2vec(M_Euc_eq, xy, Bxy, b);

    % direct method
    sol_vec = linsolve(T_mat, cq_vec, struct('SYM', true));
    NTdir= lincomb(M_Euc_eq, xy, Bxy, sol_vec);
    RepresentMat = T_mat;
    RepresentMatOrder = dimProd;
end

end


%% local funs
function [A, Bx] = SelfAdj_operator2matrix(M, x, F, Bx)

% Forms a symmetric matrix representing a linear self-adjoint operator F
% to and from the tangent space T_x M.

n = numel(Bx);
A = zeros(n, n);

% Generate the upper right part.
for j = 1 : n
    FBxj = F(Bx{j});
    for i = 1 : j
        A(i, j) = M.inner(x, FBxj, Bx{i});
    end
end

% flip half of matrix over the diagonal to make a symmetric matrix
A = A + triu(A,1)';

end










