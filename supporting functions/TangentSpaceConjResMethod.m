function [vfinal, t, rel_res, info] = TangentSpaceConjResMethod(Aw, T, b, v0, M, x, tol, maxiter, condet)
% Conjugate Residual Method for solving linear operator euqation: A(v)=b,
% where A is some self-adjoint operator to and from some linear space E,
% b is an element in E. Assume the existence of solution v.

% Yousef Saad - Iterative Methods for Sparse Linear Systems,
% Second Edition-SIAM (2003) P203. ALGORITHM 6.20

%% Case 1. has_eq_cost = false.

if ~condet.has_eq_cost
    A = Aw; % In this case, A is an operator Aw to and from TxM.
    v = v0.dx; % initialization
    r = M.lincomb(x, 1, b.dx, -1, A(v)); % r are residuals.
    p = r; % p are conjugate directions.
    b_norm = M.norm(x, b.dx);
    r_norm = M.norm(x, r);
    rel_res = r_norm/b_norm;
    Ar = A(r);
    Ap = A(p);
    rAr = M.inner(x, r, Ar);
    t = 0; % at t-th iteration
    info = zeros(maxiter+1,2);
    while 1
        info(t+1,:) = [t, rel_res];
        t = t + 1;
        a = rAr/(M.inner(x, Ap, Ap)); % step length
        v = M.lincomb(x, 1, v, a, p); % update x      % v + a*p
        r = M.lincomb(x, 1, r, -a, Ap); % residual      % r - a*Ap
        r_norm = M.norm(x, r);
        rel_res = r_norm/b_norm;
        if rel_res < tol || t == maxiter
            break
        end
        %r = M.proj(x,r);
        Ar = A(r);
        old_rAr = rAr;
        rAr = M.inner(x, r, Ar);
        beta = rAr/old_rAr; % improvement this step
        p = M.lincomb(x, 1, r, beta, p); % search direction % r + beta*p
        Ap = M.lincomb(x, 1, Ar, beta, Ap);   % Ar + beta*Ap
    end
    vfinal.dx = v;
    vfinal.dy = [];
end


%% Case 2. has_eq_cost = ture.

if condet.has_eq_cost
    A = T; % In this case, A is an operator T to and from TxM x Euc_eq.
    v = v0; % initialization
    r = ProdLincomb(M, x, 1, b, -1, A(v)); % r are residuals.
    p = r; % p are conjugate directions.
    b_norm = ProdNorm(M, x, b);
    r_norm = ProdNorm(M, x, r);
    rel_res = r_norm/b_norm;
    Ar = A(r);
    Ap = A(p);
    rAr = ProdInner(M, x, r, Ar);
    t = 0; % at t-th iteration
    info = zeros(maxiter+1,2);
    while 1
        info(t+1,:) = [t, rel_res];
        t = t + 1;
        a = rAr/(ProdInner(M, x, Ap, Ap)); % step length
        v = ProdLincomb(M, x, 1, v, a, p); % update x % v + a*p
        r = ProdLincomb(M, x, 1, r, -a, Ap); % residual % r - a*Ap
        r_norm = ProdNorm(M, x, r);
        rel_res = r_norm/b_norm;
        if rel_res < tol || t == maxiter
            break
        end
        %r = M.proj(x,r);
        Ar = A(r);
        old_rAr = rAr;
        rAr = ProdInner(M, x, r, Ar);
        beta = rAr/old_rAr; % improvement this step
        p = ProdLincomb(M, x, 1, r, beta, p); % search direction % r + beta*p
        Ap = ProdLincomb(M, x, 1, Ar, beta, Ap);   % Ar + beta*Ap
    end
    vfinal=v;
end

end







