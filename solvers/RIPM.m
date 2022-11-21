function [xfinal, costfinal, residual, info, options] = RIPM(problem, x0, options)
% Riemannian Interior Point Methods solver.
%
% function [xfinal, costfinal, residual, info, options] = RIPM(problem)
% function [xfinal, costfinal, residual, info, options] = RIPM(problem, x0)
% function [xfinal, costfinal, residual, info, options] = RIPM(problem, x0, options)
% function [xfinal, costfinal, residual, info, options] = RIPM(problem, [], options)
%
% This is a primal-dual interior point methods solver for nonlinear
% optimization problems on Riemannian manifolds, which aims to minimize the
% cost function in the given problem structure with (in)equality
% constraints.
%
% It requires access to the gradient and the Hessian of the cost function
% and the constraints.
%
% The initial iterate is x0 if it is provided. Otherwise, a random point on
% the manifold is picked. To specify options whilst not specifying an
% initial iterate, give x0 as [] (the empty matrix).
%
% For a description of the algorithm and theorems offering convergence
% guarantees, see the references below:
%   Z. Lai and A. Yoshise. Riemannian Interior Point Methods for
%   Constrained Optimization on Manifolds,
%   https://arxiv.org/abs/2203.09762
%   An earlier version of this article has been circulated under the title
%   ``Superlinear and Quadratic Convergence of Riemannian Interior Point
%   Methods.''
%
%% Input to the problem structure
% In addition, we have a separate document with detailed instructions for
% the implementation. The documentation is available on github.
%
% * problem.M
%
% * problem.cost = @(x) % f(x)
% * problem.egrad = @(x) % egrad_{x} f(x)
% * problem.ehess = @(x, dx) % ehess_{x} f(x)
%
% Lagrangian = f(x) + <z,g(x)> + <y,h(x)>. Here, we recognize two
% real-valued functions on M: MAP_ineq: x maps to <z,g(x)> by fixing z;
% MAP_eq: x maps to <y,h(x)> by fixing y.
%
% Inequality constraints is necessary.
% * problem.Euc_ineq % Codomain of g(x), manifold structure, multiplier z
% also is in Euc_ineq.
% * problem.ineq_con = @(x) % g(x): M to Euc_ineq.
% * problem.barGx = @(x, z) % egrad_{x} MAP_ineq(x)
% * problem.barGxaj = @(x, dx) % Adjoint of problem.barGx; if we restrict
% domain to TxM, then Gxaj = barGxaj.
% * problem.ehess_barGx = @(x, z, dx) % ehess_{x} MAP_ineq(x)
%
% Equality constraints is optional.
% * problem.Euc_eq % Codomain of h(x), manifold structure,  multiplier y
% also is in Euc_eq.
% * problem.eq_con = @(x) % h(x): M to Euc_eq.
% * problem.barHx = @(x, y) % egrad_{x} MAP_eq
% * problem.barHxaj = @(x, dx) % Adjoint of problem.barHx; if we restrict
% domain to TxM, then Hxaj = barHxaj.
% * problem.ehess_barHx = @(x, y, dx) % ehess_{x} MAP_eq
%
% Be cautious: Equality constraints in matrix form often contain duplicate
% expressions, e.g., h(x)=X'*X-I, when h(x)=0, only the diagonal and above
% component functions are independent and valid constraints, while the rest
% is redundant. To avoid this problem, it is sufficient to set
% problem.Euc_eq to be symmetricfactory(n). problem.Euc_eq is not only the
% domain of h(x), but also its dimmension is the actual number of equality
% constraints.
%
%% Outputs
% The outputs 'xfinal', 'costfinal', 'residual' are the last reached point
% on the manifold, its cost, and its KKT residual. The output 'info' is a
% struct-array which contains information about the iterations:
%   iter (integer)
%       The iteration number, i.e., number of steps considered so far. The
%       initial guess is 0.
%   xcurrent
%       The current point x.
%   xCurCost (double)
%       The corresponding cost value.
%   xCurPhi (double)
%       The corresponding merit value.
%   KKT_residual (double)
%       The corresponding residual of the KKT conditions.
%   sigma (double)
%       The barrier parameter for the perturbed Newton equation.
%   rho (double)
%       The barrier parameter for the perturbed Newton equation.
%   time (double)
%       The total elapsed time in seconds to reach the corresponding cost.
%   stepsize (double, <=1)
%       The size of the steplength determined by the backtracking
%       linesearch.
%   ls_iters (integer)
%       The number of loops executed by line search.
%   ls_max_steps_break (boolean)
%       Whether the line search ends due to the excess of the maximal
%       number of the backtracking (ls_max_steps):
%           0: the backtracking ends normally. 1: the backtracking ends due
%           to the excess of the number.
%   dist (double)
%       The (Riemannian) distance between the previous and the new
%       iterates.
%
% For example, type [info.KKT_residual] to obtain a vector of the
% successive the KKT residual reached at each iteration.
%
% The output 'info' will contains more informations if
% options.checkNTequation (boolean) is true.
%   NTdir_error1 (double)
%       The residual the non-condensed Newton equation. It should be zero.
%   NTdir_error2 (double)
%       If computed Newton dirction is correct, then NTdir_error2 should be
%       zero.
%   NTdir_norm (double)
%       The norm of Newton dirction.
%   NTdir_angle (double)
%       The angle between -grad phi and Newton dirction.
%
%% options
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%   maxiter (500)
%       The algorithm terminates if maxiter iterations have been executed.
%   maxtime (3600)
%       The algorithm terminates if maxtime seconds elapsed.
%   tolKKTres (1e-6)
%       The algorithm terminates if the KKT residual drops below this.
%   KrylovIterMethod (false)
%       How to solve the condensed Newton equation:
%           0: Use representing matrix method. 
%           1: Use Krylov subspace iterative method, in particular,
%           we use Conjugate Residual (CR) Method. 
%           See TangentSpaceConjResMethod.m
%   KrylovTolrelres (1e-9)
%       CR method terminates when relative residual drops below this.
%   KrylovMaxiter (1000)
%       CR method terminates if maxiter iterations have been executed.
%   checkNTequation (false)
%       Whether to check the Newton equation. If true, more informations
%       will be contained in output 'info'.
%   gamma (0.9)
%       Initial parameter for central functions. The value must be greater
%       than 0.5 and less than 1.
%   ls_execute_fun2 (false)
%       Whether to excute the backtracking for seconde central function.
%   ls_beta (1e-4)
%       Constant for the backtracking line search. The value must belong to
%       the interval (0, 1).
%   ls_theta (0.5)
%       Reduction rate of the backtracking to find an appropriate
%       steplength. The value must belong to the interval (0, 1).
%   ls_max_steps (50)
%       The algorithm breaks the backtracking if the ls_max_steps trial
%       have been executed.
%   y (zeros)
%       Initial Lagrange multiplier vector for equlity constraints.
%   z (rand)
%       Initial Lagrange multiplier vector for inequlity constraints.
%   s (rand)
%       Initial slack vector.
%   heuristic_z_s (false)
%       Whether to use a heuristic initial z and s.
%           0: Random initial z and s.
%           1: A heuristic way to construct z and s.
%   desired_tau_1 (.5)
%       Construct the particular initial z,s such that tau1 equals
%       desired_tau_1. The value must belong to the interval (0, 1).
%       desired_tau_1 is used only when heuristic_z_s is true.
%   important (1e0)
%       important controls tau_2 and initial rho. important is used only
%       when heuristic_z_s is true.
%   verbosity (3)
%       Integer number used to tune the amount and the frequency of output
%       the algorithm generates during execution (mostly as text in the
%       command window). The higher, the more output. 0 means silent.
%   rankviopena (1e+8)
%       If considering optimization on fixed-rank manifolds, the value is
%       used as the penalty when violating the fixed-rank constraints.
% 
% 
% 
% 
% Original author: Zhijian Lai, November, 20, 2022. Contributors: Change
% log:
%           November, 20, 2022: Write the code.
% 
% 
% 

%% Set localdefaults, a struct to be combined with argument options for
% declaring hyperparameters.

% ...\Manopt_7.0\manopt\manopt\tools\multitransp.m Line 36 has been
% changed!!!!

condet = RIPM_constraintsdetail(problem);

if ~condet.has_ineq_cost
    error('RIPM requires that the problem must contain inequality constraints.');
end

% Stopping criteria
localdefaults.maxiter = 500;
localdefaults.maxtime = 3600; % sec.
localdefaults.tolKKTres = 1e-6;
% Solve Newton equation
localdefaults.KrylovIterMethod = false;
localdefaults.KrylovTolrelres = 1e-9;
localdefaults.KrylovMaxiter = 1000;
localdefaults.checkNTequation = false;
% Line search
localdefaults.gamma = .9;
localdefaults.ls_execute_fun2 = false;
localdefaults.ls_beta = 1e-4;
localdefaults.ls_theta = 0.5;
localdefaults.ls_max_steps = 50;
% Other parameters
localdefaults.heuristic_z_s = false; % initial z and s
localdefaults.desired_tau_1 = .5; % we let tau_1 be desired_tau_1 in (0,1).
localdefaults.important = 1e0;% important controls tau_2 and initial rho 有时候如果不收束, 就让important小一点
% Display
localdefaults.verbosity = 3;
% Only when using fixed-rank manifolds
localdefaults.rankviopena = 1e+8;

% Merge global and local defaults, then merge user options, if any.
localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
options = mergeOptions(localdefaults, options);

%% Various definitions.

M = problem.M;

Gx = @(x, z) M.egrad2rgrad(x, problem.barGx(x, z));
Gxaj = @problem.barGxaj;

if condet.has_eq_cost
    Hx = @(x, y) M.egrad2rgrad(x, problem.barHx(x, y));
    Hxaj = @problem.barHxaj;
    egradLagrangian = @(x, y, z) problem.egrad(x) ...
        + problem.barHx(x,y) ...
        + problem.barGx(x,z);
    ehessLagrangian = @(x, y, z, dx) problem.ehess(x, dx) ...
        + problem.ehess_barHx(x, y, dx) ...
        + problem.ehess_barGx(x, z, dx);
else
    egradLagrangian = @(x, y, z) problem.egrad(x) + problem.barGx(x,z);
    ehessLagrangian = @(x, y, z, dx) problem.ehess(x, dx) + problem.ehess_barGx(x, z, dx);
end

gradLagrangian = @(x, y, z) M.egrad2rgrad(x, egradLagrangian(x, y, z));

hessLagrangian = @(x, y, z, dx) M.ehess2rhess(x, egradLagrangian(x, y, z), ...
    ehessLagrangian(x, y, z, dx), dx);

% KKT vector field
    function F = KKTVecField(x, y, z, s)
        F.dx = gradLagrangian(x, y, z);
        if condet.has_eq_cost
            F.dy = problem.eq_con(x);
        else
            F.dy = [];
        end
        F.dz = problem.ineq_con(x) + s;
        F.ds = z.*s;
    end

%% Initialization

% If no initial point x0 is given by the user, generate one at random. But
% we strongly recommend using a strictly feasible point,i.e, g(x)<0, if any.
if ~exist('x0', 'var') || isempty(x0)
    xCur = M.rand();
else
    xCur = x0;
end

% Initial y. In problems constrained only inequalities constrains, y is
% still preserved but empty.
if condet.has_eq_cost
    eq_size = size(problem.eq_con(M.rand()));
    y = zeros(eq_size);
else
    y = [];
end

% Initial z and s.
ineq_num = condet.n_ineq_constraint_cost;
ineq_size = size(problem.ineq_con(M.rand()));
if options.heuristic_z_s
    z  = ones(ineq_size);
    z(1,1) = realsqrt((ineq_num - 1)/(ineq_num/options.desired_tau_1 - 1));
    s = options.important*z;
else
    z  = rand(ineq_size);
    s  = rand(ineq_size);
end

% Construct constant perturbed item Ehat.
problem.dimM = M.dim();
problem.ZeroTxM = M.zerovec();
if condet.has_eq_cost
    problem.dimEuc_eq = problem.Euc_eq.dim();
    problem.ZeroTyEuc = problem.Euc_eq.zerovec(); % Do not use zeros(eq_size) because of fixed rank minifold.
else
    problem.dimEuc_eq = [];
    problem.ZeroTyEuc = [];
end
problem.ZeroTzEuc= problem.Euc_ineq.zerovec();
ehat = ones(ineq_size);
Ehat = struct('dx', problem.ZeroTxM, 'dy', problem.ZeroTyEuc, 'dz', problem.ZeroTzEuc, 'ds', ehat);

% Initial point for Krylov Subspace Iterative Method.
v0 = struct('dx', problem.ZeroTxM, 'dy', problem.ZeroTyEuc);

% Only for calculating the KKT residual on fixed-rank manifolds
if contains(M.name(),'rank')
    if isfield(options, 'rank')
        rankval = options.rank;
    else
        tmpx = M.rand();
        rankval = rank(tmpx.S);
    end
end

% Declare some variables for the initial savestats.
iter = 0;
KKTvec = KKTVecField(xCur, y, z, s);
PhiCur = ProdNorm(M, xCur, KKTvec)^2;
CostCur = getCost(problem, xCur);
KKTResidualCur = KKT_residual(xCur, y, z);
timetic = tic();

% Set constants for centrality conditions
tau_1 = min(z.*s, [], 'all')/(F_inner(z, s)/ineq_num);
tau_2 = F_inner(z, s)/realsqrt(PhiCur);

% Construct parameters sigma to controls the final convergence rate.
sigma = min(.5, realsqrt(PhiCur)^.5);
rho = F_inner(z, s)/ineq_num;
gamma = options.gamma;

% Save stats in a struct array info, and preallocation.
stats = savestats();
info(1) = stats;
info(min(10000, options.maxiter+1)).iter = [];

% Stopping flag, finally it should be true
stop = false;
totaltime = tic();

% Main loop
while true
    % Display iteration information.
    if options.verbosity >= 2
        fprintf('Iter: %d, Cost: %f, phi: %.16e, KKT residual: %.16e \n', iter, CostCur, PhiCur, KKTResidualCur);
    elseif options.verbosity >= 1
        if mod(iter, 100) == 0 && iter ~= 0
            fprintf('Iter: %d, Cost: %f, phi: %.16e, KKT residual: %.16e \n', iter, CostCur, PhiCur, KKTResidualCur);
        end
    end

    % Start timing this iteration.
    iter = iter + 1;
    timetic = tic();

    %% Get Newton direction NTdict by solving condensed NT equation.

    % Right hand side (cq) of condensed NT equation is (c, q).
    cq.dx = M.lincomb(xCur, -1, KKTvec.dx, ...
        -1, Gx(xCur, (z.*KKTvec.dz + sigma*rho*ehat - KKTvec.ds)./s));
    cq.dy = - KKTvec.dy; % if no eq_con, KKTvec.dy=[].

    % Define some operators.
    OpratorHessLag = @(dx) hessLagrangian(xCur, y, z, dx);
    OpratorTHETA = @(dx) Gx(xCur, Gxaj(xCur, dx).*(z./s));
    OpratorAw = @(dx) M.lincomb(xCur, 1, OpratorHessLag(dx), 1, OpratorTHETA(dx));
    OpratorHx = @(dy) Hx(xCur, dy);
    OpratorHxaj = @(dx) Hxaj(xCur, dx);
    OpratorT = @(dxdy) struct( ...
        'dx', M.lincomb(xCur, 1, OpratorAw(dxdy.dx), 1, OpratorHx(dxdy.dy)), ...
        'dy', OpratorHxaj(dxdy.dx));

    % Solve condensed NT equation T(dx,dy) = (c, q).
    if options.KrylovIterMethod
        NTdir = TangentSpaceConjResMethod( ...
            OpratorAw, ...
            OpratorT, ...
            cq, v0, M, xCur, options.KrylovTolrelres, options.KrylovMaxiter, condet);
    else % here, we use the general but not efficient method
        NTdir = RepresentMatMethod( ...
            OpratorAw, ...
            OpratorHxaj, ...
            cq, M, xCur, y, condet, problem);
    end

    % Rcovery dz and ds.
    NTdir.dz = (z.*(Gxaj(xCur, NTdir.dx)+KKTvec.dz) + sigma*rho*ehat - KKTvec.ds)./s;
    NTdir.ds = (sigma*rho*ehat - KKTvec.ds - s.*NTdir.dz)./z;

    %% Check Newton direction NTdir.
    % DEBUG ONLY.
    if options.checkNTequation
        nablaF = @(dw) CovarDerivKKT(xCur, y, z, s, dw);
        nablaFaj = @(dw) CovarDerivKKTaj(xCur, y, z, s, dw);

        % Check Item #1: The residual the non-condensed Newton equation.
        % Note that the right hand side of original NT euqation is
        % -KKTvec+sigma*rho*Ehat.
        % NTdir_error1 should be zero.
        NTeq_rhs = ProdLincomb(M, xCur, -1, KKTvec, sigma*rho, Ehat);
        nablaF_NTdir = nablaF(NTdir);
        NTdir_error1 = ProdDist(M, xCur, nablaF_NTdir, NTeq_rhs);

        % Check Item #2: If NTdir is correct solution, then
        % <grad phi, NTdir> = 2(|F(w)|^{2}+sigma*rho*InnerProduct(z,s))
        % holds, where grad phi = 2*nablaFaj(KKTvec).
        % NTdir_error2 should be zero.
        gradphi = ProdLincomb(M, xCur, 2, nablaFaj(KKTvec));
        val_innerprodut = ProdInner(M, xCur, gradphi, NTdir); % <grad phi,NTdir>
        NTdir_error2 = abs(val_innerprodut - 2*(sigma*rho*F_inner(z,s)-PhiCur));

        % Record Item: record norm of NTdir; angle between - grad phi and NTdir.
        Norm_gradphi = ProdNorm(M, xCur, gradphi);
        NTdir_norm =  ProdNorm(M, xCur, NTdir);
        NTdir_angle = - val_innerprodut/(Norm_gradphi * NTdir_norm);
    end

    %% Backtracking line search and update.

    % Central functions
    fun_1 = @(z, s) min(z.*s, [], "all") - gamma*tau_1*(F_inner(z, s)/ineq_num);
    fun_2 = @(z, s, Phi) F_inner(z, s) - gamma*tau_2*realsqrt(Phi);

    % Note that <grad phi, NTdir> = ls_RightItem, if NTdir is a correct solution.
    ls_RightItem = 2*( sigma*rho*F_inner(z,s) - PhiCur);

    stepsize=1;
    ls_max_steps_flag = false;
    r=0; % backtracking counter
    while 1
        xNew = M.retr(xCur, NTdir.dx, stepsize);
        yNew = y + stepsize*NTdir.dy;
        zNew = z + stepsize*NTdir.dz;
        sNew = s + stepsize*NTdir.ds;
        KKTvec = KKTVecField(xNew, yNew, zNew, sNew);
        PhiNew = ProdNorm(M, xNew, KKTvec)^2;
        if PhiNew - PhiCur <= options.ls_beta*stepsize*ls_RightItem ...
                && fun_1(zNew, sNew) >= 0 ...
                && (~options.ls_execute_fun2 || fun_2(zNew, sNew, PhiNew) >= 0)
            break
        end
        r = r+1;
        if r > options.ls_max_steps
            ls_max_steps_flag = true;
            break;
        end
        stepsize = stepsize*options.ls_theta;
    end

    % Update information for savestats
    if contains(M.name(),'Stiefel')
        dist = norm(xCur - xNew, 'fro');
    elseif contains(M.name(),'rank')
        % Only assuming for 'fixedrankembeddedfactory'
        if ~exist('xCurmat', 'var')
            xCurmat = xCur.U * xCur.S * xCur.V';
        end
        newxmat = xNew.U * xNew.S * xNew.V';
        dist = norm(xCurmat - newxmat, 'fro');
        xCurmat = newxmat;
    else
        dist = M.dist(xCur, xNew);
    end

    xCur = xNew; y = yNew; z = zNew; s = sNew;
    PhiCur = PhiNew;
    CostCur = getCost(problem, xCur);
    KKTResidualCur = KKT_residual(xCur, y, z);

    % update parameters
    sigma = min(.5, realsqrt(PhiCur)^.5);
    rho = F_inner(z,s)/ineq_num;
    gamma = 0.5*(gamma + 0.5);

    % Savestats
    stats = savestats();
    info(iter+1) = stats;

    % Refer to stopping criteria
    if iter >= options.maxiter
        fprintf('Max iter %d count reached\n', options.maxiter);
        options.reason = "Max iter count reached";
        stop = true;
    elseif toc(totaltime) >= options.maxtime
        fprintf('Max time %f [s] exceeded\n', options.maxtime);
        options.reason = "Max time exceeded";
        stop = true;
    elseif KKTResidualCur <= options.tolKKTres
        fprintf('KKT Residual %e tolerance reached\n', options.tolKKTres);
        options.reason = "KKT Residual tolerance reached";
        stop = true;
    elseif dist == 0
        % which will never occur because of the Armijo rule
        fprintf('Any parameter was not updated');
        options.reason = 'Any parameter was not updated';
        stop = true;
    elseif isnan(PhiCur)
        fprintf('KrylovIterMethod cannot deal with the ill-posed Newton equation.');
        options.reason = 'KrylovIterMethod cannot deal with the ill-posed Newton equation.';
        stop = true;
    end

    if stop
        options.totaltime = toc(totaltime);
        break
    end
end

% Retrun final results.
xfinal = xCur;
costfinal = CostCur;
residual  = KKTResidualCur;

info = info(1:iter+1); % Crop blank of info.

if options.verbosity >= 1
    fprintf('Total time is %f [s] (excludes statsfun)\n', info(end).time);
end

%% Sub funtions
    function stats = savestats()
        stats.iter = iter;
        stats.xcurrent = xCur;
        stats.xCurCost = CostCur;
        stats.xCurPhi = PhiCur;
        stats.KKT_residual = KKTResidualCur;
        stats.sigma=sigma;
        stats.rho=rho;
        if iter == 0
            stats.time = toc(timetic);
            stats.stepsize = NaN;
            stats.ls_iters = NaN;
            stats.ls_max_steps_flag = NaN;
            stats.dist =  NaN;
        else
            stats.time = info(iter).time + toc(timetic);
            stats.stepsize = stepsize;
            stats.ls_iters = r;
            stats.ls_max_steps_flag = ls_max_steps_flag;
            stats.dist = dist;
        end
        if options.checkNTequation && iter == 0
            stats.NTdir_error1 = NaN;
            stats.NTdir_error2 = NaN;
            stats.NTdir_norm = NaN;
            stats.NTdir_angle = NaN;
        elseif options.checkNTequation
            stats.NTdir_error1 = NTdir_error1;
            stats.NTdir_error2 = NTdir_error2;
            stats.NTdir_norm = NTdir_norm;
            stats.NTdir_angle = NTdir_angle;
        end
        stats = applyStatsfun(problem, xCur, [], [], options, stats);
    end

% Covariant Derivative of KKT vector field
% DEBUG ONLY.
    function nablaF = CovarDerivKKT(x, y, z, s, dw)
        nablaF.dx = M.lincomb(xCur, 1, hessLagrangian(x, y, z, dw.dx), 1, Gx(x, dw.dz));
        if condet.has_eq_cost
            nablaF.dx = M.lincomb(xCur, 1, nablaF.dx, 1, Hx(x, dw.dy));
            nablaF.dy = Hxaj(x, dw.dx);
        else
            nablaF.dy = [];
        end
        nablaF.dz = Gxaj(x, dw.dx) + dw.ds;
        nablaF.ds = z.*dw.ds + s.*dw.dz;
    end

% Adjoint of Covariant Derivative of vector field
% DEBUG ONLY.
    function nablaFaj = CovarDerivKKTaj(x, y, z, s, dw)
        nablaFaj.dx = M.lincomb(xCur, 1, hessLagrangian(x, y, z, dw.dx), 1, Gx(x, dw.dz));
        if condet.has_eq_cost
            nablaFaj.dx = M.lincomb(xCur, 1, nablaFaj.dx, 1, Hx(x, dw.dy));
            nablaFaj.dy = Hxaj(x, dw.dx);
        else
            nablaFaj.dy = [];
        end
        nablaFaj.dz = Gxaj(x, dw.dx) + s.*dw.ds;
        nablaFaj.ds = z.*dw.ds + dw.dz;
    end

% Calculating the KKT residual at (XCur, mus, lambdas)
% Note that mus = z, and lambdas = y here.
    function val = KKT_residual(xCur, y, z)
        % Test the violation of Riemannian KKT conditions.
        % Note that the square root of KKT_residual approximates to Phi.

        % Manifold violation: x in M.
        manvio = manifoldViolation(xCur);

        % Stationarity violation: grad_{x} Lag = 0.
        xGrad = gradLagrangian(xCur, y, z);
        xGradpowvio = M.norm(xCur, xGrad)^2;

        % Complementary slackness violation: g(x).*z = 0.
        violation = z.*problem.ineq_con(xCur);
        compowvio = norm(violation,"fro")^2;

        % #Unnecessary because central condtion.
        % % Dual feasibility violation: z >= 0.
        %violation = min(z, 0);
        %muspowvio = norm(violation,"fro")^2;

        % Primal feasibility inequality constraints violation: h(x) = 0; g(x) <=0.
        violation = max(problem.ineq_con(xCur), 0);
        ineqpowvio = norm(violation,"fro")^2;
        % equality constraints violation, if any.
        if condet.has_eq_cost
            violation = problem.eq_con(xCur);
            eqpowvio = norm(violation,"fro")^2;
        else
            eqpowvio = 0;
        end

        val = xGradpowvio + compowvio + ineqpowvio + eqpowvio;
        val = realsqrt(val) + manvio;
    end

% manifold violation
    function manvio = manifoldViolation(xCur) % # No square!
        % According to the type of manifold, calculate the violation from
        % constraints seen as the manifold. (NOTE: the procedure should be
        % out of the solver because the calculation depends on the choice
        % of the manifold.)
        manvio = 0;
        if contains(M.name(),'Sphere')
            p = xCur(:);
            manvio = abs(p.'*p - 1);
        elseif contains(M.name(),'Oblique')
            [~,N] = size(xCur);
            for col = 1:N
                manvio = manvio + abs(xCur(:,col).' * xCur(:,col) - 1);
            end
        elseif contains(M.name(),'Stiefel') %% #
            [~,k] = size(xCur);
            violation = triu(xCur.'*xCur-eye(k));
            manvio = norm(violation(:),1);
        elseif contains(M.name(),'rank') %% #
            try
                xCurrank = rank(xCur.S);
            catch
                xCurrank = rankval; % for NaN matrix
            end
            if xCurrank ~= rankval
                manvio = options.rankviopena; % as if Inf;
            end
        end
    end
end
