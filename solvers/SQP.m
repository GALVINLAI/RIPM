function [xfinal, costfinal, residual, info, options] = SQP(problem0, x0, options)
% Sequential Quadratic Programming solver for smooth objective functions 
% on Riemannian manifolds.
%
% function [xfinal, costfinal, residual, info, options] = sqponmani(problem0)
% function [xfinal, costfinal, residual, info, options] = sqponmani(problem0, x0)
% function [xfinal, costfinal, residual, info, options] = sqponmani(problem0, x0, options)
% function [xfinal, costfinal, residual, info, options] = sqponmani(problem0, [], options)
%
% This is a Sequential Qudratic Programming solver for nonlinear optimization problems
% on Riemannian manifolds, which aims to minimize the cost function
% in the given problem structure with (in)equality constraints.
%
% It requires access to the gradient and the Hessian of the cost function
% and the constraints.
%
% For a description of the algorithm and theorems offering convergence
% guarantees, see the references below:
%   M.Obara, T. Okuno, and A. Takeda. Sequential quadratic optimization for nonlinear optimization problems
%   on Riemannian manifolds, https://arxiv.org/abs/2009.07153.
%
% The initial iterate is x0 if it is provided. Otherwise, a random point on
% the manifold is picked. To specify options whilst not specifying an
% initial iterate, give x0 as [] (the empty matrix).
%
% The two outputs 'xfinal', 'costfinal', 'residual' are the last reached point on the manifold,
% its cost, and its KKT residual.
% 
% The output 'info' is a struct-array which contains information about the
% iterations:
%   iter (integer)
%       The (outer) iteration number, i.e., number of steps considered
%       so far. The initial guess is 0.
%   cost (double)
%       The corresponding cost value
%   gradnorm (double)
%       The (Riemannian) norm of the gradient of the Lagrangian
%   time (double)
%       The total elapsed time in seconds to reach the corresponding cost
%   stepsize (double, <=1)
%       The size of the steplength determined by the backtracking
%       linesearch
%   ls_max_steps_break (boolean)
%       Whether the linesearch ends due to the excess of the maximal number
%       of the backtracking (ls_max_steps):
%           0: the backtracking ends normally 
%           1: the backtracking ends due to the excess of the number
%   dist (double)
%       The (Riemannian) distance between the previous and the new iterates
%   qpexitflag (integer)
%       Reason quadprog stopped, returned as an integer. It should be 1:
%           1: Function converged to the solution x.
%           0: Number of iterations exceeded options.MaxIterations.
%           -2: Problem is infeasible. Or, for 'interior-point-convex', the
%               step size was smaller than options.StepTolerance, but 
%               constraints were not satisfied.
%           2: Step size was smaller than options.StepTolerance, 
%               constraints were satisfied. (only when using
%               'interior-point-convex')
%           -6: Nonconvex problem detected. (only when using 'interior-point-convex')
%           -8: Unable to compute a step direction. (only when 'interior-point-convex')
%       (Cf.: https://www.mathworks.com/help/optim/ug/quadprog.html#d123e131181)
%   rho (double)
%       The penalty parameter for the ell-1 penalty function
%   KKT_residual (double)
%       The sum of the residual of the KKT conditions
%   maxviolation (double)
%       The maximal value of the violations of (in)equality constraints
%   meanviolation (double)
%       The mean value of the violations of (in)equality constraints
%
% For example, type [info.gradnorm] to obtain a vector of the successive
% the norms of the gradient of the Lagrangian reached at each iteration.
%
% The options structure is used to overwrite the default values. All
% options have a default value and are hence optional. To force an option
% value, pass an options structure with a field options.optionname, where
% optionname is one of the following and the default value is indicated
% between parentheses:
%   maxiter (300)
%       The algorithm terminates if maxiter (outer) iterations have been
%       executed.
%   maxtime (3600)
%       The algorithm terminates if maxtime seconds elapsed.
%   tolKKTres (1e-8)
%       The algorithm terminates if the KKT residual drops below this.
%   modify_hessian ('mineigval_matlab')
%       The algorithm sets the linear operators on subproblems in the
%       following manner:
%           'eye': the identity matrix,
%           'mineigval_matlab': a regularized Hessian of the Lagrangian.
%           We regularize the Hessian matrix if it is not positive-definite
%           by replacing negative engenvalues with positive epsilon,
%           'mineigval_manopt' (default): a regularized Hessian of the Lagrangian.
%           We regularize the Hessian matrix if it is not positive-definite
%           by adding the absolute value of the minimal eigenvalue and an 
%           positive epsilon to the diagonal elements in the matrix. We get
%           the eigenvalue by solving the Rayleigh quotient minization.
%   mineigval_correction (1e-8)
%       The algorithm adds this value for the linear matrix on the
%       subproblem to be positive-definite when using 'mineigval_matlab' or
%       'mineigval_manopt' procedures.
%   mineigval_threshold (1e-3)
%       Threshold value used in correcting the diagonal
%       elements by 'mineigval_manopt'. If the corrected elements is lower
%       than this value, we further modify the value to attain the
%       threshold.
%   tau (0.5)
%       Constant for updating the penalty parameter. The value
%       must be positive.
%   rho (1)
%       Initial penarty parameter for the ell-1 penalty
%       function. The value must be positive.
%   beta (0.9)
%       Magnification of the backtracking to find an appropriate
%       steplength. The value must belong to the interval (0, 1).
%   gamma (0.25)
%       Constant for the backtracking line search. The value
%       must belong to the interval (0, 1).
%   mus (ones)
%       Initial Lagrange multiplier vector for inequlity constraints.
%   lambdas (ones)
%       Initial Lagrange multiplier vector for equlity constraints.
%   ls_max_steps (10000)
%       The algorithm breaks the backtracking if the ls_max_steps trial
%       have been executed.
%   ls_threshold (1e-8)
%       The algorithm finishes the line search if a steplength satisfies 
%       the condition with the tolerance of this threshold.
%   verbosity (1)
%       Integer number used to tune the amount and the frequency of output 
%       the algorithm generates during execution (mostly as text in the 
%       command window). The higher, the more output. 0 means silent.
%   qp_verbosity (0)
%       Integer number used to tune the amount and the frequency of output 
%       the algorithm generates during execution of the quadratic optimization
%       subproblems (mostly as text in the command window). The higher, the
%       more output. 0 means silent.
%   --
%   rankviopena (1e+8)
%       If considering optimization on fixed-rank manifolds, the value is
%       used as the penalty when violating the fixed-rank constraints.
%       Note: the procedures to check the satisfiablity of manifold constraints 
%           should be moved to outside of the solver.
%
% Original author: Mitsuaki Obara, January 20, 2020.
% Contributors: 
% Change log: 
%           June, 7, 2021: Clean the code (just delete the comments) and
%           add the new introduction.
%           January, 20, 2020: Write the code.

    % % Verify that the problem description is sufficient for the solver.
    % if ~canGetCost(problem0)
    %     warning('manopt:getCost', ...
    %         'No cost provided. The algorithm will likely abort.');
    % end
    % if ~canGetGradient(problem0) && ~canGetApproxGradient(problem0)
    %     % Note: we do not give a warning if an approximate gradient is
    %     % explicitly given in the problem description, as in that case the user
    %     % seems to be aware of the issue.
    %     warning('manopt:getGradient:approx', ...
    %        ['No gradient provided. Using an FD approximation instead (slow).\n' ...
    %         'It may be necessary to increase options.tolgradnorm.\n' ...
    %         'To disable this warning: warning(''off'', ''manopt:getGradient:approx'')']);
    %     problem0.approxgrad = approxgradientFD(problem0);
    % end
        
    condet = constraintsdetail(problem0);
    M = problem0.M;
    
    % Set localdefaults, a struct to be combined with argument options for
    % declaring hyperparameters.
    
    % Stopping criteria
    localdefaults.maxiter = 300;
    localdefaults.maxtime = 3600;
    localdefaults.tolKKTres = 1e-8;
    % StoreDB (TODO: implement)
    %   localdefaults.storedepth = 3;
    % Regularization of the hessian matrix to be positive-definete.
    localdefaults.modify_hessian = 'mineigval_matlab';
    localdefaults.mineigval_correction = 1e-8; 
    localdefaults.mineigval_threshold = 1e-3;
    % Initial parameters for the merit function and the Lagrangian
    localdefaults.tau = 0.5;  % as long as tau > 0
    localdefaults.rho = 1;  % as long as rho > 0
    localdefaults.beta = 0.9;  % as long as 1 > beta > 0
    localdefaults.gamma = 0.25; % as long as 1 > gamma > 0  
    localdefaults.mus = ones(condet.n_ineq_constraint_cost, 1);
    localdefaults.lambdas = ones(condet.n_eq_constraint_cost, 1);    
    % Linesearch
    localdefaults.ls_max_steps  = 10000;
    localdefaults.ls_threshold = 1e-8;
    % Display
    localdefaults.verbosity = 3;
    localdefaults.qp_verbosity = 0;
    % Only when using fixed-rank manifolds
    localdefaults.rankviopena = 1e+8;
    
    % Merge global and local defaults, then merge w/ user options, if any.
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    % Set the quadprog verbosities
    if options.qp_verbosity == 0
        qpoptions = optimset('Display','off');
    else
        qpoptions = [];
    end
    
    % Create a random starting point if no starting point is provided.
    if ~exist('x0', 'var')|| isempty(x0)
        xCur = problem0.M.rand(); 
    else
        xCur = x0;
    end
    
    % Only for calculating the KKT residual on fixed-rank manifolds 
    if contains(problem0.M.name(),'rank')
        if isfield(options, 'rank')
            rankval = options.rank;
        else
            tmpx = problem0.M.rand();
            rankval = rank(tmpx.S);
        end
    end
    
    % Create a store database and get a key for the current x (TODO)
    % storedb = StoreDB(options.storedepth);
    % key = storedb.getNewKey();
    
    % Create some initial variables which will be used in the following
    % loop.
    mus = options.mus;  % init. mus and lambdas for the Lagrangian
    lambdas = options.lambdas;
    rho = options.rho;  % init. rho for merit function
    
    % Declare some variables for the initial savestats
    iter = 0;
    xCurCost = getCost(problem0, xCur);
    xCurLagGrad = gradLagrangian(xCur, mus, lambdas);
    xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);
    xCurResidual = KKT_residual(xCur, mus, lambdas);
    [xCurMaxViolation, xCurMeanViolation] = const_evaluation(xCur);
    timetic = tic();
    
    % Save stats in a struct array info, and preallocate.
    stats = savestats();
    info(1) = stats;
    info(min(10000, options.maxiter+1)).iter = [];
    
    % Stopping flag, finally it should be true
    stop = false;
    totaltime = tic();
    
    % Main loop where we solve subproblems iteratively
    while true
        if options.verbosity >= 2
            fprintf('Iter: %d, Cost: %f, KKT residual: %.16e \n', iter, xCurCost, xCurResidual);
        elseif options.verbosity >= 1
            if mod(iter, 100) == 0 && iter ~= 0
                fprintf('Iter: %d, Cost: %f, KKT resiidual: %.16e \n', iter, xCurCost, xCurResidual);
            end
        end
        
        iter = iter + 1;
        timetic = tic();

        % Flags about updating parameters for a stopping criterion
        updateflag_rho = false;
        updateflag_Lagmult = false;
        
        % Get current Hessian and gradient of the cost function.
        % Also, make a "qpinfo" structure stading for the subproblem
        % at the current point.
        costLag = @(X) costLagrangian(X, mus, lambdas);
        gradLag = @(X) gradLagrangian(X, mus, lambdas); % in the tangent space
        hessLag = @(X, d) hessLagrangian(X, d, mus, lambdas); % in the tangent space
        auxproblem.M = problem0.M;
        auxproblem.cost = costLag;
        auxproblem.grad = gradLag;
        auxproblem.hess = hessLag;
        qpinfo = struct();

        %% Make H, basis, and n and modify H to be positive-definite in a 
        % predescribed manner when necessary.
        % (The difference between mineiegval_matlab and mineigval_manopt is
        % correct negative eigenvalues respectively or altogther.)
        if strcmp(options.modify_hessian, "eye")
            qpinfo.basis = tangentorthobasis(auxproblem.M, xCur, auxproblem.M.dim());
            qpinfo.n = numel(qpinfo.basis);
            qpinfo.H = eye(qpinfo.n);  % the identity matrix as replacement to Hessian.
        elseif strcmp(options.modify_hessian, 'mineigval_matlab') 
            [qpinfo.H, qpinfo.basis] = hessianmatrix(auxproblem, xCur);
            qpinfo.n = numel(qpinfo.basis);
            [U,T] = schur(qpinfo.H);
            for i = 1 : qpinfo.n
                if T(i,i) < 1e-5  
                    T(i,i) = options.mineigval_correction;
                end
            end
            qpinfo.H = U * T * U';
        elseif strcmp(options.modify_hessian, 'mineigval_manopt')
            [qpinfo.H, qpinfo.basis] = hessianmatrix(auxproblem, xCur);
            qpinfo.n = numel(qpinfo.basis);
            [~ ,qpinfo.mineigval] = hessianextreme(auxproblem, xCur);
            if qpinfo.mineigval < 0
                qpinfo.mineigval_diagcoeff = max(options.mineigval_threshold,...
                    abs(qpinfo.mineigval)) + options.mineigval_correction;
                qpinfo.H = qpinfo.H + qpinfo.mineigval_diagcoeff * eye(qpinfo.n);
            end
        else
            [qpinfo.H,qpinfo.basis] = hessianmatrix(auxproblem, xCur);  % may not be positive-definite
            qpinfo.n = numel(qpinfo.basis);
        end
        qpinfo.H = 0.5 * (qpinfo.H.'+qpinfo.H);
        
        % Make f
        f = zeros(qpinfo.n, 1);
        xCurGrad = getGradient(problem0, xCur);
        for fidx =1:qpinfo.n
            f(fidx) = problem0.M.inner(xCur, xCurGrad, qpinfo.basis{fidx});
        end
        qpinfo.f = f;

        % Make inequality constraints 
        if condet.has_ineq_cost
            row = condet.n_ineq_constraint_cost;
            col = qpinfo.n;
            A = zeros(row, col);
            b = zeros(row, 1);
            for ineqrow = 1:row
                ineqcosthandle = problem0.ineq_constraint_cost{ineqrow};
                b(ineqrow) = - ineqcosthandle(xCur);

                ineqgradhandle = problem0.ineq_constraint_grad{ineqrow};
                ineqconstraint_egrad = ineqgradhandle(xCur);
                ineqconstraint_grad = problem0.M.egrad2rgrad(xCur, ineqconstraint_egrad);

                for ineqcol = 1:col
                    base = qpinfo.basis{ineqcol};
                    A(ineqrow,ineqcol) = problem0.M.inner(xCur, ineqconstraint_grad, base);
                end
            end
        else
            A = [];
            b = [];
        end
        qpinfo.A = A;
        qpinfo.b = b;

        % Make equality constraints
        if condet.has_eq_cost
            row = condet.n_eq_constraint_cost;
            col = qpinfo.n;
            Aeq = zeros(row, col);
            beq = zeros(row, 1);
            for eqrow = 1:row
                eqcosthandle = problem0.eq_constraint_cost{eqrow};
                beq(eqrow) = - eqcosthandle(xCur);

                eqgradhandle = problem0.eq_constraint_grad{eqrow};
                eqconstraint_egrad = eqgradhandle(xCur);
                eqconstraint_grad = problem0.M.egrad2rgrad(xCur, eqconstraint_egrad);

                for eqcol = 1:col
                    base = qpinfo.basis{eqcol};
                    Aeq(eqrow,eqcol) = problem0.M.inner(xCur, eqconstraint_grad, base);
                end
            end
        else
            Aeq = [];
            beq = [];
        end
        qpinfo.Aeq = Aeq;
        qpinfo.beq = beq;

        %% Compute the direction and Lagrange multipliers by solving QP with
        % quadprog, a matlab solver for QP.
        [coeff, ~, qpexitflag, ~, Lagmultipliers] = quadprog(qpinfo.H, qpinfo.f,...
                qpinfo.A, qpinfo.b, qpinfo.Aeq, qpinfo.beq, [], [], [], qpoptions);     
        
        deltaXast = problem0.M.zerovec(xCur);
        for i = 1:qpinfo.n
            deltaXast = problem0.M.lincomb(xCur, 1, deltaXast, coeff(i), qpinfo.basis{i});
        end
        
        % Update rho, a penalty parameter, if necessary.
        newacc = rho;
        if condet.has_ineq_cost
            for iterineq = 1 : condet.n_ineq_constraint_cost
                newacc = max(newacc, Lagmultipliers.ineqlin(iterineq));
            end
        end
        if condet.has_eq_cost
            for itereq = 1 : condet.n_eq_constraint_cost
                newacc = max(newacc, abs(Lagmultipliers.eqlin(itereq)));
            end
        end
        if rho < newacc
           rho = newacc + options.tau;
           updateflag_rho = true;
        end
        
        %% Comnstruct a problem structure and some variables for the line
        % search.
        meritproblem.M = problem0.M;
        meritproblem.cost = @(x) loneMeritFunction(x, rho);
        f0 = meritproblem.cost(xCur);
        
        % Compute df0 according to options.modify_hessian
        df0 = (coeff.') * (qpinfo.H) * (coeff);

        % Compute the stepsize with the ell-1 type merit function and the Armijo rule
        stepsize = 1;
        newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
        newf = meritproblem.cost(newx);
        gammadf0 = df0 * options.gamma;
        r = 0; % backtracking counter
        ls_max_steps_flag = false;
        
        % DEBUG only
        % descriptCost(meritproblem, xCur, deltaXast);
        
        while newf > ( f0 - gammadf0) && abs(newf - ( f0 - gammadf0)) > options.ls_threshold
            if r > options.ls_max_steps
                ls_max_steps_flag = true;
                break;
            end
            r = r + 1;
            stepsize = stepsize * options.beta;
            gammadf0 =  gammadf0 * options.beta;
            newx = meritproblem.M.retr(xCur, deltaXast, stepsize);
            newf = meritproblem.cost(newx);
        end
        
        %% Information for savestats (TODO: the procedure should be out of
        % the solver because the distance calculation depends on manifolds.)
        if contains(problem0.M.name(),'Stiefel') 
            dist = norm(xCur - newx, 'fro');
        elseif contains(problem0.M.name(),'rank')
            % Only assuming for 'fixedrankembeddedfactory'
            if ~exist('xCurmat', 'var')
                xCurmat = xCur.U * xCur.S * xCur.V';
            end    
            newxmat = newx.U * newx.S * newx.V';
            dist = norm(xCurmat - newxmat, 'fro');
            xCurmat = newxmat;
        else
            dist = problem0.M.dist(xCur, newx);
        end
        
        % Update variables to new iterate
        xCur = newx;
        
        if ~(isequal(mus, Lagmultipliers.ineqlin)) ...
                || ~(isequal(lambdas, Lagmultipliers.eqlin))
            updateflag_Lagmult = true;
        end
        mus = Lagmultipliers.ineqlin;
        lambdas =  Lagmultipliers.eqlin;
        
        % Information for savestats (Cont'd)
        xCurCost = getCost(problem0, xCur);
        xCurLagGrad = gradLagrangian(xCur, mus, lambdas);
        xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);        
        xCurResidual = KKT_residual(xCur, mus, lambdas);
        [xCurMaxViolation, xCurMeanViolation] = const_evaluation(xCur);
        
        % Savestats
        %key = storedb.getNewKey();
        stats = savestats();
        info(iter+1) = stats;
        
        % Refer to stopping criteria        
        if iter >= options.maxiter
            fprintf('Max iter count reached\n');
            options.reason = "Max iter count reached";
            stop = true;
        elseif toc(totaltime) >= options.maxtime
            fprintf('Max time exceeded\n');
            options.reason = "Max time exceeded";
            stop = true;
        elseif xCurResidual <= options.tolKKTres
            fprintf('KKT Residual tolerance reached\n');
            options.reason = "KKT Residual tolerance reached";
            stop = true;
        elseif dist == 0 && ~(updateflag_rho) && ~(updateflag_Lagmult)
            % which will never occur because of the Armijo rule
            fprintf('Any parameter was not updated');
            options.reason = 'Any parameter was not updated';
            stop = true;
        end
        
        if stop
            options.totaltime = toc(totaltime);
            break
        end
        
    end
    
    xfinal = xCur;
    
    residual  = xCurResidual;

    costfinal = problem0.cost(xfinal);
    
    info = info(1:iter+1); % Crop blank of info

    if options.verbosity >= 1
        fprintf('Total time is %f [s] (excludes statsfun)\n', info(end).time);
    end

    %% Routine in charge of collecting the current iteration stats
    function stats = savestats()
        stats.xcurrent = xCur; % DEBUG only
        stats.iter = iter;
        stats.cost = xCurCost;
        stats.gradnorm = xCurLagGradNorm;
        if iter == 0
            stats.time = toc(timetic);
            stats.stepsize = NaN;
            stats.ls_max_steps_break = NaN;
            stats.dist =  NaN;
            stats.qpexitflag = NaN;
        else
            stats.time = toc(timetic);
            stats.time = info(iter).time + toc(timetic);
            stats.stepsize = stepsize;
            stats.ls_max_steps_break = ls_max_steps_flag;
            stats.dist = dist;
            stats.qpexitflag = qpexitflag;
        end
        stats.rho = rho;
        stats.KKT_residual = xCurResidual;
        stats.maxviolation = xCurMaxViolation;
        stats.meanviolation = xCurMeanViolation;
        stats = applyStatsfun(problem0, xCur, [], [], options, stats);
    end
    %% cost Lagrangian
    function val = costLagrangian(x, mus, lambdas)
        val = getCost(problem0, x);

        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + mus(numineq) * cost_numineq;
            end
        end

        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + lambdas(numeq) * cost_numeq;
            end
        end
    end
    %% grad Lagrangian
    function gradLag = gradLagrangian(x, mus, lambdas)
        gradLag = getGradient(problem0, x); % getGradient!!

        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                gradhandle = problem0.ineq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);

                gradLag = problem0.M.lincomb(x, 1, gradLag, mus(numineq), constraint_grad);
            end
        end

        if condet.has_eq_cost
            for numeq = 1:condet.n_eq_constraint_cost
                gradhandle = problem0.eq_constraint_grad{numeq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem0.M.egrad2rgrad(x, constraint_grad);

                gradLag = problem0.M.lincomb(x, 1, gradLag, lambdas(numeq), constraint_grad);
            end
        end
    end
    %% hess Lagrangian
    function hessLag = hessLagrangian(x, dir, mus, lambdas)
        hessLag = getHessian(problem0, x, dir); % getGradient!!

        if condet.has_ineq_cost
            for numineq = 1 : condet.n_ineq_constraint_cost
                gradhandle = problem0.ineq_constraint_grad{numineq};
                constraint_egrad = gradhandle(x);

                hesshandle = problem0.ineq_constraint_hess{numineq};
                constraint_ehess = hesshandle(x, dir);

                constraint_hess = problem0.M.ehess2rhess(x, constraint_egrad,...
                                                         constraint_ehess, dir);
                hessLag = problem0.M.lincomb(x, 1, hessLag,...
                    mus(numineq), constraint_hess);
            end
        end
        if condet.has_eq_cost
            for numeq = 1 : condet.n_eq_constraint_cost
                gradhandle = problem0.eq_constraint_grad{numeq};
                constraint_egrad = gradhandle(x);

                hesshandle = problem0.eq_constraint_hess{numeq};
                constraint_ehess = hesshandle(x, dir);

                constraint_hess = problem0.M.ehess2rhess(x, constraint_egrad,...
                                                        constraint_ehess, dir);

                hessLag = problem0.M.lincomb(x, 1, hessLag,...
                    lambdas(numeq), constraint_hess);
            end
        end
    end
    %% merit function (3.2) in paper
    function val = loneMeritFunction(x, rho)
        val = getCost(problem0, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + rho * max(0, cost_numineq);
            end
        end

        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + rho * abs(cost_numeq);
            end
        end
    end
    %%
    function [maxviolation, meanviolation] = const_evaluation(xCur)
        maxviolation = 0;
        meanviolation = 0;
        
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(xCur);
                maxviolation = max(maxviolation, cost_at_x);
                meanviolation = meanviolation + max(0, cost_at_x);
            end
        end
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_at_x = abs(costhandle(xCur));
                maxviolation = max(maxviolation, cost_at_x);
                meanviolation = meanviolation + cost_at_x;
            end
        end
        if condet.has_ineq_cost || condet.has_eq_cost
            meanviolation = meanviolation / (condet.n_ineq_constraint_cost + condet.n_eq_constraint_cost);
        end
    end
    
    %% Calculating the KKT residual at (XCur, mus, lambdas)
    function val = KKT_residual(xCur, mus, lambdas)
        xGrad = gradLagrangian(xCur, mus, lambdas);
        val = problem0.M.norm(xCur, xGrad)^2;
             
        compowvio = complementaryPowerViolation(xCur, mus);
        muspowvio = musposiPowerViolation(mus);
                
        val = val + compowvio;
        val = val + muspowvio;
        
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(xCur);
                violation = max(0, cost_at_x);
                val = val + violation^2;
            end
        end
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_at_x = abs(costhandle(xCur));
                val = val + cost_at_x^2;
            end
        end        

        val = sqrt(val);

        manpowvio = manifoldViolation(xCur);
        val = val + manpowvio;
                     
    end
    %%
    function compowvio = complementaryPowerViolation(xCur, mus)
        compowvio = 0;
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(xCur);
                violation = mus(numineq) * cost_at_x;
                compowvio = compowvio + violation^2;
            end
        end
    end
    %% mu>=0
    function musvio = musposiPowerViolation(mus)
        musvio = 0;
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                violation = max(-mus(numineq), 0);
                musvio = musvio + violation^2;
            end
        end
    end
     %% #Computing manifold violation
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