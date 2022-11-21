function [xfinal, costfinal, residual, info, options] = exactpenaltyViaSmoothinglse (problem0, x0, options)
    % The following code is based on that on https://github.com/losangle/Optimization-on-manifolds-with-extra-constraints
    % Symbol # in comments is added on newly added parts or modifications.
    
    condet = constraintsdetail(problem0);
    
    % Outer Loop Setting
    localdefaults.rho = 1;
    localdefaults.thetarho = 0.3;
    localdefaults.maxOuterIter = 300;
    localdefaults.numOuterItertgn = 30;
    localdefaults.startingepsilon = 1e-1;
    localdefaults.endingepsilon = 1e-6;
    localdefaults.outerverbosity = 1;  % #verbosity for outer loops
    localdefaults.tolKKTres = 1e-8; % #a stopping criterion
    localdefaults.minstepsize = 1e-8;  % #
    % Inner Loop Setting
    localdefaults.maxInnerIter = 200;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-6;
    localdefaults.rankviopena = 1e+8;  % #rank check in KKT_residual(), only for fixed-rank manifolds
    
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    % #For the KKT_residual(), only when considering a fixed-rank manifold
    if contains(problem0.M.name(),'rank')
        if isfield(options, 'rank')
            rankval = options.rank;
        else
            tmpx = problem0.M.rand();
            rankval = rank(tmpx.S);
        end
    end
    %
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    theta_epsilon = nthroot(options.endingepsilon/options.startingepsilon, options.numOuterItertgn);
    
    M = problem0.M;
    xCur = x0;
    xPrev = xCur;
    epsilon = options.startingepsilon;
    rho = options.rho;
    
    
    % #For savestats
    xCurMaxLagMult = maxabsLagrangemultipliers(xCur, problem0, epsilon);
    gradfun = @(X) grad_exactpenalty(X, problem0, rho);
    xCurResidual = KKT_residual();
    %
    
    OuterIter = 0;
    stats = savestats(x0);
    info(1) = stats;
    info(min(10000, options.maxOuterIter+1)).iter = [];
    
    totaltime = tic();
    
    for OuterIter = 1 : options.maxOuterIter
        timetic = tic();
        
        % #Verbosity modified
        if options.outerverbosity >= 2
            fprintf('Iteration: %d    ', OuterIter);
        elseif options.outerverbosity == 1 && mod(OuterIter, 100) == 0 
            fprintf('Iteration: %d    ', OuterIter);
        end
        %
        
        costfun = @(X) cost_exactpenalty(X, problem0, rho);
        gradfun = @(X) grad_exactpenalty(X, problem0, rho);
        problem.cost = costfun;
        problem.grad = gradfun;
        problem.M = M;
        
        inneroptions.tolgradnorm = tolgradnorm;
        inneroptions.verbosity = 0;
        inneroptions.maxiter = options.maxInnerIter;
        inneroptions.minstepsize = options.minstepsize;
        
        [xCur, cost, innerInfo, Oldinneroptions] = rlbfgs(problem, xCur, inneroptions);
        
        % #For savestats
        xCurMaxLagMult = maxabsLagrangemultipliers(xCur, problem0, epsilon);
        xCurResidual = KKT_residual();
        
        % #Calculating the distance for savestats
        if contains(problem0.M.name(),'Stiefel') 
            dist = norm(xCur - xPrev, 'fro');
        elseif contains(problem0.M.name(),'rank')  % #only for 'fixedrankembeddedfactory'
            if ~exist('xCurmat', 'var')
                xCurmat = xCur.U * xCur.S * xCur.V';
            end    
            xPrevmat = xPrev.U * xPrev.S * xPrev.V';
            dist = norm(xCurmat - xPrevmat, 'fro');
            xCurmat = xPrevmat;
        else
            dist = problem0.M.dist(xCur, xPrev);
        end
        %
        
        % Save stats
        stats = savestats(xCur);
        info(OuterIter+1) = stats;
        
        % #According to Proposition 4.2 in [Liu and Boumal, 2019], we turn
        % off the following update.
        %if stats.maxviolation > epsilon
        %    rho = rho/options.thetarho;
        %end
        %
        
        % #Update check
        oldeps = epsilon;
        oldtolgradnorm = tolgradnorm;
        %
        
        epsilon  = max(options.endingepsilon, theta_epsilon * epsilon);
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm);
                
        % #Verbosity modified 
        if options.outerverbosity >= 2
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        elseif options.outerverbosity == 1 && mod(OuterIter, 100) == 0 
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        end
        %
        
        % #Stopping criteria
        if xCurResidual < options.tolKKTres && tolgradnorm <= options.endingtolgradnorm
            fprintf("KKT residual tolerance reached\n")            
            break;
        elseif dist == 0 && (epsilon == oldeps) && (tolgradnorm == oldtolgradnorm)  % #nothing changed, meaning that the alg. keeps producing the same point hereafter.          
            fprintf("Any parameter did not change\n")
            break;
        end
        %
        
        xPrev = xCur;
        
        if toc(totaltime) > options.maxtime
            break;
        end
    end

    info = info(1: OuterIter+1);

    residual  = KKT_residual(); % #

    xfinal = xCur;

    costfinal = getCost(problem, xCur); % #

    if options.verbosity >= 1
    fprintf('Total time is %f [s] (excludes statsfun)\n', info(end).time);
    end
    
    %%
    function stats = savestats(x) 
        stats.iter = OuterIter;
        stats.xcurrent = xCur;
        if stats.iter == 0
            stats.time = 0;
            stats.dist = NaN;  % #
        else
            stats.time = info(OuterIter).time + toc(timetic);
            stats.dist = dist;  % #
        end
        [maxviolation, meanviolation, costCur] = evaluation(problem0, x, condet);
        stats.maxviolation = maxviolation;
        stats.meanviolation = meanviolation;
        stats.cost = costCur;
        % #Information on Lagrange multipliers and KKT residual
        stats.maxabsLagMult = xCurMaxLagMult; 
        stats.KKT_residual = xCurResidual; 
        %
    end
    %%
    function val = cost_exactpenalty(x, problem, rho)
        val = getCost(problem, x);
        % Adding ineq constraint cost
        if condet.has_ineq_cost
            for numineq = 1 : condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(x);
                s = max(0, cost_at_x);
                additional_cost = s + epsilon * log( exp((cost_at_x - s)/epsilon) + exp(-s/epsilon));
                val = val + rho * additional_cost;
            end
        end
        % Eq constraint cost
        if condet.has_eq_cost
            for numeq = 1 : condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numeq};
                cost_at_x = costhandle(x);
                s = max(-cost_at_x, cost_at_x);
                additional_cost = s + epsilon * log( exp((cost_at_x - s)/epsilon) + exp((-cost_at_x-s)/epsilon));
                val = val + rho * additional_cost;
            end
        end
    end
    %%
    function val = grad_exactpenalty(x, problem, rho)
        val = getGradient(problem, x);
        if condet.has_ineq_cost
            for numineq = 1 : condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                s = max(0, cost_at_x);
                gradhandle = problem.ineq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                coef = rho * exp((cost_at_x-s)/epsilon)/(exp((cost_at_x-s)/epsilon)+exp(-s/epsilon));
                val = problem.M.lincomb(x, 1, val, coef, constraint_grad);
            end
        end
        if condet.has_eq_cost
            for numineq = 1 : condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                s = max(-cost_at_x, cost_at_x);
                gradhandle = problem.eq_constraint_grad{numineq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                coef = rho * (exp((cost_at_x-s)/epsilon)-exp((-cost_at_x-s)/epsilon))/(exp((cost_at_x-s)/epsilon)+exp((-cost_at_x-s)/epsilon));
                val = problem.M.lincomb(x, 1, val, coef, constraint_grad);
            end 
        end
    end
    %%
    % #Computing KKT residual at the current iterate
    function val = KKT_residual()
        grad = gradfun(xCur);
        val = problem0.M.norm(xCur, grad)^2;
        
        compowvio = complementaryPowerViolation(xCur, problem0, epsilon);
        val = val + compowvio;
        
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

    % #Calculating the violation on the complementary condition
    function compowvio = complementaryPowerViolation(x, problem, u)
        compowvio = 0;
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                s = max(0, cost_at_x);
                % Note the following code includes rho
                lambda = rho * exp((cost_at_x-s)/epsilon)/(exp((cost_at_x-s)/epsilon)+exp(-s/epsilon));
                violation = lambda * cost_at_x;
                compowvio = compowvio + violation^2;
            end
        end
    end

    % #Computing the maximum among the absolute values of Lagrange multiplies.
    function val = maxabsLagrangemultipliers(x, problem, u)
        val = -1; % #meaning no constraints
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                lambda = 1 / ( 1 + exp(- cost_at_x / u) );
                val = max(val, abs(lambda));
            end
        end
        if condet.has_eq_cost
            for numineq = 1: condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numineq};            
                cost_at_x = costhandle(x);
                exp_plus = exp(cost_at_x / u);
                exp_minus = exp(-cost_at_x / u);
                gamma = (exp_plus - exp_minus) / (exp_plus + exp_minus);
                val = max(val, abs(gamma));
            end 
        end
    end

    %% #Computing manifold violation
    function manvio = manifoldViolation(xCur)  % # No square!
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