function [xfinal, costfinal, residual, info, options] = almbddmultiplier(problem0, x0, options)
    % The following code is based on that on https://github.com/losangle/Optimization-on-manifolds-with-extra-constraints
    % Symbol # in comments is added on newly added parts or modifications.
    
    condet = constraintsdetail(problem0);
    
    % Outer Loop Setting
    localdefaults.rho = 1;
    localdefaults.lambdas = ones(condet.n_ineq_constraint_cost, 1);
    localdefaults.gammas = ones(condet.n_eq_constraint_cost, 1);
    localdefaults.bound = 20;
    localdefaults.tau = 0.8;
    localdefaults.thetarho = 0.3;
    localdefaults.maxOuterIter = 300;
    localdefaults.numOuterItertgn = 30;
    localdefaults.outerverbosity = 3;  % #verbosity for outer loops
    localdefaults.tolKKTres = 1e-8;  % #a stopping criterion
    localdefaults.minstepsize = 1e-8;  % #
    % Inner Loop Setting
    localdefaults.maxInnerIter = 200;
    localdefaults.startingtolgradnorm = 1e-3;
    localdefaults.endingtolgradnorm = 1e-6;
    localdefaults.rankviopena = 1e+8;  % #rank check in KKT_residual(), only for fixed-rank manifolds,
    
    localdefaults = mergeOptions(getGlobalDefaults(), localdefaults);
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    options = mergeOptions(localdefaults, options);
    
    tolgradnorm = options.startingtolgradnorm;
    thetatolgradnorm = nthroot(options.endingtolgradnorm/options.startingtolgradnorm, options.numOuterItertgn);
    
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
    
    lambdas = options.lambdas;
    gammas = options.gammas;
    rho = options.rho;
    oldacc = Inf;
    M = problem0.M;
    xCur = x0;
    xPrev = xCur;
    OuterIter = 0;
    
    % #For savestats
    gradLagfun = @(X) grad_alm(X, problem0, rho, lambdas, gammas);
    complviofun = @(x) complementaryPowerViolation(x, rho, lambdas);
    xCurLagGrad = gradLagfun(xCur);
    xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);
    xCurMaxLagMult = maxabsLagrangemultipliers(lambdas, gammas);
    xCurResidual = KKT_residual();
    %
    
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
        
        costfun = @(X) cost_alm(X, problem0, rho, lambdas, gammas);
        gradfun = @(X) grad_alm(X, problem0, rho, lambdas, gammas);
        problem.cost = costfun;
        problem.grad = gradfun;
        problem.M = M;
        inneroptions.tolgradnorm = tolgradnorm;
        inneroptions.verbosity = 0;
        inneroptions.maxiter = options.maxInnerIter;
        inneroptions.minstepsize = options.minstepsize;
        
         
        [xCur, cost, innerinfo, Oldinneroptions] = rlbfgs(problem, xCur, inneroptions);
        
        % #For KKT_residual()
        gradLagfun = gradfun;
        complviofun = @(x) complementaryPowerViolation(x, rho, lambdas);
        %
        
        % #Updating flag for stopping algorithm
        updateflag_rho = false;
        updateflag_Lagmult = false;
        %
        
        % Update Multipliers
        newacc = 0;
        for iterineq = 1 : condet.n_ineq_constraint_cost
            costhandler = problem0.ineq_constraint_cost{iterineq};
            cost_iter = costhandler(xCur);
            newacc = max(newacc, abs(max(-lambdas(iterineq)/rho, cost_iter)));
            
            % #Update check 1
            newlambda = min(options.bound, max(lambdas(iterineq) + rho * cost_iter, 0));
            if lambdas(iterineq) ~= newlambda
                lambdas(iterineq) = newlambda;
                updateflag_Lagmult = true;
            end
            %
        end
        
        for itereq = 1 : condet.n_eq_constraint_cost
            costhandler = problem0.eq_constraint_cost{itereq};
            cost_iter = costhandler(xCur);
            newacc = max(newacc, abs(cost_iter));
            
            % #Update check 2
            newgamma = min(options.bound, max(-options.bound, gammas(itereq) + rho * cost_iter));
            if gammas(itereq) ~= newgamma
                gammas(itereq) = newgamma;
                updateflag_Lagmult = true;
            end
            %
        end
        
        if OuterIter == 1 || newacc > options.tau * oldacc
            rho = rho/options.thetarho;
            updateflag_rho = true;
        end
        oldacc = newacc;
        oldtolgradnorm = tolgradnorm;
        tolgradnorm = max(options.endingtolgradnorm, tolgradnorm * thetatolgradnorm);
        
        
        % #Information for savestats
        xCurLagGrad = gradLagfun(xCur);
        xCurLagGradNorm = problem0.M.norm(xCur, xCurLagGrad);
        xCurMaxLagMult = maxabsLagrangemultipliers(lambdas, gammas);
        xCurResidual = KKT_residual();
        %
        
        % #Calculating the distance for savestats
        if contains(problem0.M.name(),'Stiefel') 
            dist = norm(xCur - xPrev, 'fro');
        elseif contains(problem0.M.name(),'rank')  % #only assuming for 'fixedrankembeddedfactory'
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
        
        % #Verbosity modified
        if options.outerverbosity >= 2
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        elseif options.outerverbosity == 1 && mod(OuterIter, 100) == 0 
            fprintf('KKT Residual: %.16e\n', xCurResidual)
        end
        %
        
        % #Stopping criteria
        if xCurResidual < options.tolKKTres && tolgradnorm <= options.endingtolgradnorm  % #stopping criterion based on KKT_residual()
            fprintf("KKT Residual tolerance reached\n")
            break;
        elseif dist == 0 && ~(updateflag_rho) && ~(updateflag_Lagmult) ...
                && (tolgradnorm == oldtolgradnorm)  % #nothing changed, meaning that the algo. keeps producing the same point hereafter.
            fprintf("Any parameter did not change\n")
            break; 
        end
        %
        
        % The original one remains here, just in case
        % if norm(xPrev-xCur, 'fro') < options.minstepsize && tolgradnorm <= options.endingtolgradnorm
        %     break;
        % end
        
        if toc(totaltime) > options.maxtime
            break;
        end
        
        xPrev = xCur;
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
            stats.dist = NaN; % #
        else
            stats.time = info(OuterIter).time + toc(timetic);
            stats.dist = dist;
        end
        stats.rho = rho;
        [maxviolation, meanviolation, costCur] = evaluation(problem0, x, condet);
        stats.maxviolation = maxviolation;
        stats.meanviolation = meanviolation;
        stats.cost = costCur;
        % #Information on Lagrange multipliers and KKT residual
        stats.maxabsLagMult = xCurMaxLagMult;
        stats.KKT_residual = xCurResidual;
        %
        stats.LagGradNorm = xCurLagGradNorm;
    end
    %%
    function val = cost_alm(x, problem, rho, lambdas, gammas)
        val = getCost(problem, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                val = val + (rho/2) * (max(0, lambdas(numineq)/rho + cost_numineq)^2);
            end
        end
        
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                val = val + (rho/2) * (gammas(numeq)/rho + cost_numeq)^2;
            end
        end
    end
    %%
    function val = grad_alm(x, problem, rho, lambdas, gammas)
        val = getGradient(problem, x);
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(x);
                if (cost_numineq + lambdas(numineq)/rho > 0)
                    gradhandle = problem.ineq_constraint_grad{numineq};
                    constraint_grad = gradhandle(x);
                    constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                    val = problem.M.lincomb(x, 1, val, cost_numineq * rho + lambdas(numineq), constraint_grad);
                end
            end
        end
        
        if condet.has_eq_cost
            for numeq = 1:condet.n_eq_constraint_cost
                costhandle = problem.eq_constraint_cost{numeq};
                cost_numeq = costhandle(x);
                gradhandle = problem.eq_constraint_grad{numeq};
                constraint_grad = gradhandle(x);
                constraint_grad = problem.M.egrad2rgrad(x, constraint_grad);
                val = problem.M.lincomb(x, 1, val, cost_numeq * rho + gammas(numeq), constraint_grad);
            end
        end
    end

    %% #Computing KKT residual at the current iterate
     function val = KKT_residual()
        xGrad = gradLagfun(xCur);
        val = (problem0.M.norm(xCur, xGrad))^2;
       
        compowvio = complviofun(xCur);
        %lampowvio = lamposiPowerViolation(lambdas);  % #unnecessary
            
        val = val + compowvio;
        %val = val + lampowvio; % #unncecessary
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_at_x = costhandle(xCur);
                violation = max(0, cost_at_x);
                val = val + (violation)^2;
            end
        end
        if condet.has_eq_cost
            for numeq = 1: condet.n_eq_constraint_cost
                costhandle = problem0.eq_constraint_cost{numeq};
                cost_at_x = abs(costhandle(xCur));
                val = val + (cost_at_x)^2;                
            end
        end
        val = sqrt(val);

        manpowvio = manifoldViolation(xCur);
        val = val + manpowvio;
        
     end
 
    %%
    function compowvio = complementaryPowerViolation(xCur, rho, lambdas)
        compowvio = 0;
        if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                costhandle = problem0.ineq_constraint_cost{numineq};
                cost_numineq = costhandle(xCur);
                if (rho * cost_numineq + lambdas(numineq) > 0)
                    violation = (cost_numineq * rho + lambdas(numineq)) * cost_numineq;
                    compowvio = compowvio + violation^2;
                end
            end
        end
    end
    
    % #Unnecessary because we take max(0, \lambda^{k-1)_{i} +
    % \rho_{k-1}g_{i}(x_{k})) as lambdas, which are always nonnegative.
    %function musvio = lamposiPowerViolation(lambdas)
    %    musvio = 0;
    %end 
 
    % #Computing the maximum among the absolute values of Lagrange multiplies.
    function val = maxabsLagrangemultipliers(lambdas, gammas)
       val = -1; % #meaning no constraints
       if condet.has_ineq_cost
            for numineq = 1: condet.n_ineq_constraint_cost
                val = max(val, abs(lambdas(numineq)));
            end
        end
        
        if condet.has_eq_cost
            for numeq = 1:condet.n_eq_constraint_cost
                val = max(val, abs(gammas(numeq)));
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