function data = client_Stiefel_NonnegativeProjection(n, k, Xstar, C, options, specifier, setting)

% name: linear function on nonngegative orthogonality constraints
% f = -2*trace(X'*C)
% M = Stiefel(n,k)
% g = nonngegative
% ture solution is Xstar

%%

problem.M = stiefelfactory(n,k);

%% Generating an initial point x0.

x0 = Proj_Stiefel(C);

%% Set-up objective function for all methods.

problem.cost = @(X) -2*trace(X'*C);
problem.egrad = @(x) -2*C;
ZERO = zeros(n,k);
problem.ehess = @(x,dx) ZERO;

%% Set-up constraints for RALM, REPMs, RSQP.

ineqnum = n*k;
constraints_cost = cell(ineqnum, 1);
constraints_grad = cell(ineqnum, 1);
constraints_hess = cell(ineqnum, 1);

nnconst_idx = 1;
for row = 1: n
    for col = 1: k
        % cost
        constraints_cost{nnconst_idx} = @(Y) -Y(row, col);
        % egrad
        constraintgrad = zeros(n, k);
        constraintgrad(row, col) = -1;
        constraints_grad{nnconst_idx} = @(U) constraintgrad;
        % ehess
        constrainthess = zeros(n, k);
        constraints_hess{nnconst_idx} = @(X, U) constrainthess;
        nnconst_idx = nnconst_idx + 1;
    end
end

problem.ineq_constraint_cost = constraints_cost;
problem.ineq_constraint_grad = constraints_grad;
problem.ineq_constraint_hess = constraints_hess;

%% Set-up constraints for RIPM.

problem.Euc_ineq = euclideanfactory(n, k);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
ZERO = zeros(n,k);
problem.ehess_barGx = @(x, z, dx) ZERO;

%% Calculating by solvers.

method_set = cell(2,5);
method_set(1,:) = {'RALM','REPM(LQH)','REPM(LSE)','RSQP','RIPM'};
method_set(2,:) = {@almbddmultiplier, @exactpenaltyViaSmoothinglqh, @exactpenaltyViaSmoothinglse, ...
    @SQP, @RIPM};

data = NaN(4, 5);

for i = 1:5
    if specifier.ind(i)

        [method_name, method_handle]=method_set{:,i};
        fprintf('\n');
        fprintf('********** Experiment: %s\n', setting.ExperimentName);
        fprintf('********** Position: Row %d\tColumn %d\ttolKKTres %d\tRepeat No. %d\n', ...
            setting.row_dim,setting.col_dim,setting.tolKKTres,setting.repeat);
        fprintf('********** Starting Method: %s\n', method_name);

        timetic = tic();
        [xfinal, ~ , residual, info, ~] = method_handle(problem, x0, options);
        time = toc(timetic);

%         % DEBUG ONLY.
%         filename = sprintf('RC_NLRM_%s_%s.csv',method_name, setting.filepath);
%         info = rmfield(info,'xcurrent');
%         struct2csv(info, filename);

%         filename = sprintf('InfoDate_NLRM_%s_%s.mat',method_name, setting.filepath);
%         save(filename,"info","Xstar","setting","options");  % save data as *.mat files.

        iternum = length(info);

        % for
        DistToSolution = norm(xfinal- Xstar,'fro');

        data(:,i) = [residual; time; iternum; DistToSolution];
    end
end

end

