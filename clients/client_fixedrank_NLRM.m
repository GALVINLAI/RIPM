function data = client_fixedrank_NLRM(m, n, r, A, options, specifier, setting)

M = fixedrankembeddedfactory(m, n, r); % m: # of rows, n: # of cols, r: rank
problem.M = M;

%% Generating a random initial point x0.
L = rand(m, r);
R = rand(r, n);
x0_matrix = L*R;
x0 = M.matrix2triplet(x0_matrix);

%% Set-up objective function for all methods.

problem.cost=@(x) .5*norm(M.triplet2matrix(x)-A,"fro")^2;
problem.egrad=@(x) M.triplet2matrix(x)-A;
problem.ehess=@(x,dx) M.triplet2matrix(M.tangent2ambient(x, dx)); %in fact, indentity map.

%% Set-up constraints for RALM, REPMs, RSQP.

ineqnum = m*n;
nn_constraints_cost = cell(ineqnum, 1);
nn_constraints_grad = cell(ineqnum, 1);
nn_constraints_hess = cell(ineqnum, 1);

nnconst_idx = 1;
for row = 1:m
    for col = 1:n
        % cost
        nn_constraints_cost{nnconst_idx} = @(Y) nncostfun(Y, row, col);
        % egrad
        constraintgrad = zeros(m, n);
        constraintgrad(row, col) = -1;
        nn_constraints_grad{nnconst_idx} = @(U) constraintgrad;
        % ehess
        constrainthess = zeros(m, n);
        nn_constraints_hess{nnconst_idx} = @(X, U) constrainthess;
        nnconst_idx = nnconst_idx + 1;
    end
end

    function val = nncostfun(Y, row, col)
        Vt = Y.V.';
        val = - Y.U(row,:) * Y.S * Vt(:,col);
    end

problem.ineq_constraint_cost = nn_constraints_cost;
problem.ineq_constraint_grad = nn_constraints_grad;
problem.ineq_constraint_hess = nn_constraints_hess;

%% Set-up constraints for RIPM.

problem.Euc_ineq = euclideanfactory(m, n);
problem.ineq_con=@(x) -M.triplet2matrix(x);
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -M.triplet2matrix(M.tangent2ambient(x, dx));%in fact, - indentity map.
ZERO = zeros(m,n);
problem.ehess_barGx = @(x, z, dx) ZERO;

%% Calculating by solvers.

method_set = cell(2,5);
method_set(1,:) = {'RALM','REPM(LQH)','REPM(LSE)','RSQP','RIPM'};
method_set(2,:) = {@almbddmultiplier, @exactpenaltyViaSmoothinglqh, @exactpenaltyViaSmoothinglse, ...
    @SQP, @RIPM};

data = NaN(4, 5);

norm_A = norm(A,'fro');

for i = 1:5
    if specifier.ind(i)

        [method_name, method_handle]=method_set{:,i};
        fprintf('\n');
        fprintf('********** Experiment: %s\n', setting.ExperimentName);
        fprintf('********** Position: sd %.3f\tRow %d\tColumn %d\tRank %d\ttolKKTres %d\tRepeat No. %d\n', ...
            setting.sd,setting.row_dim,setting.col_dim,setting.rank,setting.tolKKTres,setting.repeat);
        fprintf('********** Starting Method: %s\n', method_name);

        timetic = tic();
        [xfinal, ~ , residual, info, ~] = method_handle(problem, x0, options);
        time = toc(timetic);
    
        % DEBUG ONLY.
%         filename = sprintf('RC_NLRM_%s_%s.csv',method_name, setting.filepath);
%         info = rmfield(info,'xcurrent');
%         struct2csv(info, filename);

        iternum = length(info);

        % for NLRM
        xfinal_matrix = M.triplet2matrix(xfinal);
        NLRMrelres = norm(xfinal_matrix- A,'fro')/ norm_A;

        data(:,i) = [residual; time; iternum; NLRMrelres];
    end
end


end