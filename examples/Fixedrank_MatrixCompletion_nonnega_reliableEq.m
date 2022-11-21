function success = Fixedrank_MatrixCompletion_nonnega_reliableEq

% solver RIPM
% name: matrix completion on nonnegative fixed-rank matrices space with
% some reliable sampled data application: nonnegative low-rank matrix
% completion with some reliable sampled data

% MITSUAKI OBARA

%% data setting

m=4; %%%%%%%%%%%%%%%%%%% q
n=2*m; %s
r=2; % rank p

L = rand(m, r);
R = rand(r, n);
A =  L*R;

% Generate a random mask for observed entries: P(i, j) = 1 if the entry
% (i, j) of A is observed, and 0 otherwise.
maskratio = 0.2; %%%%%%%%%%%%%%%%%%%
initP = zeros(m, n);
observed_indices = randperm(m*n, ceil(maskratio*m*n));
for obs = observed_indices
    initP(obs) = 1;
end

P=initP;
PA = P.*A;

%%% For nonnegativity inequality and equality
eqcont_ratio = 0.5; %%%%%%%%%%%%%%%%%%%
nonzero_num = nnz(P); % nz = NNZ(S) is the number of nonzero elements in S.
eqnum = ceil(nonzero_num*eqcont_ratio);
nonzeroidcs = find(P);

eqindices = sort(randsample(nonzeroidcs,eqnum));
P_eqcont=zeros(m, n);
for eqindex = eqindices'
    P_eqcont(eqindex) = 1;
end

%
P_eqcontA= P_eqcont.*A;

%% 
M = fixedrankembeddedfactory(m, n, r);

problem.M =M;
problem.cost=@(X) .5*norm(P.*M.triplet2matrix(X)-PA,"fro")^2;
problem.egrad=@(X) P.*M.triplet2matrix(X)-PA;
problem.ehess=@(X,dx) P.*M.triplet2matrix(M.tangent2ambient(X, dx));

problem.Euc_ineq = euclideanfactory(m, n);
problem.ineq_con=@(x) -M.triplet2matrix(x);
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -M.triplet2matrix(M.tangent2ambient(x, dx));%in fact, - indentity
ZERO = zeros(m,n);
problem.ehess_barGx = @(x, z, dx) ZERO;

problem.Euc_eq = euclideansparsefactory(P_eqcont);
problem.eq_con=@(X) P_eqcont.*M.triplet2matrix(X)-P_eqcontA;
problem.barHx=@(X,Y) P_eqcont.*Y; 
problem.barHxaj=@(X,dx) P_eqcont.*M.triplet2matrix(M.tangent2ambient(X, dx));
ZERO = zeros(m,n);
problem.ehess_barHx = @(x, z, dx) ZERO;

% RIPM_checkupto2ndorder(problem);

%% initial point X on manifold M

L = rand(m, r);
R = rand(r, n);
X_matrix = L*R;

x0 = M.matrix2triplet(X_matrix);

%%

%% 效果似乎不好
% options.z=rand(m,n);
% options.s=rand(m,n);

options.KrylovIterMethod=0;
options.tolKKTres=1e-15;

options.heuristic_z_s=1; % 很难收束如果不用 heuristic_z_s
options.important=0.13;
options.desired_tau_1=0.1;

options.maxtime = 15;
options.maxiter = 200;
options.verbosity = 3;

%RIPM_checkGradientHessian_ineq(problem);
[xfinal, costfinal, residual, info, options] = RIPM(problem, x0, options);

%X_matrix = M.triplet2matrix(W.x);

if residual <= options.tolKKTres
    success=1;
else
    success =0;
end

%% Display some statistics.

figure;
subplot(1,2,1)
semilogy([info.xCurPhi], '.-');
xlabel('Iteration number');
ylabel('phi');

subplot(1,2,2)
semilogy([info.KKT_residual], '.-');
xlabel('Iteration number');
ylabel('KKT residual');
end

