function Fixedrank_MatrixCompletion_nonnega() 

% solver RIPM_ineq
% name: matrix completion on nonnegative fixed-rank matrices space
% application: nonnegative low-rank matrix completion

% f=.5*norm(P.*M.triplet2matrix(X)-PA,"fro")^2;
% M=fixed rank manifold
% g=nonngegative


%% data setting

m=4; %q
n=2*m; %s
r=2; % rank p

L = rand(m, r);
R = rand(r, n);
A =  L*R;

% Generate a random mask for observed entries: P(i, j) = 1 if the entry
% (i, j) of A is observed, and 0 otherwise.
maskratio = 0.1;
initP = zeros(m, n);
observed_indices = randperm(m*n, ceil(maskratio*m*n));
for obs = observed_indices
    initP(obs) = 1;
end

P=initP;
PA = P.*A;

%%
problem.A = A; % for computing relres
problem.Anorm=norm(problem.A,"fro"); % for computing relres


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

%% initial point X on manifold M

L = rand(m, r);
R = rand(r, n);
X_matrix = L*R;

x0 = M.matrix2triplet(X_matrix);


%%

options.checkNTequation=0;
options.KrylovIterMethod=1;
options.KrylovTolrelres=1e-10;
options.heuristic_z_s=1;
options.important = 0.1;
options.maxiter = 200;
%RIPM_checkGradientHessian_ineq(problem);
[xfinal, costfinal, phifinal, info, options] = RIPM(problem, x0, options);


%% Display

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
