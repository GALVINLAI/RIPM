function Fixedrank_projection_nonnega() 
% old name is client_NLRM.m

% solver RIPM_ineq
% name: projection onto nonnegative fixed-rank matrices space
% application: nonnegative low rank matrix approximation (NLRM)

% Notice that X is tripet form of svd in fixed-rank manifold.

% f = .5*norm(M.triplet2matrix(X)-A,"fro")^2;
% M = fixed rank manifold
% g = nonngegative

% ture solution is original A

%% data matrix A
m = 20; % row
n = 0.8*m; % column
r = max(2, 0.1*m); % rank

%% Experiment #1. Known low rank minimizer.
B = rand(m, r);
C = rand(r, n);
origanl_A = B*C; % ture solution is original A

% Gaussian Noise
sd = 0.001; % sd is standard deviation, for Gaussian Noise.
mu = 0; % mean
vr = sd.^2; % variance
GaussianNoise = mu + sqrt(vr)*randn(m,n);

A = origanl_A + GaussianNoise; % data matrix A becomes full rank.

%% Experiment #2. Unknown low rank minimizer.
% A = rand(m, n);

%%
problem.A = A; % for computing relres
problem.Anorm=norm(problem.A,"fro"); % for computing relres

%% initial point x0 on manifold M
L = rand(m, r);
R = rand(n, r);
X0_matrix = L*R';
M = fixedrankembeddedfactory(m, n, r);
x0 = M.matrix2triplet(X0_matrix);

%%
M = fixedrankembeddedfactory(m, n, r);
problem.M =M;
problem.cost=@(x) .5*norm(M.triplet2matrix(x)-A,"fro")^2;
problem.egrad=@(x) M.triplet2matrix(x)-A;
problem.ehess=@(x,dx) M.triplet2matrix(M.tangent2ambient(x, dx)); %in fact, indentity

problem.Euc_ineq = euclideanfactory(m, n);
problem.ineq_con=@(x) -M.triplet2matrix(x);
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -M.triplet2matrix(M.tangent2ambient(x, dx));%in fact, - indentity
ZERO = zeros(m,n);
problem.ehess_barGx = @(x, z, dx) ZERO;

%%

options.checkNTequation=0;
options.KrylovIterMethod=1;

[xfinal, costfinal, phifinal, info, options] = RIPM(problem, x0, options);


%% Display
figure;
subplot(1,3,1)
semilogy([info.xCurPhi], '.-');
xlabel('Iteration number');
ylabel('phi');

subplot(1,3,2)
semilogy([info.KKT_residual], '.-');
xlabel('Iteration number');
ylabel('KKT residual');

len=size(info, 2);
relative_residuals=zeros(len,1);
for i=1:len
    X_matrix= problem.M.triplet2matrix(info(i).xcurrent);
    relative_residuals(i)=norm(X_matrix-problem.A,'fro')/problem.Anorm;
end

subplot(1,3,3)
semilogy([relative_residuals], '.-');
xlabel('Iteration number');
ylabel('relative residuals');


end



