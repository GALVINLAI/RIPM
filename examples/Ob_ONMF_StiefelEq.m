function  Ob_ONMF_StiefelEq()

% solver: RIPM
% name: Orthogonal nonnegative matrix factorization.

%% If no input is provided, we generate a quick demo.
if nargin == 0
    n=4;m=2*n;k=2;
end

%% data setting
% generate V
V=ones(k,1);
V=V/norm(V,"fro");
VVt=V*V';

%check V
%omega=min(VVt,[],'all');
%norm(V,"fro")

tilde_X=rand_Nonneg_St(n,k);% random B
C = rand(k,m); 
D = rand(n,m);
A = tilde_X*C + 0.01*D; % A: n x m.
AAt = A*A';

%%
problem.M=obliquefactory(n,k);

problem.cost=@(X) -trace(X'*AAt*X);
problem.egrad=@(X) -2*AAt*X;
problem.ehess = @(x,dx) -2*AAt*dx;

problem.Euc_ineq = euclideanfactory(n, k);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
ZERO = zeros(n,k);
problem.ehess_barGx = @(x, z, dx) ZERO;

problem.Euc_eq = euclideanfactory(1);
problem.eq_con=@(X) norm(X*V,'fro')^2-1;
problem.barHx=@(X,Y) (2*Y)*(X*VVt);
problem.barHxaj=@(X,dx) F_inner(2*X*VVt,dx);
problem.ehess_barHx=@(X,Y,dx) (2*Y)*(dx*VVt);

% RIPM_checkupto2ndorder(problem)

%% initial point X on manifold M

x0 = rand_Nonneg_St(n,k);

%

options.KrylovIterMethod=0;
options.heuristic_z_s=1;
options.important=0.5;
options.desired_tau_1=0.5;
options.maxiter=200;
options.tolKKTres=1e-14;

[xfinal, costfinal, phifinal, info, options]= RIPM(problem, x0, options);

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

%%
