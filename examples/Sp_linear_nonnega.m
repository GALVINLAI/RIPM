function  Sp_linear_nonnega()

% solver: RIPM
% name: linear function on nonngegative sphere

% f=linear function
% M=sphere
% g=nonngegative

% guess solution is  [1;zeros(n-1,1)]


%% data setting
n=50;
a = [-1;abs(rand(n-1,1))];

% guess solution
problem.X_sol = [1;zeros(n-1,1)];

%% 

problem.M = spherefactory(n);
problem.cost=@(x) a.'*x;
problem.egrad=@(x) a;
ZERO = zeros(n,1);
problem.ehess = @(x,dx) ZERO;

problem.Euc_ineq = euclideanfactory(n);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
problem.ehess_barGx = @(x, z, dx) ZERO;

% initial point X on manifold M
x0=abs(problem.M.rand());

%options.maxtime=0.01;

options.checkNTequation=0;
options.KrylovIterMethod=0;
options.maxiter=50;

%RIPM_checkGradientHessian_ineq(problem);
[xfinal, costfinal, phifinal, info, options]= RIPM(problem, x0, options);

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
dist_to_sol=zeros(len,1);
for i=1:len
dist_to_sol(i)=norm(info(i).xcurrent-problem.X_sol,'fro');
end

subplot(1,3,3)
semilogy([dist_to_sol], '.-');
xlabel('Iteration number');
ylabel('dist to sol');


end
