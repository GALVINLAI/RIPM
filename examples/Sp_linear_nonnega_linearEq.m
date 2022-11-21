function Sp_linear_nonnega_linearEq()

% sovler: RIPM
% name: the first example of RIPM (dim = 3)

% f=linear function
% M=sphere
% h=linear function
% g=nonngegative

% only considered for dim = 3.

% true solution is [sqrt(2)/2;sqrt(2)/2;0];

%% data setting

n = 3;
a = [-1 2 1]';
b=[-1 1 -1]';

% true solution
problem.X_sol = [sqrt(2)/2;sqrt(2)/2;0];

%% 
problem.M = spherefactory(n);
problem.cost = @(x) a.'*x;
problem.egrad = @(x) a;
ZERO = zeros(n,1);
problem.ehess = @(x,dx) ZERO;

problem.Euc_ineq = euclideanfactory(n,1);
problem.ineq_con = @(x) -x;
problem.barGx = @(x,z) -z;
problem.barGxaj = @(x,dx) -dx;
problem.ehess_barGx = @(x, z, dx) ZERO;

problem.Euc_eq = euclideanfactory(1,1);
problem.eq_con=@(x) b.'*x;
problem.barHx=@(x,y) y*b;
problem.barHxaj=@(x,dx) b.'*dx;
problem.ehess_barHx = @(x, y, dx) ZERO; % ehess_{x} <y,h(x)>

% initial point X on manifold M
x0=abs(problem.M.rand());

%options.tolphi=1e-5;
%options.maxtime=0.01;
%options.maxiter=5;

options.checkNTequation=1;
options.KrylovIterMethod=0; 
options.heuristic_z_s=1;
options.maxiter = 200;
%RIPM_checkGradientHessian(problem);
[xfinal, costfinal, residual, info, options] = RIPM(problem, x0, options);


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
