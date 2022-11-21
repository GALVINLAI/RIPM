function [xfinal, costfinal, phifinal, info, options]= Euc_linear_nonnega_sphereEq()
% solver: RIPM
% name: linear function on nonngegative sphere

% f=linear function
% M=Euc(n,1)
% h=sphere x'*x=1
% g=nonngegative

% guess solution is  [1;zeros(n-1,1)]


%% data setting

n=50;
a = [-1;abs(rand(n-1,1))];

% a guess solution
problem.X_sol = [1;zeros(n-1,1)];

%%
problem.M = euclideanfactory(n);
problem.cost=@(x) a.'*x;
problem.egrad=@(x) a;
ZERO = zeros(n,1);
problem.ehess = @(x,dx) ZERO;

problem.Euc_ineq = euclideanfactory(n);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
problem.ehess_barGx = @(x, z, dx) ZERO;

problem.Euc_eq = euclideanfactory(1,1);
problem.eq_con=@(x) x'*x-1;
problem.barHx=@(x,y) y*(2*x);
problem.barHxaj=@(x,dx) 2*(x.'*dx);
problem.ehess_barHx = @(x, y, dx) y*2*dx;


% initial point X on manifold M
x0=abs(problem.M.rand());

options.checkNTequation=0;
options.KrylovIterMethod=1;
options.heuristic_z_s=1;
%RIPM_checkGradientHessian(problem);
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

