
clc; clear; close all;
%% Sep. 18, 2023 Zhijian

% tuple.x1 = euclideanfactory(n), with h(x)=x'*x=1

%% data setting

n=20;

a1 = [-1;abs(rand(n-1,1))];
% a1 = rand(n,1);
a2 = rand(n,1);

%%
% tuple.x1 = spherefactory(n);% x1 
tuple.x1 = euclideanfactory(n);% x1 
tuple.x2 = euclideanfactory(n);% x2
problem.M = productmanifold(tuple);
%%
problem.cost=@(x) a1.'*x.x1 + a2.'*x.x2;
problem.egrad=@(x) struct('x1',a1,'x2',a2);
ZERO = zeros(n,1);
problem.ehess = @(x,dx) struct('x1', ZERO, 'x2', ZERO);
%% g(x)
Eineq_tuple.x1 = euclideanfactory(n);% x1 
Eineq_tuple.x2 = euclideanfactory(n);% x2
problem.Euc_ineq = productmanifold(Eineq_tuple);

problem.ineq_con=@(x) struct('x1',-x.x1,'x2',-x.x2+0.5);
problem.barGx=@(x,z) struct('x1',-z.x1,'x2',-z.x2);
problem.barGxaj=@(x,dx) struct('x1',-dx.x1,'x2',-dx.x2);
problem.ehess_barGx = @(x, z, dx) struct('x1',ZERO,'x2',ZERO);

%% h(x)
% Since there is only one eq_con, we don't need a product structure
% of codomain of h.
problem.Euc_eq = euclideanfactory(1,1);

problem.eq_con=@(x) (x.x1).'*(x.x1)-1;
problem.barHx=@(x, y) struct('x1',y*2*(x.x1),'x2', ZERO);
problem.barHxaj=@(x,dx) 2*(x.x1).'*(dx.x1);
problem.ehess_barHx = @(x, y, dx) struct('x1',y*2*(dx.x1),'x2',ZERO);

RIPM_checkupto2ndorder(problem)

%%

% initial point X on manifold M
Sp = spherefactory(n,1);
x0.x1=abs(Sp.rand());
x0.x2=ones(n,1);

options.checkNTequation=0;   
options.KrylovIterMethod=1;
options.important=0.5; %


[xfinal, costfinal, phifinal, info, options]= RIPM_ProductManifold(problem, x0, options);



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


% len=size(info, 2);
% dist_to_sol=zeros(len,1);
% for i=1:len
% dist_to_sol(i)=norm(info(i).xcurrent-problem.X_sol,'fro');
% end

% subplot(1,3,3)
% semilogy([dist_to_sol], '.-');
% xlabel('Iteration number');
% ylabel('dist to sol');



% end

