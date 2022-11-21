function success = Ob_linear_or_projection_nonnega_StiefelEq(n,k)

% solver: RIPM
% name: linear/projection function on nonngegative orthogonality constraints


% f = -2*trace(X'*C) % f is linear function, if linear_modle=ture
% or, f = norm(X-C,"fro").^2  % f is projection function, if linear_modle=false

% M = oblique(n,k)
% h = norm(X*V)^2-1  % StiefelEq
% g = nonngegative

% ture solution is Xstar


%% If no input is provided, we generate a quick demo.
if nargin == 0
    n=50;k=0.2*n;
end

%% data setting
% generate V
V=ones(k,1);
V=V/norm(V,"fro");
VVt=V*V';

%check V
%omega=min(VVt,[],'all');
%norm(V,"fro")

% construct data C
% step 1. generate a Nonnegative_Stiefel Xstar.
B=rand_Nonneg_St(n,k);% random B
X1= (B>0).*(1+rand(n,k)); % Generate the same style as B
Xstar = X1./sqrt(sum(X1.*X1));% normalize every columns; now, Xstar is Nonnegative_Stiefel
% check: % isNonnegSt(Xstar)

% step 2. generate matrix L as Proposition 1 in Jiang2022
% dk =.5+3*rand(k,1); % vector d
% xi=.8;
% L = xi*((dk*dk').^.5).*rand(k,k);
% L(sub2ind([k,k],1:k,1:k))=dk; % reset L's diagonal as vector d

L = rand(k,k);
L = L + k*eye(k);

% step 3.
C=Xstar*L';

% guess solution
problem.X_sol=Xstar;
problem.X_sol_q=norm(Xstar-C,"fro");
problem.C=C;

%%
problem.M=obliquefactory(n,k);

problem.cost=@(X) -2*trace(X'*C);
problem.egrad=@(x) -2*C;
ZERO = zeros(n,k);
problem.ehess = @(x,dx) ZERO;

problem.Euc_ineq = euclideanfactory(n, k);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
ZERO = zeros(n,k);
problem.ehess_barGx = @(x, z, dx) ZERO;

problem.Euc_eq = euclideanfactory(1, 1);
problem.eq_con=@(X) norm(X*V,'fro')^2-1;
problem.barHx=@(X,Y) (2*Y)*(X*VVt);
problem.barHxaj=@(X,dx) 2*F_inner(X*VVt,dx);
problem.ehess_barHx=@(X,Y,dx) (2*Y)*(dx*VVt);

% RIPM_checkupto2ndorder(problem);

%% initial point X on manifold M

x0=Proj_Stiefel(C);

%
% options.important=0.25;

options.KrylovIterMethod = 1;
options.maxiter = 300;

[xfinal, costfinal, residual, info, options]= RIPM(problem, x0, options);

if residual <= options.tolKKTres
    success=1;
else
    success =0;
end


%% Display some statistics.

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
    dist_to_sol(i)= abs( norm(info(i).xcurrent-problem.C,"fro")/problem.X_sol_q - 1);
end

subplot(1,3,3)
semilogy([dist_to_sol], '.-');
xlabel('Iteration number');
ylabel('dist to sol');


end

%%
