function Stiefel_linear_or_projection_nonnega(n,k,linear_modle,initialflag)
% solver: RIPM_ineq
% name: linear/projection function on nonngegative orthogonality constraints

% f = -2*trace(X'*C) % f is linear function, if linear_modle=ture
% or, f = norm(X-C,"fro").^2  % f is projection function, if linear_modle=false

% M = Stiefel(n,k)
% g = nonngegative

% ture solution is Xstar


%% If no input is provided, we generate a quick demo.
if nargin == 0
    n=50;k=10;
    linear_modle=1;
    initialflag=2;
end

%% data setting


% construct data C
% step 1. generate a Nonnegative_Stiefel Xstar.
B=rand_Nonneg_St(n,k);% random B
X1= (B>0).*(1+rand(n,k)); % Generate the same style as B
Xstar = X1./sqrt(sum(X1.*X1));% normalize every columns; now, Xstar is Nonnegative_Stiefel
% check: % isNonnegSt(Xstar)

L = rand(k,k);
L = L + k*eye(k);

% step 3.
C=Xstar*(L');

% guess solution
problem.X_sol=Xstar;
problem.X_sol_q=norm(Xstar-C,"fro");
problem.C=C;

%%
problem.M=stiefelfactory(n,k);

% choose f(x)
if linear_modle
    % f is linear function, if linear_modle=ture
    problem.cost=@(X) -2*trace(X'*C);
    problem.egrad=@(x) -2*C;
    ZERO = zeros(n,k);
    problem.ehess = @(x,dx) ZERO;
else
    % f is projection function, if linear_modle=false
    problem.cost=@(X) norm(X-C,"fro").^2;
    problem.egrad=@(X) 2*(X-C);
    problem.ehess=@(X,dx) 2*dx;
end

problem.Euc_ineq = euclideanfactory(n, k);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
ZERO = zeros(n,k);
problem.ehess_barGx = @(x, z, dx) ZERO;


%% initial point X on manifold M

switch initialflag
    case 1
        x0=rand_Nonneg_St(n,k);

    case 2 %%%%%%%
        x0=Proj_Stiefel(C);

    case 3
        [x0,flag]=rounding(C);
        if ~flag
            return;
        end
end



%%

options.checkNTequation=0;
options.KrylovIterMethod=1;
options.maxiter=150;

[xfinal, costfinal, phifinal, info, options]= RIPM(problem, x0, options);


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
    % dist_to_sol(i)= norm(info(i).xcurrent-problem.C,"fro");
end

subplot(1,3,3)
semilogy([dist_to_sol], '.-');
xlabel('Iteration number');
ylabel('dist to sol');



end




