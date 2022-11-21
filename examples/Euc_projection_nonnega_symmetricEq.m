function Euc_projection_nonnega_symmetricEq()

% solver: RIPM
% name: projection onto nonnegative symmetric matrices space

% f=.5*norm(X-A,"fro")^2
% M=Euc(n,n)
% h=symmetric X-X'=0
% g=nonngegative

% solution is [sym(A)]_{+}


%% data setting

n=3;
A=randn(n,n); % has both positive and negative entries.

% guess solution
problem.X_sol=Sym(A);
problem.X_sol(problem.X_sol<0)=0;

%%

problem.M =  euclideanfactory(n,n);
problem.cost=@(X) .5*norm(X-A,"fro")^2;
problem.egrad =@(X) X-A;
problem.ehess = @(x,dx) dx;

problem.Euc_ineq = euclideanfactory(n, n);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
ZERO = zeros(n,n);
problem.ehess_barGx = @(x, z, dx) ZERO;

problem.Euc_eq = skewsymmetricfactory(n);
problem.eq_con=@(X) X-X';
problem.barHx=@(X,Y) 2*Y;
problem.barHxaj=@(X,dx) dx-dx';
ZERO = zeros(n,n);
problem.ehess_barHx = @(x, y, dx) ZERO;


% initial point X on manifold M
x0=abs(problem.M.rand());

options.checkNTequation=1;
options.KrylovIterMethod=0;


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
dist_to_sol(i)=norm(info(i).xcurrent-problem.X_sol,'fro');
end

subplot(1,3,3)
semilogy([dist_to_sol], '.-');
xlabel('Iteration number');
ylabel('dist to sol');

end
