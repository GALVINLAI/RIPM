function  Euc_projection_nonnega()

% solver: RIPM_ineq
% name: projection onto nonnegative orthant

% f=.5*norm(X-A,"fro")^2;
% M=Euclidean space
% g=nonnegative

% true solution is [A]_{+}

%% data setting

m=5;
n=10;
A=randn(m,n); % has both positive and negative entries.

% true solution
problem.X_sol=A;
problem.X_sol(problem.X_sol<0)=0; % solution is [A]_{+}


%% 

problem.M=euclideanfactory(m, n);
problem.cost=@(X) .5*norm(X-A,"fro")^2;
problem.egrad=@(X) X-A;
problem.ehess = @(x,dx) dx;

problem.Euc_ineq = euclideanfactory(m, n);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
ZERO = zeros(m,n);
problem.ehess_barGx = @(x, z, dx) ZERO;

% initial point X on manifold M
x0=abs(problem.M.rand());

options.checkNTequation=0;
options.KrylovIterMethod=1;


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
