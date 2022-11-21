function [w, cost, info , options] = Sp_quadratic_nonnega()

% solver: RIPM_ineq
% name: quadratic function on nonngegative sphere
% or, Non-negative PCA

% f=quadratic function
% M=sphere
% g=nonngegative

% we do not know the true solution.


%% data setting

dim_set = [10, 50, 200, 500, 1000, 2000];  % Dimension of "the Cov Matrix"
snrset = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0]; % Signal Strength
deltaset = [0.1, 0.3, 0.7, 0.9];           % Sparsity


dim = 100; % d=n
snr = 0.5;
delta = 0.1;

%_______Set up data______
T = dim; % T is n
samplesize = floor(delta*dim); % cardinality |S| 
S = randsample(dim, samplesize); % support S 
v = zeros(dim,1); % ture principal direction v_{0}
v(S) = 1/sqrt(samplesize);
X = sqrt(snr) * v * (v.'); % first part of Z
Z = randn(dim)/sqrt(T); % random symmetric noise matrix N % but why not nonsymetry?
for ii = 1: dim
    Z(ii,ii) = randn * 2/sqrt(T);
end

A=-2*Z;

%% 
n=size(A,1);
problem.M = spherefactory(n);
problem.cost=@(x) .5*x'*A*x;
problem.egrad=@(x) .5*A*x+.5*A'*x;% NOTICE, egrad is A*x if A is symmetric;
problem.ehess =@(x,dx) .5*A*dx+.5*A'*dx;

problem.Euc_ineq = euclideanfactory(n);
problem.ineq_con=@(x) -x;
problem.barGx=@(x,z) -z;
problem.barGxaj=@(x,dx) -dx;
ZERO = zeros(n,1);
problem.ehess_barGx = @(x, z, dx) ZERO;


% initial point X on manifold M
x0=abs(problem.M.rand());

options.checkNTequation=0;
options.KrylovIterMethod=0;
options.maxiter=50;


%RIPM_checkGradientHessian_ineq(problem);
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
