function  [x, xcost, all_info, outer_info, Success]=Rie_Barrier_simple(a,mu,x)

% Barrier Methods on Riemannian Optimization for a simple problem.
% For now, this method is tested only for inqualities constrains.
% min _{x \in \mathbb{S}^{n-1}} a^{T}x s.t. x >= 0.
% See Section 7.1 An Intuitive Barrier Method on Manifolds

% If no input is provided, generate random data for a quick demo
if nargin == 0
    % problem setting
    a=[-1 2 1]';
    n=size(a,1);
    mu=10;
    x=ones(n,1)*(sqrt(n)/n);
    sol_x = [1, 0, 0]'; % [1 0 0] is a solution.
elseif nargin ~= 3
    error('Please provide 6 inputs (or none for a demo).');
end

% Create the problem structure.
n=size(a,1);
problem.M = spherefactory(n);

options.verbosity = 0; % Change this number for more or less output

all_info=[];

outeriter = 0;
outer_stats.outeriter = outeriter;
outer_stats.current_x = x;
outer_info(outeriter+1) = outer_stats;

while 1

    problem.cost  = @(x) a'*x-mu*sum(reallog(x)); % logarithmic barrier function
    problem.egrad = @(x) a-mu.*diag(1./x)*ones(n,1); % euclidean gradient of cost
    problem.ehess = @(x, u) (mu.*diag(1./(x.^2)))*u; % euclidean hessian of cost

    % record at an iteration
    metrics.val_mu= @(problem, x) mu; % value of mu
    metrics.barrier_mulitplier = @(problem, x) mu./x; % barrier mulitplier
    metrics.dist = @(problem, x) norm(x-sol_x); %
    metrics.current_x = @(problem, x) x; % x itself
    options.statsfun = statsfunhelper(metrics);

    %Stopping criteria
    options.tolgradnorm = 1e-03;
    options.maxiter = 10; % max iters for each mu_k

    % Solve.
    [x, xcost, info, options] = trustregions(problem,x,options);
    all_info=[all_info, info];

    outeriter = outeriter + 1;
    outer_stats.outeriter = outeriter;
    outer_stats.current_x = x;
    outer_info(outeriter+1) = outer_stats;

    if mu < 1e-10 %stop condtion
        sprintf('Success!!')
        Success=1;
        break
    else
        mu=mu/1.5; % shrink mu
    end

    if size(all_info,2) >= 4000
        sprintf('Failed!!')
        Success=0;
        break
    end

end

if nargin == 0
    %Plot the details of each iteration.
    subplot(1,2,1)
    semilogy([all_info.gradnorm], '.-');
    xlabel('Iteration number');
    ylabel('Norm of the gradient of barrier function.');

    subplot(1,2,2)
    semilogy([all_info.dist], '.-');
    xlabel('Iteration number');
    ylabel('Distance to solution x^{*}.');
end

end

