clear; % Clear all variables from the workspace.
clc; % Clear the Command Window.
close all; % Close all open figures.

%% algorithm and problem
% This is the Local Newton RIPM in our paper.
%
% solver: Algorithm 2 - Prototype Algorithm of RIPM
% Because it is a local convergent algorithm, it will diverge if inital
% point not close to solution point sufficiently.
%
% problem: linear function on nonnegative sphere
%   f = a'*x (a linear objective)
%   M = sphere manifold
%   g = nonnegative constraint, i.e., g(x) = -x <= 0.

maxiter=300; % Set the maximum number of iterations for the algorithm to 300.
tol=1e-15; % Set the tolerance level for convergence to 1e-15 (a very small number close to machine precision).

%% setting
% Define the dimension of the problem as 100.
n=100; 

% Create a vector 'a' with its first element as -1 and the rest as absolute values of random numbers.
% note that objective function f = a'*x
a= [-1; abs(rand(n-1,1))]; 

% Define the true solution as a vector with 1 as the first element and zeros for the rest.
x_sol = [1; zeros(n-1, 1)]; 

% Create a sphere manifold of dimension 'n' using the 'spherefactory' function.
M = spherefactory(n); 

% Generate a random point on the manifold 'M' and take its absolute value to ensure non-negativity.
x = abs(M.rand()); 

% initial Lagrange multipliers
s = ones(n, 1); % Initialize the Lagrange multiplier 's' as a vector of ones.
z = ones(n, 1); % Initialize the Lagrange multiplier 'z' as a vector of ones.

% some constant matrix and vector
I = eye(n); % Create an identity matrix of size 'n'.
e = ones(n, 1); % Create a vector of ones of size 'n'.
gamma_hat = 0.5; % Define a constant 'gamma_hat' with a value of 0.5.

% Define the KKT vector field 'F' as a function of x, s, and z.
F = @(x, s, z) ...
    [(I - x * x') * (a - z); 
    x - s;
    diag(z) * diag(s) * e];

NormF = @(x, s, z) norm(F(x, s, z), 'fro'); % Define a function to compute the Frobenius norm of F.

Fnorm = NormF(x, s, z); % Calculate the norm of F at the initial values of x, s, and z.

% parameters
mu = 1; % Initialize parameter 'mu'.
gamma = 0.5; % Initialize parameter 'gamma'.

Record_x = []; % Initialize an empty array to record values of 'x'.
Record_s = []; % Initialize an empty array to record values of 's'.
Record_z = []; % Initialize an empty array to record values of 'z'.
Record_Fnorm = []; % Initialize an empty array to record norms of 'F'.
Record_mu = []; % Initialize an empty array to record values of 'mu'.
Record_gamma = []; % Initialize an empty array to record values of 'gamma'.
Record_dist2sol = []; % Initialize an empty array to record distances to the solution.
k = 0; % Initialize iteration counter 'k'.

%% iteration

while true
    % Record status
    Record_x = [Record_x; x'];
    Record_dist2sol = [Record_dist2sol; M.dist(x, x_sol)];
    Record_s = [Record_s; s'];
    Record_z = [Record_z; z'];
    Record_Fnorm = [Record_Fnorm; Fnorm];
    Record_mu = [Record_mu; mu];
    Record_gamma = [Record_gamma; gamma];  

    % Construct the left hand side of the perturbed Newton equation
    LeftSide=[
        x',         zeros(1,n), zeros(1,n);
        (z-a)'*x*I, zeros(n,n), -(I-x*x');
        I,          -I,         zeros(n,n);
        zeros(n,n), diag(z),    diag(s)
    ];

    % Construct the right hand side of the perturbed Newton equation
    RightSide = [
        0;
        (I - x * x') * (z - a);
        s - x;
        mu * e - diag(z) * diag(s) * e % mu*e is the perturbed term
    ];

    % Solve the perturbed Newton equation for dw
    dw = LeftSide \ RightSide; % Solve for dw such that LeftSide*dw=RightSide.

    % Extract increments for x, s, and z
    dx = dw(1:n);
    ds = dw(n + 1:2 * n);
    dz = dw(2 * n + 1:3 * n);

    % Compute the step size to ensure non-negativity of (z, s)
    InfBox = ones(1, n) * inf;
    for i = 1:n
        if ds(i) < 0
            InfBox(i) = -s(i) / ds(i);
        end
    end
    minbox_s = min(InfBox);

    InfBox = ones(1, n) * inf;
    for i = 1:n
        if dz(i) < 0
            InfBox(i) = -z(i) / dz(i);
        end
    end
    minbox_z = min(InfBox);

    % Determine the actual step size
    stepsize = min([1, gamma * minbox_s, gamma * minbox_z]);

    % Update the iterate points
    x = M.retr(x, dx, stepsize);
    s = s + stepsize * ds;
    z = z + stepsize * dz;
    Fnorm = NormF(x, s, z); % Update the Frobenius norm of F.

    % Check for convergence or reaching the maximum number of iterations
    if Fnorm < tol || k == maxiter
        break;
    end

    % Update parameters for the next iteration
    mu = min([mu / 1.5, 0.5 * Fnorm^2]);
    gamma = max(gamma_hat, 1 - Fnorm);
    % Alternatively, gamma could be updated as: gamma = 0.5 * (1 + max(0, 1 - Fnorm));
    k = k + 1; % Increment the iteration counter.
end

%% plot the iterations

% Define the position [left, bottom, width, height] in pixels
figure('Position', [100, 100, 800, 600]); % Adjust the values as needed

% Subplot 1
subplot(2, 2, 1);
semilogy(Record_dist2sol);
title('Distances to the solution');

% Subplot 2
subplot(2, 2, 2);
semilogy(Record_Fnorm);
title('Norm of KKT vector field');

% Subplot 3
subplot(2, 2, 3);
semilogy(Record_mu);
title('mu');

% Subplot 4
subplot(2, 2, 4);
plot(Record_gamma);
title('gamma');

% Adjust the layout of the subplots
sgtitle('Four Values as per Iterations');



