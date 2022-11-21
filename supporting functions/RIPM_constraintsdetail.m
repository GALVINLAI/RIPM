function condet = RIPM_constraintsdetail(problem)
% It takes in an problem and returns a struct with the
% following fields:
%
% Booleans:
% has_ineq_cost: True if problem has inequality cost functions. (e.g)
% has_eq_cost
%
% Int:
% n_ineq_constraint_cost: Number of cost function in inequality constraints.
% n_eq_constraint_cost
%
% It displays warning if the number of cost and grad functiosn do not match.
%

condet.has_ineq_cost = isfield(problem, 'ineq_con');
condet.has_eq_cost = isfield(problem, 'eq_con');

if condet.has_ineq_cost
    condet.n_ineq_constraint_cost  = problem.Euc_ineq.dim();
else
    condet.n_ineq_constraint_cost = 0;
end

if condet.has_eq_cost
    condet.n_eq_constraint_cost  = problem.Euc_eq.dim();
else
    condet.n_eq_constraint_cost = 0;
end

end