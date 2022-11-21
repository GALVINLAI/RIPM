function condet = constraintsdetail(problem)
% It takes in an problem and returns a struct with the
% following fields:
% 
% Booleans:
% has_ineq_cost: True if problem has inequality cost functions. (e.g)
% has_ineq_grad
% has_ineq_hess
% has_eq_cost
% has_eq_grad
% has_eq_hess
% 
% Int:
% n_ineq_constraint_cost: Number of cost function in inequality constraints.
% n_ineq_constraint_grad
% n_ineq_constraint_hess
% n_eq_constraint_cost
% n_eq_constraint_grad
% n_eq_constraint_hess
% 
% It displays warning if the number of cost and grad functiosn do not match.
%
%
% This file is part of Manopt: www.manopt.org.
% Original author: Changshuo Liu, September 3, 2017.
% Modified author: Mitsuaki Obara, January 31, 2020.
%
% Change Log: Add hessian parts to the program (MO)


    condet.has_ineq_cost = isfield(problem, 'ineq_constraint_cost');
    condet.has_ineq_grad = isfield(problem, 'ineq_constraint_grad');
    condet.has_ineq_hess = isfield(problem, 'ineq_constraint_hess');
    condet.has_eq_cost = isfield(problem, 'eq_constraint_cost');
    condet.has_eq_grad = isfield(problem, 'eq_constraint_grad');
    condet.has_eq_hess = isfield(problem, 'eq_constraint_hess');
    
    if condet.has_ineq_cost
        condet.n_ineq_constraint_cost  = length(problem.ineq_constraint_cost);
    else
        condet.n_ineq_constraint_cost = 0;
    end
    if condet.has_ineq_grad
        condet.n_ineq_constraint_grad  = length(problem.ineq_constraint_grad);
    else
        condet.n_ineq_constraint_grad  = 0;
    end
    if condet.has_ineq_hess
        condet.n_ineq_constraint_hess  = length(problem.ineq_constraint_hess);
    else
        condet.n_ineq_constraint_hess  = 0;
    end
    
    if condet.has_eq_cost
        condet.n_eq_constraint_cost  = length(problem.eq_constraint_cost);
    else
        condet.n_eq_constraint_cost = 0;
    end
    if condet.has_eq_grad
        condet.n_eq_constraint_grad  = length(problem.eq_constraint_grad);
    else 
        condet.n_eq_constraint_grad = 0;
    end
    if condet.has_eq_hess
        condet.n_eq_constraint_hess  = length(problem.eq_constraint_hess);
    else 
        condet.n_eq_constraint_hess = 0;
    end
    
    if (condet.n_ineq_constraint_cost ~= condet.n_ineq_constraint_grad)
        warning('checkconstraints:number',['the number of cost functions of'...
            'inequality constraints do not match the number of gradient functions']);
    end
    if (condet.n_ineq_constraint_grad ~= condet.n_ineq_constraint_hess)
        warning('checkconstraints:number',['the number of grad functions of'...
            'inequality constraints do not match the number of hessian functions']);
    end
    if (condet.n_ineq_constraint_cost ~= condet.n_ineq_constraint_hess)
        warning('checkconstraints:number',['the number of cost functions of'...
            'inequality constraints do not match the number of hessian functions']);
    end    
    
    if (condet.n_eq_constraint_cost ~= condet.n_eq_constraint_grad)
        warning('checkconstraints:number',['the number of cost functions of'...
            'equality constraints do not match the number of gradient functions']);
    end
    if (condet.n_eq_constraint_grad ~= condet.n_eq_constraint_hess)
        warning('checkconstraints:number',['the number of grad functions of'...
            'equality constraints do not match the number of hessian functions']);
    end    
    if (condet.n_eq_constraint_cost ~= condet.n_eq_constraint_hess)
        warning('checkconstraints:number',['the number of cost functions of'...
            'equality constraints do not match the number of hessian functions']);
    end 
end