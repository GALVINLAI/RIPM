function checkconstraints_upto2ndorder(problem, x0, d)

    has_d = exist('d', 'var') && ~isempty(d);
    has_x = exist('x0', 'var')&& ~isempty(x0);
    has_ineq_cost = isfield(problem, 'ineq_constraint_cost');
    has_ineq_grad = isfield(problem, 'ineq_constraint_grad');
    has_ineq_hess = isfield(problem, 'ineq_constraint_hess');
    has_eq_cost = isfield(problem, 'eq_constraint_cost');
    has_eq_grad = isfield(problem, 'eq_constraint_grad');
    has_eq_hess = isfield(problem, 'eq_constraint_hess');
    
        if has_ineq_cost
            n_ineq_constraint_cost  = length(problem.ineq_constraint_cost);
        else
            n_ineq_constraint_cost = 0;
        end
        if has_ineq_grad
            n_ineq_constraint_grad  = length(problem.ineq_constraint_grad);
        else
            n_ineq_constraint_grad  = 0;
        end
        if has_ineq_hess
            n_ineq_constraint_hess  = length(problem.ineq_constraint_hess);
        else
            n_ineq_constraint_hess  = 0;
        end

        if has_eq_cost
            n_eq_constraint_cost  = length(problem.eq_constraint_cost);
        else
            n_eq_constraint_cost = 0;
        end
        if has_eq_grad
            n_eq_constraint_grad  = length(problem.eq_constraint_grad);
        else 
            n_eq_constraint_grad = 0;
        end
        if has_eq_hess
            n_eq_constraint_hess  = length(problem.eq_constraint_hess);
        else 
            n_eq_constraint_hess = 0;
        end

        if (n_ineq_constraint_cost ~= n_ineq_constraint_grad)
            warning('checkconstraints:number',['the number of cost functions of'...
                'inequality constraints do not match the number of gradient functions']);
        end
        if (n_ineq_constraint_grad ~= n_ineq_constraint_hess)
            warning('checkconstraints:number',['the number of grad functions of'...
                'inequality constraints do not match the number of hessian functions']);
        end
        if (n_ineq_constraint_cost ~= n_ineq_constraint_hess)
            warning('checkconstraints:number',['the number of cost functions of'...
                'inequality constraints do not match the number of hessian functions']);
        end    

        if (n_eq_constraint_cost ~= n_eq_constraint_grad)
            warning('checkconstraints:number',['the number of cost functions of'...
                'equality constraints do not match the number of gradient functions']);
        end
        if (n_eq_constraint_grad ~= n_eq_constraint_hess)
            warning('checkconstraints:number',['the number of grad functions of'...
                'equality constraints do not match the number of hessian functions']);
        end    
        if (n_eq_constraint_cost ~= n_eq_constraint_hess)
            warning('checkconstraints:number',['the number of cost functions of'...
                'equality constraints do not match the number of hessian functions']);
        end 

    
    for iter = 1:n_ineq_constraint_cost
        newproblem.M = problem.M;
        newproblem.cost = problem.ineq_constraint_cost{iter};
        newproblem.egrad = problem.ineq_constraint_grad{iter};
        newproblem.ehess = problem.ineq_constraint_hess{iter};
        if has_x
            if has_d
                figure;
                checkgradient(newproblem, x0, d);
                figure;
                checkhessian(newproblem, x0, d);
            else
                figure;
                checkgradient(newproblem, x0);
                figure;
                checkhessian(newproblem, x0);
            end
        else
            figure;
            checkgradient(newproblem);
            figure;
            checkhessian(newproblem);
        end
    end
    
    for iter = 1:n_eq_constraint_cost
        newproblem.M = problem.M;
        newproblem.cost = problem.eq_constraint_cost{iter};
        newproblem.egrad = problem.eq_constraint_grad{iter};
        newproblem.ehess = problem.eq_constraint_hess{iter};
        if has_x
            if has_d
                figure;
                checkgradient(newproblem, x0, d);
                figure;
                checkhessian(newproble, x0, d);
            else
                figure;
                checkgradient(newproblem, x0);
                figure;
                checkhessian(newproblem, x0);
            end
        else
            figure;
            checkgradient(newproblem);
            figure;
            checkhessian(newproblem);
        end
    end  
end