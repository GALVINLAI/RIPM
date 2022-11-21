function RIPM_checkupto2ndorder(problem)

condet = RIPM_constraintsdetail(problem);

%f(x)
fprintf('\n **************Starting check gradient of Objective function**************\n');
checkgradient(problem);
pause;
fprintf('\n **************Starting check hessian of Objective function**************\n');
checkhessian(problem);
pause;

%Inequality constraints
z = problem.Euc_ineq.rand();
problem_ineq_con.M=problem.M;
problem_ineq_con.cost = @(x) F_inner(problem.ineq_con(x), z);
problem_ineq_con.egrad = @(x) problem.barGx(x, z);
problem_ineq_con.ehess = @(x, dx) problem.ehess_barGx(x, z, dx);
fprintf('\n **************Starting check gradient of Inequality constraints**************\n');
checkgradient(problem_ineq_con);
pause;
fprintf('\n **************Starting check gradient of Inequality constraints**************\n');
checkhessian(problem_ineq_con);
pause;


if condet.has_ineq_cost
    %Equality constraints
    problem_eq_con.M=problem.M;
    y = problem.Euc_eq.rand();
    problem_eq_con.cost = @(x) F_inner(problem.eq_con(x), y);
    problem_eq_con.egrad = @(x) problem.barHx(x, y);
    problem_eq_con.ehess = @(x, dx) problem.ehess_barHx(x, y, dx);
    fprintf('\n **************Starting check gradient of Equality constraints**************\n');
    checkgradient(problem_eq_con);
    pause;
    fprintf('\n **************Starting check gradient of Equality constraints**************\n');
    checkhessian(problem_eq_con);
    pause;
end


end