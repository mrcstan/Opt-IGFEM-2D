function test_fmincon_irregular_point()
options = optimoptions(@fmincon, ...
        'FunValCheck', 'off',... % check for complex, inf or nan obj fun val and terminate if true   
        'GradObj', 'on',...     % supply gradient for the objective funcitons
        'GradCon', 'on', ...    % supply gradient for the constraint funcitons
        'DerivativeCheck', 'off', ... %Compare user-supplied derivatives (gradients of objective or constraints) to finite-differencing derivatives. 
        'FinDiffType','central', ...
        'Algorithm', 'sqp',...  % Choose the algorithm to solve the optimization
        'TolCon', 1e-6,...   % tolerance for constraint
        'TolX',1e-6,...      % tolerance for parameter value (multiplied by velocity) 
        'TolFun', 1e-6, ... % tolerence on the funciton we optimize
        'InitBarrierParam',0.01, ...% valid only for interior point algorithm 
        'Display', 'iter-detailed');   

iniGuess = -ones(2,1);
lb = -inf(2,1);
ub = inf(2,1);
[param_val, fval, exitflag, output, lambda, grad] ...
        = fmincon (@testObjFun,...    % Pass the objective function
                   iniGuess,...    % optimization parameter initial value (x0) (multiplied by velocity). x0 can be a scalar, vector, or matrix.
                   [], [],...      % the inequalities: subject to the linear inequalities A*x <= b. Here x is param_val.
                   [], [],...      % the linear equalities Aeq*x = beq. Here x is param_val.
                   lb, ub,...      % bound for the param_val
                   @testNlConFun,...   % the nonlinear constriant: the nonlinear inequalities c(x) or equalities ceq(x) defined in nonlcon. 
                   options);       % "sensitivity.options" which was defined earlier
param_val
fval
exitflag
%output
%lambda
%grad
end

function [val,grad] = testObjFun(X)
val =   X(1) + X(2);
if (nargout > 1)
    grad = [1;1];
end
end

function [g,geq,gradg,gradgeq] = testNlConFun(X)
g = [];
geq = nan(2,1);
geq(1) = (X(1) - 1)^2 + X(2)^2 - 1;
geq(2) = (X(1) - 2)^2 + X(2)^2 - 4;
if (nargout > 2)
    gradg = [];
    gradgeq = nan(2,2);
    gradgeq(:,1) = [2*(X(1) - 1); 2*X(2)];
    gradgeq(:,2) = [2*(X(1) - 2); 2*X(2)];
end
end

