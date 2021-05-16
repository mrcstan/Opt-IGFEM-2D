function test_fminimax
    global fcount 
    fcount = 0;
    
    %{
    x0 = [0.1; 0.1];       % Make a starting guess at solution
    options = optimoptions(@fminimax, ...
                            'FunValCheck', 'off',... % check for complex, inf or nan obj fun val and terminate if true   
                            'GradObj', 'on',...     % supply gradient for the objective funcitons
                            'GradConstr', 'off', ...    % supply gradient for the constraint funcitons
                            'DerivativeCheck', 'off', ... %Compare user-supplied derivatives (gradients of objective or constraints) to finite-differencing derivatives. 
                            'TolCon', 1e-6,...   % tolerance for constraint
                            'TolX',1e-6,...      % tolerance for parameter value (multiplied by velocity) 
                            'TolFun', 1e-6, ... % tolerence on the funciton we optimize
                            'Display', 'iter-detailed');   
    [x,fval] = fminimax(@myfun,x0,[],[],[],[],[],[],[],options);
    %}
    %
    x0 = [0.1; 0.1; 0.0];       % Make a starting guess at solution
    options = optimoptions(@fmincon, ...
                            'FunValCheck', 'off',... % check for complex, inf or nan obj fun val and terminate if true   
                            'GradObj', 'on',...     % supply gradient for the objective funcitons
                            'GradConstr', 'on', ...    % supply gradient for the constraint funcitons
                            'DerivativeCheck', 'off', ... %Compare user-supplied derivatives (gradients of objective or constraints) to finite-differencing derivatives. 
                            'TolCon', 1e-6,...   % tolerance for constraint
                            'TolX',1e-6,...      % tolerance for parameter value (multiplied by velocity) 
                            'TolFun', 1e-6, ... % tolerence on the funciton we optimize
                            'Algorithm', 'sqp',...  % Choose the algorithm to solve the optimization
                            'Display', 'iter-detailed');   
    [x,fval] = fmincon(@extraParam,x0,[],[],[],[],[],[],@nlconfun,options);
    %
    x
    fval
end

function [f,gradf] = myfun(x)
global fcount
fcount = fcount + 1;
% fminimax fails if gradobj = on and nan values are returned
%if fcount > 3 && fcount < 10
%    f = nan(1,5);
%    gradf = zeros(2,5);
%    return
%end
f(1)= 2*x(1)^2+x(2)^2-48*x(1)-40*x(2)+304;     % Objectives
f(2)= -x(1)^2 - 3*x(2)^2;
f(3)= x(1) + 3*x(2) -18;
f(4)= -x(1)- x(2);
f(5)= x(1) + x(2) - 8;
if nargout > 1
    gradf = nan(2,5);
    gradf(:,1) = [4*x(1)-48;2*x(2)-40];
    gradf(:,2) = [-2*x(1);-6*x(2)];
    gradf(:,3) = [1;3];
    gradf(:,4) = [-1;-1];
    gradf(:,5) = [1,1];
else
    gradf = [];
end
end

function [f,gradf] = extraParam(x)
f = x(3);
if nargout > 1
    gradf = [0;0;1];
else
    gradf = [];
end
end

function [f,feq,gradf,gradfeq] = nlconfun(x)
global fcount
fcount = fcount + 1;
% fminimax fails if gradobj = on and nan values are returned
if fcount > 3 && fcount < 5
    f = nan(1,5);
    feq = [];
    gradf = nan(3,5);
    gradfeq = [];
    return
end
f(1)= 2*x(1)^2+x(2)^2-48*x(1)-40*x(2)+304-x(3);     % Objectives
f(2)= -x(1)^2 - 3*x(2)^2-x(3);
f(3)= x(1) + 3*x(2) -18 -x(3);
f(4)= -x(1)- x(2) -x(3);
f(5)= x(1) + x(2) - 8  -x(3);
feq = [];
gradfeq = [];
if nargout > 2
    gradf = nan(3,5);
    gradf(:,1) = [4*x(1)-48;2*x(2)-40;-1];
    gradf(:,2) = [-2*x(1);-6*x(2);-1];
    gradf(:,3) = [1;3;-1];
    gradf(:,4) = [-1;-1;-1];
    gradf(:,5) = [1;1;-1];
else
    gradf = [];
end
end



