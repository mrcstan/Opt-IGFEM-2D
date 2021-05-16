%%% Created by Marcus Tan on 10/11/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. 
%%% All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function returns the value and gradient of the cost function
%%% of the nonlinear program equivalent to the minimax problem
function [f,gradf] = minimax_param(x,scale)
f = x(end)/scale;
if nargout > 1
    gradf = [zeros(numel(x)-1,1);1/scale];
else
    gradf = [];
end
end