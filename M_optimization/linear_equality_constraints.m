%%% Created by Marcus Tan on 1/19/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function creates the matrix Aeq and the RHS vector beq for 
%%% imposing the constraints that pairs of design parameters are equal 
%%% Aeq and beq are meant for the linear equalities in fmincon
function [Aeq,beq] = linear_equality_constraints(restrictedParams,...
                                                 nDesignParams)
if (restrictedParams.nParams == 0)
    Aeq = [];
    beq = [];
    return
end
nRows = size(restrictedParams.paramPairs,1);                                             
beq = zeros(nRows,1);

nCols = nDesignParams + restrictedParams.nParams;
Aeq = full(sparse([1:nRows,1:nRows]',restrictedParams.paramPairs(:),...
                  [ones(nRows,1);-ones(nRows,1)],nRows,nCols));
end