%%% Created by Marcus Tan on 1/19/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lb,ub] = params_bounds(designParams,restrictedParams)
nParams = designParams.nParams + restrictedParams.nParams;
lb = -inf(nParams,1);
ub = inf(nParams,1);
lb(1:designParams.nParams) = designParams.bounds(:,1) - designParams.iniVals;
ub(1:designParams.nParams) = designParams.bounds(:,2) - designParams.iniVals;
if (restrictedParams.nParams)
    for i = 1:size(restrictedParams.paramPairs,1)
        lb(restrictedParams.paramPairs(i,2)) = lb(restrictedParams.paramPairs(i,1));
        ub(restrictedParams.paramPairs(i,2)) = ub(restrictedParams.paramPairs(i,1));
    end
end
end