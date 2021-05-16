%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 10/22/2014
%%% Modified on 9/8/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function provides the objective function and also its gradient wrt
% design parameters if 'GradObj' is set to 'on' in the optimization options
% INPUT:
%   costFunction.type:
%       'DOMAIN_P_NORM_TEMP'
%       'VARIANCE'
%       'PRESSURE'
%   costFunction.area: area of domain (only valid if 
%                       costFunction.type == 'VARIANCE'
%   costFunction.objOpt: for mx_assemble_pseudo_adjoint_forces.cpp
%   costFunction.objOpt.normp: the p value in the p-norm
%   costFunction.objOpt.nDesignParams: number of design parameters
%   UUR: the solution at all nodes. this may also include Lagrange
%        multipliers at the end
%   UUR_update: the updated solution at all nodes, i.e., offset+backgroun
%               values. this may also include Lagrange multiplers at the
%               end
%   nLagrangeMult: number of Lagrange multipliers
%   if gradient is requested,
%       gauss:
%       KPP, KPF, KFP, KFF
%       fGloInd: logical indices of the free nodes and the Lagrange
%                multiplier
%       pGloInd: logical indices of the prescribed nodes
function [theta, gradTheta, oneNorm, gradOneNorm]...
    = objective_n_gradient(mesh, ...
                           channels, ...
                           costFunction, ...                                                   
                           UUR, ...
                           UUR_update, ...
                           nLagrangeMult, ...
                           gauss, ...
                           supg, ...
                           Kff, ...
                           fGloInd, ...
                           calcGrad, ...
                           calcOneNorm)
%{
if (costFunction.type == 0 && ~costFunction.isVariance) % DOMAIN_P_NORM_TEMP                                           
    normalization = max(UUR_update); % use max temperature as normalization factor
    UUR = UUR/normalization;
    mesh.elem.heatSource = mesh.elem.heatSource/normalization;
    mesh.convect.Tref = mesh.convect.Tref/normalization;
else
    normalization = 1;       
end
%}

normalization = 1;  

switch costFunction.type
    case {'DOMAIN_P_NORM_TEMP','VARIANCE'}
    if (calcGrad)
        [Fpseudo,Fadj,gradThetaRem,objVal,Fadj1norm,grad1normRem,oneNorm] ...
            = mx_assemble_pseudo_adjoint_forces(costFunction.objOpt,...
                                                mesh.node.coords',...
                                                mesh.elem.elem_node',...
                                                mesh.elem.heatSource,...
                                                mesh.convect,...
                                                gauss,...
                                                mesh.elem.parent,...
                                                channels,...
                                                mesh.elem.Neumann,...
                                                UUR, ...
                                                supg);
        fGloInd(end-nLagrangeMult+1:end) = []; % remove lagrange multipler indices
        Fpseudof = Fpseudo(:,fGloInd);
        Fadjf = Fadj(fGloInd);
        % these are zero since there is no boundary
        % term in the objective function
        Lambdaf = Kff' \ [Fadjf;zeros(nLagrangeMult,1)]; % column vector
        Lambdaf((end-nLagrangeMult+1):end) = [];           
        if (strcmpi(costFunction.type,'VARIANCE'))
            theta = objVal/costFunction.area - (oneNorm/costFunction.area)^2;
            Fadj1normf = Fadj1norm(fGloInd);
            Lambda1normf = Kff' \ [Fadj1normf;zeros(nLagrangeMult,1)]; % column vector
            Lambda1normf((end-nLagrangeMult+1):end) = [];
            gradTheta = (Fpseudof*Lambdaf + gradThetaRem)/costFunction.area ...
                       - 2/costFunction.area^2*oneNorm*(Fpseudof*Lambda1normf + grad1normRem);
        else
            % DOMAIN_P_NORM_TEMP
            theta = normalization...
                    *objVal^(1./costFunction.normp);
            gradTheta = (Fpseudof*Lambdaf + gradThetaRem)...
                        *theta/(objVal*costFunction.normp);
        end
    
        if (calcOneNorm)
            Fadjf = Fadj1norm(fGloInd);
            % these are zero since there is no boundary
            % term in the objective function
            Lambdaf = Kff' \ [Fadjf;zeros(nLagrangeMult,1)]; % column vector
            Lambdaf((end-nLagrangeMult+1):end) = [];
            oneNorm = normalization * oneNorm;
            gradOneNorm = (Fpseudof*Lambdaf + grad1normRem);

        else
            oneNorm = [];
            gradOneNorm = [];
        end

    else
        [~,~,~,objVal,~,~,oneNorm] ...
            = mx_assemble_pseudo_adjoint_forces(costFunction.objOpt,...
                                                mesh.node.coords',...
                                                mesh.elem.elem_node',...
                                                mesh.elem.heatSource,...
                                                mesh.convect,...
                                                gauss,...
                                                mesh.elem.parent,...
                                                channels,...
                                                mesh.elem.Neumann,...
                                                UUR, ...
                                                supg);                                    
        %UUR_update(end-nLagrangeMult+1:end) = [];
        if (costFunction.type == 0) 
            %theta = field_p_norm(mesh.elem,mesh.node.coords,UUR,...
            %                      costFunction.normp,gauss);
            if (costFunction.isVariance)
                theta = objVal/costFunction.area - (oneNorm/costFunction.area)^2;
            % DOMAIN_P_NORM_TEMP    
            else
                theta = normalization...
                        *objVal^(1./costFunction.normp);  
            end
        elseif (costFunction.type == 1) % DOMAIN_AVE_TEMP
            theta = average_temp(mesh.elem,mesh.node.coords,UUR_update);
        elseif (costFunction.type == 2) % COMPLIANCE
            theta = compliance(mesh.elem,mesh.node.coords,UUR,gauss);
        else
            error('unknown cost function type')
        end
        if (calcOneNorm)
            oneNorm = normalization * oneNorm;       
        else
            oneNorm = [];
        end
        gradOneNorm = [];
        gradTheta = [];
    end
    otherwise 
        error('unknown cost function type')
end
end

