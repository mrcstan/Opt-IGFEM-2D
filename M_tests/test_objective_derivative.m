%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/25/2014
%%% Last modified on 9/8/2015
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%function Main
clear all;
close all;
format long e;
% global variables used all over the code
global G_funcEvalCounter
global G_simDirectory
global G_filePrefix
G_simDirectory = './';
G_filePrefix = '2be_deleted';
G_funcEvalCounter = 0;

path(path, '../../NURBS/nurbs_toolbox')
path(path, '../../IGFEM-Curves-2D/M_geom_toolbox')
path(path, '../../IGFEM-Curves-2D/M_preFEM')
path(path, '../../IGFEM-Curves-2D/M_channels')
path(path, '../../IGFEM-Curves-2D/M_FEM')
path(path, '../../IGFEM-Curves-2D/mx_FEM')
path(path, '../../IGFEM-Curves-2D/M_postprocessing')
path(path, '../M_optimization')
path(path, '../M_opt_postprocessing')
path(path, '../mx_sensitivity')
path(path, '../../SISL/mx_SISL')
path(path, '../../IGFEM-Curves-2D/ChannelFiles')
path(path, '../InputFiles')
%% MESH AND USER INPUT
inputFile = 'DK_DP_check_linear.in';
%channelFile = 'DK_DP_check_nonlinear.channel';
channelFile = 'DK_DP_check_branching_8params_linear.channel';
%costFuncType = 'P_NORM_CHANNEL_W_OFFSET';
%costFuncType = 'NODAL_T_IN';
%costFuncType = 'P_NORM';
costFuncType = 'PRESSURE';
costFuncNormp = 8;

[mesh,femParams.gauss,...
 femParams.tol,...
 femParams.refine,...
 femParams.otherFlags,~, ...
 femParams.moveNode] = read_inputs(inputFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other user inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh.heatSourceFunc = [];
minPolyAngle = 5; % in degrees
minPolyArea = 0.001 * mesh.domainArea;
minSidePolyAngle = 0.5; % in degrees
minSidePolyArea = 0.001* mesh.domainArea;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENSITIVITY ANALYSIS INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
senstvtyTriNpt1D = 10; % number of gauss points for line quadrature in sensitivity analysis
senstvtyTriNpt2D = 16; % number of gauss points for triangular quadrature in sensitivity analysis
senstvtyQuadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM sensitivity analysis 
%% NO USER INPUT BELOW THIS LINE
%% Initialize required FEM data that does not change from one simulation to the next
sensitivity.gauss.line = gauss_points_and_weights(true,senstvtyTriNpt1D,1,'combined');
sensitivity.gauss.elem = gauss_points_and_weights(true,senstvtyTriNpt2D,2,'combined');
if (femParams.otherFlags.polyIGFEM)
    sensitivity.gauss.quadElem = gauss_points_and_weights(false,senstvtyQuadNpt1D,2,'combined');
end

% label edges for refinement later. Important!!!!!!!!!!
if femParams.refine.maxRefineLevel
    mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 
end
femParams.tol.halfLineWidth = mesh.edge.minLength*femParams.tol.halfLineWidthFrac;
femParams.moveNode.dist = mesh.edge.minLength*femParams.moveNode.distFrac;

%sensitivity.options.TypicalX = 0.1*ones(sensitivity.designParams.nParams,1);
% change sensitivity parameters
% check optimization inputs
sensitivity.costFunction.objOpt.node = []; % outlet node not requested
% intDomainType indicates whether the primary objective function is
% integrated over the whole domain. Note that this does not apply to the
% second possible objective function, i.e, the 1-norm.
% 0 = whole domain, 1 = channels only
sensitivity.costFunction.objOpt.intDomainType = 0;
switch costFuncType
    case 'P_NORM'     
        sensitivity.costFunction.objOpt.normp = costFuncNormp;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
    case 'P_NORM_CHANNEL'     
        sensitivity.costFunction.objOpt.normp = costFuncNormp;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType; 
        sensitivity.costFunction.objOpt.intDomainType = 1; 
    case 'P_NORM_CHANNEL_W_OFFSET'     
        sensitivity.costFunction.objOpt.normp = costFuncNormp;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType; 
        sensitivity.costFunction.objOpt.intDomainType = 1; 
        sensitivity.costFunction.objOpt.offset = 273.15;
        %sensitivity.costFunction.objOpt.offset = 0.0;
    case 'VARIANCE'       
        sensitivity.costFunction.objOpt.normp = 2;
        sensitivity.costFunction.objOpt.calcOneNorm = true;
        sensitivity.costFunction.type = costFuncType; 
        sensitivity.costFunction.area = mesh.domainArea;         
    case 'PRESSURE'
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
        sensitivity.costFunction.objOpt.intDomainType = 1;
    case 'PRESSURE_AVE_PANEL'
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
        sensitivity.costFunction.objOpt.intDomainType = 0;   
        sensitivity.costFunction.area = mesh.domainArea;     
    case 'AREA'
        sensitivity.costFunction.type = costFuncType;
    case 'TPA'
        sensitivity.costFunction.objOpt.normp = costFuncNormp;
        sensitivity.costFunction.objOpt.calcOneNorm = true;
        sensitivity.costFunction.type = costFuncType;
        sensitivity.costFunction.area = mesh.domainArea;
        if numel(costFuncWeights) ~= 3
            error('Three weights and reference scales must be provided')
        end
        sensitivity.costFunction.weights = costFuncWeights./costFuncScales;
    case 'NODAL_T_IN'
        % NOTE: only allows one node for now
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
        sensitivity.costFunction.objOpt.node = -1;
    case 'NODAL_T_OUT'
        % NOTE: only allows one node for now
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
        sensitivity.costFunction.objOpt.node = -2;
    case 'NODAL_T_INOUT'
        % NOTE: only allows one node for now
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
        sensitivity.costFunction.objOpt.node = -3; 
    otherwise
        error('unrecognized cost function type')
end

options.Toffset = mesh.convect.Toffset;
[channels,sensitivity.designParams, ...
 sensitivity.restrictedParams]...
    = preprocess_channels(channelFile,options);

sensitivity.costFunction.objOpt.nDesignParams = sensitivity.designParams.nParams;
 %%  anonymous functions 
% create an implicit function to pass the FEM data as the objective fucntion to fmincons
% "@ (param_val)" is the parameter to be optimized, the other parameters are
% only needed to calculate the cost function and are not our primary variables
iniGuess = zeros(sensitivity.designParams.nParams,1); % initial guess for the change in design parameter

costfun = @ (delParams) fem_n_sensitivity (delParams,  ...
                                           mesh, ...
                                           channels, ...
                                           femParams, ...
                                           sensitivity, ...
                                           1, false, false);

% create an implicit function to pass the FEM data as the constraint fucntion to fmincons
%nlconfun = @ (delParams) nonlinear_constraints(delParams, ...
%                                               sensitivity, ...
%                                               channels);

del = 1e-7;
method = 'central';
if strcmpi(method,'forward')  
    thetaNew = nan(sensitivity.designParams.nParams,1);
elseif strcmpi(method,'central')
    thetaNew = nan(sensitivity.designParams.nParams,2);
else
    error('unrecognized method')
end
FDtheta = nan(sensitivity.designParams.nParams,1);

xo = zeros(sensitivity.designParams.nParams,1);
fprintf('\nreference theta \n');
costfunTimer = tic;
[thetaOld,gradTheta] = costfun(xo);
toc(costfunTimer);


disp('analytical gradient')
disp(gradTheta)
%
for i = 1:sensitivity.designParams.nParams
    fprintf('\n-----------------------------------------------------------\n')
    fprintf('design param %i/%i',i,sensitivity.designParams.nParams);
    fprintf('\n-----------------------------------------------------------\n')
    if (strcmpi(method,'forward'))
        %costfunTimer = tic;
        xdel = xo;
        xdel(i) = xdel(i) + del; 
        thetaNew(i) = costfun(xdel);
        %[~,~,thetaNew(i)] = costfun(xdel);
        FDtheta(i) = (thetaNew(i)-thetaOld)/del;
        %toc(costfunTimer);
    elseif (strcmpi(method,'central'))
        for j = 1:2
            %costfunTimer = tic;
            xdel = xo;
            xdel(i) = xdel(i) + del*(-1)^j; 
            thetaNew(i,j) = costfun(xdel);
            %[~,~,thetaNew(i,j)] = costfun(xdel);
            %toc(costfunTimer);
        end
        FDtheta(i) = 0.5*(thetaNew(i,2)-thetaNew(i,1))/del;
    end
end



disp('absolute grad diff')
disp(abs(FDtheta-gradTheta))

disp('relative grad diff')
disp(abs(1-FDtheta./gradTheta))

fprintf('analytical grad theta\n')
fprintf('%11.9e \n',gradTheta)
fprintf('FD grad theta\n')
fprintf('%11.9e \n',FDtheta)
fclose ('all');
%

      



