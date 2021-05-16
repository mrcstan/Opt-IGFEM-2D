 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/25/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%function Main
clear all;
close all;
format long e;
clc
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
inputFile = 'square_PDMS.in';
channelFile = 'DK_DP_check.channel';
costFuncType = 'P_NORM';
costFuncNormp = 1;

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
senstvtyTriNpt1D = 4; % number of gauss points for line quadrature in sensitivity analysis
senstvtyTriNpt2D = 73; % number of gauss points for triangular quadrature in sensitivity analysis
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
switch costFuncType
    case 'P_NORM'     
        sensitivity.costFunction.objOpt.normp = costFuncNormp;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
    case 'VARIANCE'       
        sensitivity.costFunction.objOpt.normp = 2;
        sensitivity.costFunction.objOpt.calcOneNorm = true;
        sensitivity.costFunction.type = costFuncType; 
        sensitivity.costFunction.area = mesh.domainArea;         
    case 'PRESSURE'
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = costFuncType;
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
    otherwise
        error('unrecognized cost function type')
end


[channels,sensitivity.designParams, ...
 sensitivity.restrictedParams]...
    = preprocess_channels(channelFile);

sensitivity.costFunction.objOpt.nDesignParams = sensitivity.designParams.nParams;
 %%  anonymous functions 
% create an implicit function to pass the FEM data as the objective fucntion to fmincons
% "@ (param_val)" is the parameter to be optimized, the other parameters are
% only needed to calculate the cost function and are not our primary variables
iniGuess = zeros(sensitivity.designParams.nParams,1); % initial guess for the change in design parameter

costfun = @ (delParams,calcGrad,UUR) check_DK_DP(delParams,  ...
                                    mesh, ...
                                    channels, ...
                                    femParams, ...
                                    sensitivity, ...
                                    calcGrad, UUR);

% create an implicit function to pass the FEM data as the constraint fucntion to fmincons
%nlconfun = @ (delParams) nonlinear_constraints(delParams, ...
%                                               sensitivity, ...
%                                               channels);

del = 1e-6;
method = 'central';


xo = zeros(sensitivity.designParams.nParams,1);
fprintf('\nreference K and P \n');
costfunTimer = tic;

[Kold,Pold,UURold] = costfun(xo,true,[]); 

toc(costfunTimer);
 

%
if strcmpi(method,'forward')
    FD_DK = nan(size(Kold));
    FD_DP = nan(size(Pold));
elseif strcmpi(method,'central')
    CD_DK = nan(size(Kold));
    FD_DP = nan(size(Pold));
else
    error('unrecognized method')
end

if (strcmpi(method,'forward'))
    costfunTimer = tic;
    xdel = xo;
    xdel = xdel + del; 
    [Knew,Pnew] = costfun(xdel,false,UURold);
    FD_DK = (Knew-Kold)/del;
    FD_DP = (Pnew-Pold)/del;
    toc(costfunTimer);
elseif (strcmpi(method,'central'))
    costfunTimer = tic;
    xdel = xo;
    xdel = xdel-del; 
    [Knew1,Pnew1] = costfun(xdel,false,UURold);
    toc(costfunTimer);
    
    costfunTimer = tic;
    xdel = xo;
    xdel = xdel+del; 
    [Knew2,Pnew2] = costfun(xdel,false,UURold);
    toc(costfunTimer);
    
    CD_DK = 0.5*(Knew2 - Knew1)/del;
    CD_DP = 0.5*(Pnew2 - Pnew1)/del;
end

format shortG

if strcmpi(method,'forward')
    disp('FD DK = ')
    disp(FD_DK)
    disp('FD DP = ')
    disp(FD_DP)
elseif strcmpi(method,'central')
    disp('CD DK = ')
    disp(CD_DK)
    disp('CD DP = ')
    disp(CD_DP)
end
%
preprocess_channels(channelFile);
fclose ('all');


      



