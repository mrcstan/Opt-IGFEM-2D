%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/25/2014
%%% Last modified on 10/9/2015
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
path(path, '../../SISL')
path(path, '../../IGFEM-Curves-2D/ChannelFiles')
path(path, '../../IGFEM-Curves-2D/blockedChannelFiles')
%% MESH AND USER INPUT
% 2 choices: I) supply an Abaqus mesh file
%            II) give information for generating the mesh in this file
% In both cases, the BC's must be specify in this file as this code is
% capable of mesh refinement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I) Abaqus mesh file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% straight/cross channel
%meshfile = 'square_p1_40x40_patran_temp_all.inp';
%meshfile = 'square_p1_10x10C_abq_temp_all.inp';
%meshfile = 'square_p1_10x10_abq_temp_top_heatflux_bottom.inp';
%mesh = read_input(meshfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II) Information for generating mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% domain boundary (in meters)
mesh.boundary.xi = 0.0; 
mesh.boundary.xf = 0.05;
mesh.boundary.yi = 0.0;
mesh.boundary.yf = 0.05;
%
% element connectivity and node coordinates
[mesh.elem.elem_node,mesh.node.coords] ...
    = generate_uniform_mesh([7,7],... % number of elements in x- and y-directions
                            [mesh.boundary.xi,...
                             mesh.boundary.xf,...
                             mesh.boundary.yi,...
                             mesh.boundary.yf],...
                             2); % 1: diagonal along NW direction
                                 % 2: diagonal along NE direction

                                
% element material and material conductivity
%mesh.material.conductivity = 0.2*0.003; % epoxy
%mesh.material.conductivity = 0.00288; % epoxy
%mesh.material.conductivity = 2.9*0.003; % old composite value
mesh.material.conductivity = 2.7*0.003; % composite value Dec 9, 2014
%mesh.material.conductivity = 0.6;
mesh.elem.material = int32(ones(size(mesh.elem.elem_node,1),1));

% distributed heat source 
mesh.elem.heatSource = 6000.0*ones(size(mesh.elem.elem_node,1),1); 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III) Other user inputs that must be specified regardless of choice I or II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%channelFile = 'serpentine_simpleB_8_params.channel';
%channelFile = 'parallel2_test_pressure_constraint_w_diams.channel';
%channelFile = 'DK_DP_check_branchingA_2params.channel';
%blockedChannelFile = 'DK_DP_check_branchingA_2params.blk';
channelFile = 'parallel2x2_test_sensitivity.channel';
blockedChannelFile = 'parallel2x2nBlk1.blk';

% cost function types available are
%       i) P_NORM: costFuncNormp must be provided
%       ii) VARIANCE
%       iii) PRESSURE
%       iv) AREA
costFuncType = 'P_NORM';
costFuncNormp = 8;
costFuncWeights = [1.0,1.0,1.0]/3.0;
costFuncScales = [50,1e5,0.1];

minPolyAngle = 2; % in degrees
domainArea = (mesh.boundary.xf - mesh.boundary.xi) ...
                   *(mesh.boundary.yf - mesh.boundary.yi);
minPolyArea = 0.1^2*domainArea;

domainThickness = 3e-3; % in meters

% randomize initial values of design parameters subject to bounds
randomizeIniVals = false;                       
% boundary conditions
% 1 : x = xi, 2 : x = xf, 3 : y = yi, 4 : x = yf
mesh.BCs.boundaries = []; % totally insulated boundaries
%mesh.BCs.boundaries = [3,4]; 
%mesh.BCs.types = [2,1];
%mesh.BCs.values_or_funcs = {10,20};
%mesh.BCs.values_or_funcs = {2000,@(x,y) 20*x/0.1};
%mesh.BCs.values_or_funcs = {1,21.5};
%mesh.BCs.values_or_funcs = {u,u,u,u};   

% specify convective heat coefficient
%mesh.convect.coef = 8+5.6;
%mesh.convect.Tref = 21;
mesh.convect.coef = 0.0; 
mesh.convect.Tref = 0.0;

% for nodes_curves_intersect.m and move_nodes.m 
femParams.tol.node = 1e-6; % tolerance for checking if original nodes coincide with channels
femParams.moveNode.distFrac = 0.05; % distance to move node when it coincides with channel
                                    % or when a branching point or kink coincides
                                    % with an element edge, as a fraction of the
                                    % minimum edge length
femParams.moveNode.maxAttemps = 10;
femParams.moveNode.randDirection = false;
femParams.tol.boundary = 1e-13; % also for set_boundary_conditions.m

% tolerances for edges_curve_intersect.cpp
femParams.tol.nurbsParam = 1.e-8;
femParams.tol.epsco = 1e-15; % computational tolerance for curve_edge_intersect.mexw64
femParams.tol.epsge = 1e-6; % geometrical tolerance for curve_edge_intersect.mexw64
femParams.tol.cosAngleTol = cos(20*pi/180.0); % cosine of the min angle of the resulting element after edge flipping
femParams.tol.halfLineWidthFrac = 1e-4; % half-width fraction of channels
femParams.tol.vert = tan(89.9*pi/180); % tan of max positive angle wrt horizontal axis beyond which a line is considered vertical
femParams.opt.maxRefineLevel = 6;
femParams.opt.refineJuncElem = false; % force refinement of element with branching
femParams.tol.intersectEdges = 1e-13; % for intersect_edges in single_edge_curves_intersect.m

femParams.tol.channelSelfIntersect = 1e-13; % tolerance for channels_self_intersections

% tolerenace for approximating slender child element with triangle or quadrilateral
% slenderTol = []; % if approximation not desired
femParams.slenderTol.minAngle = 5; % min angle of a triangular child element below which NURBS is approximated by line segments
femParams.slenderTol.maxAspectRatio = inf;  % max aspect ratio of a quadrialteral child element below which NURBS is approximated by line segments

femParams.polyIGFEM = true; % use polynomial IGFEM
femParams.supg = true; % apply SUPG

%% NO USER INPUT BELOW THIS LINE

% gauss quadrature schemes
triNpt1D = 4;
triNpt2D = 7;
senstvtyTriNpt1D = 4; % number of gauss points for line quadrature in sensitivity analysis
senstvtyTriNpt2D = 16; % number of gauss points for triangular quadrature in sensitivity analysis
quaNpt1Dt = 4; % number of gauss points per knot span along channel
quaNpt1Dn = 4; % number of gauss points per knot span normal to channel
quadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM
senstvtyQuadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM sensitivity analysis 
%% NO USER INPUT BELOW THIS LINE
%% Initialize required FEM data that does not change from one simulation to the next
femParams.gauss.line = gauss_points_and_weights(true,triNpt1D,1,'combined');
femParams.gauss.elem = gauss_points_and_weights(true,triNpt2D,2,'combined');
sensitivity.gauss.line = gauss_points_and_weights(true,senstvtyTriNpt1D,1,'combined');
sensitivity.gauss.elem = gauss_points_and_weights(true,senstvtyTriNpt2D,2,'combined');
if (femParams.polyIGFEM)
    femParams.gauss.quadElem = gauss_points_and_weights(false,quadNpt1D,2,'combined');
    sensitivity.gauss.quadElem = gauss_points_and_weights(false,senstvtyQuadNpt1D,2,'combined');
else
    femParams.gauss.qua1Dt = gauss_points_and_weights(false,quaNpt1Dt,1,'separate'); % along channel
    femParams.gauss.qua1Dn = gauss_points_and_weights(false,quaNpt1Dn,1,'separate'); % normal to channel
end

domainVol =  (mesh.boundary.xf-mesh.boundary.xi) ...
            *(mesh.boundary.yf-mesh.boundary.yi)*domainThickness;

[channels,mesh.heatSourceFunc,...
 sensitivity.designParams,sensitivity.restrictedParams]...
    = preprocess_channels(channelFile);
blockedSets = read_blocked_channels(blockedChannelFile);
%sensitivity.options.TypicalX = 0.1*ones(sensitivity.designParams.nParams,1);
% change sensitivity parameters
switch costFuncType
    case 'P_NORM'     
        sensitivity.costFunction.objOpt.normp = costFuncNormp;
        sensitivity.costFunction.type = costFuncType;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
    otherwise
        error('unrecognized cost function type')
end

sensitivity.costFunction.objOpt.nDesignParams = sensitivity.designParams.nParams;
%% Initialize required FEM data that does not change from one function evaluation to the next
% label edges for refinement later. Important!!!!!!!!!!
mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 
% the edge_node information is generated by generate_conforming_mesh
%fprintf('generate edge_node \n')
mesh.edge.edge_node = find_edge_node(mesh.elem.elem_node);

mesh.edge.length = find_edge_length(mesh.edge.edge_node,mesh.node.coords);
mesh.edge.minLength = min(mesh.edge.length);
femParams.tol.halfLineWidth = mesh.edge.minLength*femParams.tol.halfLineWidthFrac;
femParams.moveNode.dist = mesh.edge.minLength*femParams.moveNode.distFrac;
 %%  anonymous functions 
% create an implicit function to pass the FEM data as the objective fucntion to fmincons
% "@ (param_val)" is the parameter to be optimized, the other parameters are
% only needed to calculate the cost function and are not our primary variables
iniGuess = zeros(sensitivity.designParams.nParams,1); % initial guess for the change in design parameter


costfun = @ (delParams) ...
                fem_n_sensitivity_blocked_channels (delParams,  ...
                                                    blockedSets, ...
                                                    mesh, ...
                                                    channels, ...
                                                    femParams, ...
                                                    sensitivity, ...
                                                    false, false);

% create an implicit function to pass the FEM data as the constraint fucntion to fmincons
%nlconfun = @ (delParams) nonlinear_constraints(delParams, ...
%                                               sensitivity, ...
%                                               channels);

nObjs = numel(blockedSets);
del = 1e-8;
method = 'central';
if strcmpi(method,'forward')  
    thetaNew = nan(sensitivity.designParams.nParams,nObjs,1);
elseif strcmpi(method,'central')
    thetaNew = nan(sensitivity.designParams.nParams,nObjs,2);
else
    error('unrecognized method')
end
FDtheta = nan(sensitivity.designParams.nParams,nObjs);

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
        costfunTimer = tic;
        xdel = xo;
        xdel(i) = xdel(i) + del; 
        thetaNew(i,:) = costfun(xdel);
        %[~,~,thetaNew(i)] = costfun(xdel);
        FDtheta(i,:) = (thetaNew(i,:)-thetaOld')/del;
        toc(costfunTimer);
    elseif (strcmpi(method,'central'))
        for j = 1:2
            costfunTimer = tic;
            xdel = xo;
            xdel(i) = xdel(i) + del*(-1)^j; 
            thetaNew(i,:,j) = costfun(xdel);
            toc(costfunTimer);
        end
        FDtheta(i,:) = 0.5*(thetaNew(i,:,2)-thetaNew(i,:,1))/del;
    end
end



disp('absolute grad diff')
disp(abs(FDtheta-gradTheta))

disp('relative grad diff')
disp(abs(1-FDtheta./gradTheta))


disp('grad theta')
disp(gradTheta)
fclose ('all');
%

      



