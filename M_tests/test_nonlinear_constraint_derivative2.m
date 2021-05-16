%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/25/2014
%%% Last modified on 9/8/2015
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
clear all
path(path, '../../IGFEM-Curves-2D/M_geom_toolbox')
path(path, '../../IGFEM-Curves-2D/M_postprocessing')
path(path, '../../IGFEM-Curves-2D/M_preFEM')
path(path, '../../IGFEM-Curves-2D/M_channels')
path(path, '../../IGFEM-Curves-2D/M_FEM')
path(path, '../../IGFEM-Curves-2D/mx_FEM')
path(path, '../../IGFEM-Curves-2D/ChannelFiles')
path(path,'../../NURBS/nurbs_toolbox')
path(path, '../M_optimization')
path(path, '../mx_sensitivity')
path(path, '../../SISL/mx_SISL')
path(path, '../M_opt_postprocessing')

global G_simDirectory
global G_filePrefix
global G_out
global G_Tmax2TnRatio
global G_funcEvalCounter
G_simDirectory = './';
G_filePrefix = 'test_P';
G_out = [];
G_Tmax2TnRatio = 1; % initial guess for ratio
G_funcEvalCounter = 0;

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
mesh.boundary.xf = 0.15;
mesh.boundary.yi = 0.0;
mesh.boundary.yf = 0.2;
%
% element connectivity and node coordinates
[mesh.elem.elem_node,mesh.node.coords] ...
    = generate_uniform_mesh([9,12],... % number of elements in x- and y-directions
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
mesh.elem.heatSource = 500.0*ones(size(mesh.elem.elem_node,1),1); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III) Other user inputs that must be specified regardless of choice I or II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%channelFile = 'serpentine_simpleB_8_params.channel';
channelFile = 'parallel2_test_pressure_constraint_w_diams.channel';
%channelFile = 'DK_DP_check_branchingA_2params.channel';

readOptions.boundsFile = [];
readOptions.sampleFile = [];
readOptions.sampleNum = [];
%readOptions.margeTriangles = false;
%readOptions.polygonFile = 'parallelTwo_w_diams_NE.polygon';
readOptions.polygonFile = [];
readOptions.nlconfun = @nonlinear_constraints;
readOptions.nlcon.minPolyArea = 0.001*0.15*0.2;
readOptions.nlcon.sinMinPolyAngle = sin(5*pi/180);
readOptions.nlcon.minSidePolyArea = 0.001*0.15*0.2;
readOptions.nlcon.sinMinSidePolyAngle = sin(1*pi/180);
readOptions.nlcon.minDistSq = (0.001*0.15)^2;
readOptions.nlcon.areaScale = 0.15*0.2;
readOptions.nlcon.distSqScale = 0.15^2;
readOptions.nlcon.type = 'Pmin';
readOptions.nlcon.minP = 15*1000;
readOptions.nlcon.maxP = 20*1000;
readOptions.nlcon.maxTmax = 60;
readOptions.nlcon.minAorVFrac = 0.018;
readOptions.nlcon.maxAorVFrac = 0.2;

% normalized normal constraints for generating Pareto front
readOptions.nlcon.NNC.objInds = [1,2,3];
readOptions.nlcon.NNC.anchorPts = [25.49,246.5,70;
                                   46994,261,2000;
                                   0.1,0.05,0.02];
[readOptions.nlcon.NNC.normalization,...
 readOptions.nlcon.NNC.scaledAnchorPts,...
 readOptions.nlcon.NNC.UtopiaLineVec,...
 hyperPlanePts]= NNC_parameters(readOptions.nlcon.NNC.anchorPts,51);
readOptions.nlcon.NNC.hyperPlanePt = hyperPlanePts(:,26);

domainThickness = 3e-3; % in meters
minPolyAngle = 2; % in degrees

readOptions.nlcon.area = (mesh.boundary.xf - mesh.boundary.xi) ...
                        *(mesh.boundary.yf - mesh.boundary.yi);
domainVol = domainThickness*readOptions.nlcon.area;
minPolyArea = 0.01*readOptions.nlcon.area;

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
triNpt2D = 25;
senstvtyTriNpt1D = 4; % number of gauss points for line quadrature in sensitivity analysis
senstvtyTriNpt2D = 25; % number of gauss points for triangular quadrature in sensitivity analysis
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


[channels,mesh.heatSourceFunc,...
 sensitivity.designParams,sensitivity.restrictedParams]...
    = preprocess_channels(channelFile,readOptions);
channels.domainVol = domainVol;
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

del = 1e-7;
method = 'central';

sensitivity.costFunction.type = 'AREA';
func = @(x) nonlinear_constraints(x,...
                                  sensitivity.designParams,...
                                  sensitivity.restrictedParams,...
                                  channels,readOptions.nlcon,...
                                  mesh,femParams,sensitivity);
nParams = sensitivity.designParams.nParams ...
            +sensitivity.restrictedParams.nParams;
xo = zeros(nParams,1);
[gOld,~,gradg] = func(xo);
FD_g = nan(nParams,numel(gOld));
for i = 1:nParams
    fprintf('\n-----------------------------------------------------------\n')
    fprintf('design param %i/%i',i,sensitivity.designParams.nParams);
    fprintf('\n-----------------------------------------------------------\n')
    if (strcmpi(method,'forward'))
        xdel = xo;
        xdel(i) = xdel(i) + del; 
        gNew1 = func(xdel);
        FD_g(i,:) = (gNew1 - gOld)'/del;
        
    elseif (strcmpi(method,'central'))
        xdel = xo;
        xdel(i) = xdel(i) + del; 
        gNew2 = func(xdel);

        xdel = xo;
        xdel(i) = xdel(i) - del; 
        gNew1 = func(xdel);

        FD_g(i,:) = 0.5*(gNew2 - gNew1)'/del;

    else
        error('unrecognized method')
    end
end
abs_gradg_diff = abs(FD_g - gradg);
fprintf('maximum abs diff = %g \n',max(abs_gradg_diff(:)));

tol = 1e-10;
nzInd = abs(gradg) > tol;
rel_diff = abs(1 - FD_g(nzInd)./gradg(nzInd));
fprintf('maximum rel diff = %g \n',max(rel_diff(:)));
