 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/25/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%function Main
%clear all;
close all;
format long e;
%clc
% global variables used all over the code
global G_funcEvalCounter
global G_simDirectory
global G_filePrefix
G_simDirectory = './';
G_filePrefix = 'test';
G_funcEvalCounter = 0;

path(path, '../NURBS/nurbs_toolbox')
path(path, '../IGFEM-Curves-2D/M_geom_toolbox')
path(path, '../IGFEM-Curves-2D/M_preFEM')
path(path, '../IGFEM-Curves-2D/M_channels')
path(path, '../IGFEM-Curves-2D/M_FEM')
path(path, '../IGFEM-Curves-2D/mx_FEM')
path(path, '../IGFEM-Curves-2D/M_postprocessing')
path(path, './M_optimization')
path(path, './M_opt_postprocessing')
path(path, './mx_sensitivity')
path(path, '../SISL')
path(path, '../IGFEM-Curves-2D/ChannelFiles')

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
mesh.boundary.xf = 0.1;
mesh.boundary.yi = 0.0;
mesh.boundary.yf = 0.1;
%
% element connectivity and node coordinates
[mesh.elem.elem_node,mesh.node.coords] ...
    = generate_uniform_mesh([1,1],... % number of elements in x- and y-directions
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III) Other user inputs that must be specified regardless of choice I or II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify all channel files and cost function types in the cell arrays
% channelFiles and costFunctionTypes
% cost function types available are
% i) DOMAIN_AVE_TEMP
% ii) DOMAIN_P_NORM_TEMP: the power of the norm costFunction.normp must be
%                         provided
% iii) COMPLIANCE: for testing. the adjoint force nablaPI_1 
%                  should be the same as UUR
% Also specify the numerical vector costFuncPows and maxVolFRacs;
% If an element of costFuncTypes is not DOMAIN_P_NORM
% the value of the corresponding costFuncPows is not relevant
% NOTE:All 4 variables must be arrays of equal length
%channelFile = 'serpentine_simpleB_8_params.channel';
%channelFile = 'parallel2_startA_domain_bounds.channel';
channelFile = 'DK_DP_check.channel';
%costFuncType = 'DOMAIN_P_NORM_TEMP';
costFuncType = 'DOMAIN_P_NORM_TEMP';
costFuncNormp = 1;
maxVolFrac = 0.005; % only relavant when vol constraint is imposed
minPolyAngle = 2; % in degrees
minPolyArea = 0.1^2*(mesh.boundary.xf - mesh.boundary.xi) ...
                   *(mesh.boundary.yf - mesh.boundary.yi) * ones(1,1);

domainThickness = 3e-3; % in meters

% randomize initial values of design parameters subject to bounds
randomizeIniVals = false;                       
% boundary conditions
% 1 : x = xi, 2 : x = xf, 3 : y = yi, 4 : x = yf
mesh.BCs.boundaries = []; % totally insulated boundaries
%mesh.BCs.boundaries = 1;
%mesh.BCs.types = 1;
%mesh.BCs.values_or_funcs = {20};
%mesh.BCs.boundaries = [3,4]; 
%mesh.BCs.types = [2,1];
%mesh.BCs.values_or_funcs = {10,20};
%mesh.BCs.values_or_funcs = {2000,@(x,y) 20*x/0.1};
%mesh.BCs.values_or_funcs = {1,21.5};
%mesh.BCs.values_or_funcs = {u,u,u,u};   

% specify convective heat coefficient
mesh.convect.coef = 0.0; 
mesh.convect.Tref = 0.0;

% for nodes_curves_intersect.m and move_nodes.m 
femParams.tol.node = 1e-6; % tolerance for checking if original nodes coincide with channels
femParams.moveNode.del = 0.5e-3; % distance to move node when it coincides with channel
femParams.moveNode.maxAttemps = 10;
femParams.tol.boundary = 1e-13; % also for set_boundary_conditions.m

% tolerances for edges_curve_intersect.cpp
femParams.tol.nurbsParam = 1.e-6;
femParams.tol.epsco = 1e-15; % computational tolerance for curve_edge_intersect.mexw64
femParams.tol.epsge = 1e-6; % geometrical tolerance for curve_edge_intersect.mexw64
femParams.tol.cosAngleTol = cos(20*pi/180.0); % cosine of the min angle of the resulting element after edge flipping
femParams.tol.halfLineWidthFrac = 1e-4; % half-width fraction of channels
femParams.tol.vert = tan(89.9*pi/180); % tan of max positive angle wrt horizontal axis beyond which a line is considered vertical
femParams.opt.maxRefineLevel = 6;
femParams.opt.refineJuncElem = false; % force refinement of element with branching

% tolerenace for approximating slender child element with triangle or quadrilateral
% slenderTol = []; % if approximation not desired
femParams.slenderTol.minAngle = 5; % min angle of a triangular child element below which NURBS is approximated by line segments
femParams.slenderTol.maxAspectRatio = inf;  % max aspect ratio of a quadrialteral child element below which NURBS is approximated by line segments

femParams.polyIGFEM = true; % use polynomial IGFEM
femParams.supg = false; % apply SUPG

% gauss quadrature schemes
triNpt1D = 2;
triNpt2D = 3;
quaNpt1Dt = 4; % number of gauss points per knot span along channel
quaNpt1Dn = 4; % number of gauss points per knot span normal to channel
quadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM
 
%% NO USER INPUT BELOW THIS LINE
%% Initialize required FEM data that does not change from one simulation to the next
femParams.gauss.line = gauss_points_and_weights(true,triNpt1D,1,'combined');
femParams.gauss.elem = gauss_points_and_weights(true,triNpt2D,2,'combined');
if (femParams.polyIGFEM)
    femParams.gauss.quadElem = gauss_points_and_weights(false,quadNpt1D,2,'combined');
else
    femParams.gauss.qua1Dt = gauss_points_and_weights(false,quaNpt1Dt,1,'separate'); % along channel
    femParams.gauss.qua1Dn = gauss_points_and_weights(false,quaNpt1Dn,1,'separate'); % normal to channel
end


domainVol =  (mesh.boundary.xf-mesh.boundary.xi) ...
            *(mesh.boundary.yf-mesh.boundary.yi)*domainThickness;

[channels,mesh.heatSourceFunc,sensitivity.designParams]...
    = read_channels(channelFile);

channels.domainVol = domainVol;
channels.maxVolFrac = maxVolFrac;
%[channels.polygons,sensitivity.designParams.polyVertices2params] ...
%    = channel_polygons(channels,[]);
channels.sinMinPolyAngle = sin(pi/180*minPolyAngle);
channels.minPolyArea = minPolyArea;

sensitivity.options.TypicalX = 0.1*ones(sensitivity.designParams.nParams,1);
% change sensitivity parameters
switch costFuncType
    case 'DOMAIN_P_NORM_TEMP'
        sensitivity.costFunction.type = 0;
    case 'DOMAIN_AVE_TEMP'
        sensitivity.costFunction.type = 1;
    case 'COMPLIANCE'
        sensitivity.costFunction.type = 2;
    otherwise
        error('unrecognized cost function type')
end
sensitivity.costFunction.normp =  costFuncNormp;
 % for passing into mx_assemble_pseudo_adjoint_forces
sensitivity.costFunction.nDesignParams = sensitivity.designParams.nParams; 
%% Initialize required FEM data that does not change from one function evaluation to the next
% label edges for refinement later. Important!!!!!!!!!!
mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 
% the edge_node information is generated by generate_conforming_mesh
%fprintf('generate edge_node \n')
mesh.edge.edge_node = find_edge_node(mesh.elem.elem_node);

mesh.edge.length = find_edge_length(mesh.edge.edge_node,mesh.node.coords);
mesh.edge.minLength = min(mesh.edge.length);
femParams.tol.halfLineWidth = mesh.edge.minLength*femParams.tol.halfLineWidthFrac;

 %%  anonymous functions 
% create an implicit function to pass the FEM data as the objective fucntion to fmincons
% "@ (param_val)" is the parameter to be optimized, the other parameters are
% only needed to calculate the cost function and are not our primary variables
iniGuess = zeros(sensitivity.designParams.nParams,1); % initial guess for the change in design parameter

costfun = @ (delParams) check_Dmcft(delParams,  ...
                                    mesh, ...
                                    channels, ...
                                    femParams, ...
                                    sensitivity);


del = 1e-7;
method = 'central';
chanNum = 1;
lineSegNum = 1;

xo = zeros(sensitivity.designParams.nParams,1);
fprintf('\nreference mcft and Dmcft \n');
costfunTimer = tic;
[mcft,Dmcft] = costfun(xo);
mcft = mcft{chanNum}(:,lineSegNum);
Dmcft = Dmcft{chanNum}(:,:,lineSegNum);
toc(costfunTimer);

if (strcmpi(method,'forward'))
    costfunTimer = tic;
    xdel = xo;
    xdel = xdel + del; 
    mcftNew = costfun(xdel);
    mcftNew = mcftNew{chanNum}(:,lineSegNum);
    FD_Dmcft = (mcftNew - mcft)/del;
    toc(costfunTimer);
elseif (strcmpi(method,'central'))
    costfunTimer = tic;
    xdel = xo;
    xdel = xdel-del; 
    mcftNew1 = costfun(xdel);
    mcftNew1 = mcftNew1{chanNum}(:,lineSegNum);
    toc(costfunTimer);
    
    costfunTimer = tic;
    xdel = xo;
    xdel = xdel+del; 
    mcftNew2 = costfun(xdel);
    mcftNew2 = mcftNew2{chanNum}(:,lineSegNum);
    toc(costfunTimer);
    
    CD_Dmcft = 0.5*(mcftNew2 - mcftNew1)/del
end

%
read_channels(channelFile);
fclose ('all');


      



