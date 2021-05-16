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

% global variables used all over the code for outputs only
% VERY IMPORTANT TO REINITIALIZE THESE GLOBAL VARIABLES AFTER EVERY
% SIMULATION !!!!!!!!!!!!!!
global G_history
global G_filePrefix
global G_simDirectory
global G_funcEvalCounter
global G_constraint
global G_out % for output every iteration, avoid using it for calculation
G_filePrefix = '30x40_mesh';
G_simDirectory = 'parallel2_movable_mid_channel_xo_p04_yo_p1_x1_p11_y1_p1_ror1_p015_QoQ1_250_do_p185_p015/';
mkdir(G_simDirectory)

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
path(path, '../SISL/mx_SISL')
path(path, '../IGFEM-Curves-2D/ChannelFiles')
path(path, '../MatlabUsefulFunctions/export_fig')

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
domainThickness = 0.003;
%
% element connectivity and node coordinates
[mesh.elem.elem_node,mesh.node.coords] ...
    = generate_uniform_mesh([30,40],... % number of elements in x- and y-directions
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
mesh.elem.heatSource = 0.0*ones(size(mesh.elem.elem_node,1),1); 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III) Other user inputs that must be specified regardless of choice I or II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelFile = 'parallel2_movable_mid_channel.channel';
svFile = 'parallel2_movable_mid_channel_xo_p04_yo_p1_x1_p11_y1_p1_ror1_p015_QoQ1_250_do_p185_p015_30x40_mesh';
costFuncType = 'DOMAIN_P_NORM_TEMP';
costFuncNormp = 8;

minPolyAngle = 5; % in degrees
domainArea = (mesh.boundary.xf - mesh.boundary.xi) ...
                   *(mesh.boundary.yf - mesh.boundary.yi);
domainVol =  domainArea*domainThickness;

minPolyArea = 0.001 * domainArea;
minSidePolyAngle = 0.5; % in degrees
minSidePolyArea = 0.001* domainArea;
domainThickness = 3e-3; % in meters
                    
% boundary conditions
% 1 : x = xi, 2 : x = xf, 3 : y = yi, 4 : x = yf
mesh.BCs.boundaries = []; % totally insulated boundaries
%mesh.BCs.boundaries = 1;
%mesh.BCs.types = 1;
%mesh.BCs.values_or_funcs = {27};
%mesh.BCs.boundaries = [3,4]; 
%mesh.BCs.types = [2,1];
%mesh.BCs.values_or_funcs = {10,20};
%mesh.BCs.values_or_funcs = {2000,@(x,y) 20*x/0.1};
%mesh.BCs.values_or_funcs = {1,21.5};
%mesh.BCs.values_or_funcs = {u,u,u,u};   
%mesh.BCs.boundaries = [3,4]; 
%mesh.BCs.types = [2,1];
%mesh.BCs.values_or_funcs = {30,27};

% specify convective heat coefficient
mesh.convect.coef = 0.0; 
mesh.convect.Tref = 0.0;

% for nodes_curves_intersect.m and move_nodes.m 
femParams.tol.node = 1e-6; % tolerance for checking if original nodes coincide with channels
femParams.moveNode.distFrac = 0.05; % distance to move node when it coincides with channel
                                    % or when a branching point or kink coincides
                                    % with an element edge, as a fraction of the
                                    % minimum edge length
femParams.moveNode.maxAttemps = 10;
femParams.tol.boundary = 1e-13; % also for set_boundary_conditions.m

% tolerances for edges_curve_intersect.cpp
femParams.tol.nurbsParam = 1.e-8;
femParams.tol.epsco = 1e-15; % computational tolerance for curve_edge_intersect.mexw64
femParams.tol.epsge = 1e-6; % geometrical tolerance for curve_edge_intersect.mexw64
femParams.tol.cosAngleTol = cos(20*pi/180.0); % cosine of the min angle of the resulting element after edge flipping
femParams.tol.halfLineWidthFrac = 1e-4; % half-width fraction of channels
femParams.tol.vert = tan(89.9*pi/180); % tan of max positive angle wrt horizontal axis beyond which a line is considered vertical
femParams.opt.maxRefineLevel = 8;
femParams.opt.refineJuncElem = false; % force refinement of element with branching
femParams.tol.intersectEdges = 1e-13; % for intersect_edges in single_edge_curves_intersect.m

femParams.tol.channelSelfIntersect = 1e-13; % tolerance for channels_self_intersections

% tolerenace for approximating slender child element with triangle or quadrilateral
% slenderTol = []; % if approximation not desired
femParams.slenderTol.minAngle = 5; % min angle of a triangular child element below which NURBS is approximated by line segments
femParams.slenderTol.maxAspectRatio = inf;  % max aspect ratio of a quadrialteral child element below which NURBS is approximated by line segments

femParams.polyIGFEM = true; % use polynomial IGFEM
femParams.supg = true; % apply SUPG

% gauss quadrature schemes
triNpt1D = 4;
triNpt2D = 16;
senstvtyTriNpt1D = 4; % number of gauss points for line quadrature in sensitivity analysis
senstvtyTriNpt2D = 16; % number of gauss points for triangular quadrature in sensitivity analysis
quaNpt1Dt = 4; % number of gauss points per knot span along channel
quaNpt1Dn = 4; % number of gauss points per knot span normal to channel
quadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM
senstvtyQuadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM sensitivity analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IV) SENSITIVITY ANALYSIS INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sensitivity.optimization_type = 'SHAPE';
% 'GradObj', 'on': adjoint method; 'GradObj','off':finite difference
%'DerivativeCheck', 'on', ... %Compare user-supplied derivatives (gradients of objective or constraints) to finite-differencing derivatives. 
% 'GradCon'...    % turn off the gradient mode for the cosntriant funcitons
% different algorithms fmincon uses for solving the optimization problem
%'trust-region-reflective' 
%'interior-point'
%'sqp'
%'active-set'
% define the "option" parameter to later pass to the fmincon function
sensitivity.options ...
    = optimoptions(@fmincon, ...
        'FunValCheck', 'off',... % check for complex, inf or nan obj fun val and terminate if true   
        'GradObj', 'on',...     % supply gradient for the objective funcitons
        'GradCon', 'on', ...    % supply gradient for the constraint funcitons
        'DerivativeCheck', 'off', ... %Compare user-supplied derivatives (gradients of objective or constraints) to finite-differencing derivatives. 
        'FinDiffType','central', ...
        'outputfcn', @outfun,...  % output the history of the converged solution
        'Algorithm', 'sqp',...  % Choose the algorithm to solve the optimization
        'TolCon', 1e-6,...   % tolerance for constraint
        'TolX',1e-6,...      % tolerance for parameter value (multiplied by velocity) 
        'TolFun', 1e-6, ... % tolerence on the funciton we optimize
        'InitBarrierParam',0.01, ...% valid only for interior point algorithm 
        'Display', 'iter-detailed');   
% 'ScaleProblem','obj-and-constr', ... % normalize all constraints and the objective function
% 'DiffMaxChange', 1e-6, ... % maximum change in variables for FD gradients
scaleObj = true;
scaleConstraint = true;
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
% label edges for refinement later. Important!!!!!!!!!!
mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 
% the edge_node information is generated by generate_conforming_mesh
%fprintf('generate edge_node \n')
mesh.edge.edge_node = find_edge_node(mesh.elem.elem_node);

mesh.edge.length = find_edge_length(mesh.edge.edge_node,mesh.node.coords);
mesh.edge.minLength = min(mesh.edge.length);
femParams.tol.halfLineWidth = mesh.edge.minLength*femParams.tol.halfLineWidthFrac;
femParams.moveNode.dist = mesh.edge.minLength*femParams.moveNode.distFrac;


% check optimization inputs
switch costFuncType
    case 'DOMAIN_P_NORM_TEMP'
        sensitivity.costFunction.type = 0;        
        sensitivity.costFunction.normp = costFuncNormp;
        sensitivity.costFunction.isVariance = false;
    case 'DOMAIN_AVE_TEMP'
        sensitivity.costFunction.type = 1; 
        sensitivity.costFunction.isVariance = false;
    case 'COMPLIANCE'
        sensitivity.costFunction.type = 2;
        sensitivity.costFunction.isVariance = false;
    case 'VARIANCE'
        sensitivity.costFunction.type = 0;
        sensitivity.costFunction.normp = 2;
        sensitivity.costFunction.isVariance = true;  
        sensitivity.costFunction.area = domainArea;
    otherwise
        error('unrecognized cost function type')
end
%% global variables
G_history.iter = [];
G_history.fval = [];
G_history.SD = [];
G_history.x = [];
G_history.cstr = [];
G_history.inPressure = [];
G_history.Tave = [];
G_history.Tmax = [];
G_funcEvalCounter =0;
G_constraint = [];
G_history.firstorderopt = [];
G_history.Tmax2TnRatio = [];
G_history.volFrac = [];

G_out = [];
if (strcmpi(costFuncType,'VARIANCE'))
    G_out.isVariance = true;
else
    G_out.isVariance = false;
end 

sensitivity.objScale = 1; % important to reset scale !

 %%  anonymous functions 
[channels,mesh.heatSourceFunc,...
     sensitivity.designParams,sensitivity.restrictedParams]...
        = preprocess_channels(channelFile); 
channels.domainVol = domainVol;

G_out.designParams = sensitivity.designParams; % for output only

%% INITIAL VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sensitivity.designParams.iniVals = [0.185;0.015]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sensitivity.costFunction.nDesignParams = sensitivity.designParams.nParams; 

[lb,ub] = params_bounds(sensitivity.designParams, ...
                        sensitivity.restrictedParams);

% linear equality constraints
[Aeq,beq] = linear_equality_constraints(sensitivity.restrictedParams, ...
                                        sensitivity.designParams.nParams);

% initial guess for the change in design parameter
iniGuess = zeros(sensitivity.designParams.nParams + sensitivity.restrictedParams.nParams,1); 

costfun = @ (delParams) fem_n_sensitivity (delParams,  ...
                                           mesh, ...
                                           channels, ...
                                           femParams, ...
                                           sensitivity);
                                       
if (scaleObj) % must redeclare anonymous function after any of its arguments is changed
    sensitivity.objScale = costfun(iniGuess);
    costfun = @ (delParams) fem_n_sensitivity (delParams,  ...
                                               mesh, ...
                                               channels, ...
                                               femParams, ...
                                               sensitivity);
    if (strcmpi(costFuncType,'DOMAIN_P_NORM_TEMP'))
        G_out.objScale = sensitivity.objScale/(domainArea)^(1/costFuncNormp);
    else
        G_out.objScale = sensitivity.objScale;
    end
else
    if (strcmpi(costFuncType,'DOMAIN_P_NORM_TEMP'))
        G_out.objScale = 1.0./(domainArea)^(1/costFuncNormp);
    else
        G_out.objScale = 1.0;
    end
end
    
nlconfun = [];                                     
fminconTimer = tic;
[param_val, fval, exitflag, output, lambda, grad] ...
    = fmincon (costfun,...    % Pass the objective function
               iniGuess,...    % optimization parameter initial value (x0) (multiplied by velocity). x0 can be a scalar, vector, or matrix.
               [], [],...      % the inequalities: subject to the linear inequalities A*x <= b. Here x is param_val.
               Aeq, beq,...      % the linear equalities Aeq*x = beq. Here x is param_val.
               lb, ub,...      % bound for the param_val
               nlconfun,...   % the nonlinear constriant: the nonlinear inequalities c(x) or equalities ceq(x) defined in nonlcon. 
               sensitivity.options);       % "sensitivity.options" which was defined earlier
fmincontime = toc(fminconTimer);
fprintf('fmincon time = %g \n',fmincontime)  

designParamVals = bsxfun(@plus,G_history.x,sensitivity.designParams.iniVals);

save(svFile)

%{
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
contourf(designParamPt1,designParamPt2,objVals,12)
hold on 
xlabel('d_1 (m)','fontsize',30)
ylabel('d_2 (m)','fontsize',30)
set(gca,'fontsize',30,'color',[0.8,0.8,0.8])
%set(gca,'fontsize',30)
axis image
hcolor = colorbar;
set(hcolor,'YTick',linspace(60,220,9))
ylabel(hcolor,'T_8 (^{\circ}C)','rot',0,'fontsize',30)
set(hcolor,'fontsize',30)
plot(designParamVals(1,:),designParamVals(2,:),'d--','color',...
                        'w','linewidth',3,...
                        'markersize',8,'markerfacecolor','w')
%{
for i = 1:size(designParamVals,2)
    text(designParamVals(1,i),designParamVals(2,i),num2str(i),...
        'fontsize',30,'color','w')
end
%}
                        
%text(designParamVals(1,1),designParamVals(2,1),['d_o=(',num2str(designParamVals(1,1)),',',num2str(designParamVals(1,1)),')'],...
%     'fontsize',30','color','k')
xlim([0.013,0.187])
ylim([0.013,0.187])
%}
fclose ('all');





      



