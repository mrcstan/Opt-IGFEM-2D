 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 9/23/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   channelFiles
%   polygonFile,...
%   directory
%   filePrefix
%   inputFile: 
%   costFuncType:
%       i) P_NORM: costFuncNormp must be provided
%       ii) VARIANCE
%       iii) PRESSURE
%   costFuncNormp
%       If costFuncType ~= P_NORM, then this is the output p-norm
%       temperature
%   anchorPts
%   nHyperPlanePts1: number of points along the axis of the first objective
%                    in the hyperplane of the space of the objective
%                    functions
%   randomizeFIrst: randomize the configuration of the first channel
%   mergeTriangles: merge the triangles of the channel network
function [] = generate_pareto_front(channelFile,...
                                    polygonFile,...
                                    directory,...
                                    filePrefix,...
                                    inputFile,...
                                    costFuncType,...
                                    costFuncNormp,...
                                    anchorPts,...
                                    nHyperPlanePts1,...
                                    nlconType, ...
                                    nlconMinP, ...
                                    nlconMaxP, ...
                                    nlconMaxTmax, ...
                                    nlconMinAorVFrac, ...
                                    nlconMaxAorVFrac, ...
                                    mergeTriangles)
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
global G_Tmax2TnRatio % for max temperature constraint

%profile on;
% matlabpool close
% matlabpool(6)

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
path(path, '../IGFEM-Curves-2D/SampleFiles')
path(path, '../IGFEM-Curves-2D/GeomConstrFiles')
path(path, './InputFiles')
%% MESH AND USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information for generating mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
senstvtyTriNpt2D = 16; % number of gauss points for triangular quadrature in sensitivity analysis
senstvtyQuadNpt1D = 4; % number of gauss points in one direction for quad child element in polynomial IGFEM sensitivity analysis 
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
        'GradConstr', 'on', ...    % supply gradient for the constraint funcitons
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
% scaleObj = true;
scaleConstraint = true;
%% NO USER INPUT BELOW THIS LINE
%% Initialize required FEM data that does not change from one simulation to the next
sensitivity.gauss.line = gauss_points_and_weights(true,senstvtyTriNpt1D,1,'combined');
sensitivity.gauss.elem = gauss_points_and_weights(true,senstvtyTriNpt2D,2,'combined');
if femParams.otherFlags.polyIGFEM
    sensitivity.gauss.quadElem = gauss_points_and_weights(false,senstvtyQuadNpt1D,2,'combined');
end
% label edges for refinement later. Important!!!!!!!!!!
if femParams.refine.maxRefineLevel
    mesh.elem.elem_node = label(mesh.node.coords,mesh.elem.elem_node); 
end
femParams.tol.halfLineWidth = mesh.edge.minLength*femParams.tol.halfLineWidthFrac;
femParams.moveNode.dist = mesh.edge.minLength*femParams.moveNode.distFrac;


% check optimization inputs
splitCostFuncType = regexp(costFuncType,',','split');
nObjs = numel(splitCostFuncType);
if (nObjs < 2)
    error('multi-objective optimization requires two or more objectives')  
end
if (nObjs ~= size(anchorPts,1) || nObjs ~= size(anchorPts,2))
    error('number of anchor points and their coordinates must equal number of objectives')
end

sensitivity.costFunction.objOpt.node = []; % outlet node not requested
% intDomainType indicates whether the primary objective function is
% integrated over the whole domain. Note that this does not apply to the
% second possible objective function, i.e, the 1-norm.
% 0 = whole domain, 1 = channels only
sensitivity.costFunction.objOpt.intDomainType = 0;
switch splitCostFuncType{end}
    case 'T'
        sensitivity.costFunction.objOpt.normp = costFuncNormp;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = 'P_NORM';
    case 'P'
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = 'PRESSURE';
        sensitivity.costFunction.area = mesh.domainArea;
        sensitivity.costFunction.objOpt.intDomainType = 1;  
    case 'A'
        sensitivity.costFunction.type = 'AREA';
        sensitivity.costFunction.area = mesh.domainArea;
    otherwise
        error('unrecognized cost function type')
end

% Assume that there are only three possible objectives
readOptions.nlcon.NNC.objInds = zeros(3,1);
for i = 1:nObjs
    if strcmpi(splitCostFuncType{i},'T')
        readOptions.nlcon.NNC.objInds(1) = i;
    elseif strcmpi(splitCostFuncType{i},'P')
        readOptions.nlcon.NNC.objInds(2) = i;
    elseif strcmpi(splitCostFuncType{i},'A')
        readOptions.nlcon.NNC.objInds(3) = i;
    else
        error('some of the objectives are unknown')
    end
end

readOptions.nlcon.NNC.anchorPts = anchorPts;
[readOptions.nlcon.NNC.normalization,...
 readOptions.nlcon.NNC.scaledAnchorPts,...
 readOptions.nlcon.NNC.UtopiaLineVec,...
 hyperPlanePts]= NNC_parameters(anchorPts,nHyperPlanePts1);
nSimulations = size(hyperPlanePts,2);

if (ischar(nlconType))
    splitStr = regexp(nlconType,',','split');
    for i = 1:numel(splitStr)       
        if (~strcmpi(splitStr{i},'Pmin') ...
            && ~strcmpi(splitStr{i},'Pmax') ...
            && ~strcmpi(splitStr{i},'Tmax') ...
            && ~strcmpi(splitStr{i},'Amin') ...
            && ~strcmpi(splitStr{i},'Amax') ...
            && ~strcmpi(splitStr{i},'Vmin') ...
            && ~strcmpi(splitStr{i},'Vmax'))
            error('unknown nonlinear constraint type')
        end
    end
elseif (~isempty(nlconType))
    error('unknown nonlinear constraint type')
end

readOptions.nlcon.type = nlconType;
readOptions.nlcon.minP = nlconMinP;
readOptions.nlcon.maxP = nlconMaxP;
readOptions.nlcon.maxTmax = nlconMaxTmax;
readOptions.nlcon.area = mesh.domainArea;
readOptions.nlcon.minAorVFrac = nlconMinAorVFrac;
readOptions.nlcon.maxAorVFrac = nlconMaxAorVFrac;

% Initialize the global variables
directory = ['./',directory];
mkdir(directory);
G_filePrefix = filePrefix; % prefix to be added to the file names

% write simulation parameters to file
fnameSim = [directory,'/',filePrefix,'_simulation_parameters.txt'];
fidSim = fopen(fnameSim,'w');
fprintf(fidSim,'number of simulations %i\n',nSimulations);
fprintf(fidSim,'channel file %s\n',channelFile);
fprintf(fidSim,'polygon file %s\n',polygonFile);
fprintf(fidSim,'cost function type %s\n',costFuncType);
fprintf(fidSim,'cost function norm %g\n',costFuncNormp);
fprintf(fidSim,'anchor points:\n');
for i = 1:size(anchorPts,1)
    fprintf(fidSim,'%g\t',anchorPts(i,:));
    fprintf(fidSim,'\n');
end
fprintf(fidSim,'number of hyperplane points in first obj direction %g\n',nHyperPlanePts1);
fprintf(fidSim,'algorithm %s\n',sensitivity.options.Algorithm);
fprintf(fidSim,'TolFun %g\n',sensitivity.options.TolFun);
fprintf(fidSim,'TolX %g\n',sensitivity.options.TolX);
fprintf(fidSim,'TolCon %g\n',sensitivity.options.TolCon);
fprintf(fidSim,'normalized normal constraint is applied \n');
fprintf(fidSim,'area constraint scaled %i\n',scaleConstraint);
if (isempty(nlconType))
    fprintf(fidSim,'nonlinear constraint type not applied\n');
else
    fprintf(fidSim,'nonlinear constraint type %s\n',nlconType);
end
if (isempty(nlconMinP))
    fprintf(fidSim,'minimum pressure (Pa) not applied\n');
else
    fprintf(fidSim,'minimum pressure (Pa) %g\n',nlconMinP);
end
fprintf(fidSim,'maximum pressure (Pa) %g\n',nlconMaxP);
fprintf(fidSim,'maximum Tmax (Celsius) %g\n',nlconMaxTmax);
fprintf(fidSim,'minimum area or volume fraction %g\n',nlconMinAorVFrac);
fprintf(fidSim,'maximum area of volume fraction %g\n',nlconMaxAorVFrac);
fprintf(fidSim,'side triangle min angle %g\n',minSidePolyAngle);
fprintf(fidSim,'side triangle min area %g\n',minSidePolyArea);
fprintf(fidSim,'interior triangle min angle %g\n',minPolyAngle);
fprintf(fidSim,'interior triangle min area %g\n',minPolyArea);
fclose(fidSim);

% copy input file
spltStr = regexp(inputFile,'/','split');
copyfile(['./InputFiles/',inputFile],[directory,'/',spltStr{end}]);

%read channel file and design parameters
readOptions.figBounds = [mesh.boundary.xi, mesh.boundary.xf, ...
                         mesh.boundary.yi, mesh.boundary.yf];
readOptions.filePrefix = [directory,'/',filePrefix];
readOptions.polygonFile = polygonFile;
readOptions.polygonFig = [readOptions.filePrefix,'_polygons'];
readOptions.nlconfun = @nonlinear_constraints;
readOptions.nlcon.sinMinPolyAngle = sin(pi/180*minPolyAngle);
readOptions.nlcon.minPolyArea = minPolyArea;
readOptions.nlcon.sinMinSidePolyAngle = sin(minSidePolyAngle*pi/180);
readOptions.nlcon.minSidePolyArea = minSidePolyArea;

readOptions.mergeTriangles = mergeTriangles;

readOptions.boundsFile = [];
readOptions.sampleFile = [];
readOptions.sampleNum = [];
readOptions.Toffset = mesh.convect.Toffset;

[channels,...
 sensitivity.designParams, ...
 sensitivity.restrictedParams]...
    = preprocess_channels(channelFile,readOptions); 
channels.domainVol = mesh.domainVol;
channels.domainArea = mesh.domainArea;
readOptions.nlcon.PScale = nlconMaxP;


% Pareto front candidates
paretoCands = nan(nObjs,nObjs+nSimulations);
paretoCands(:,(end-nObjs+1):end) = anchorPts;
paretoCandsTmax = nan(1,nSimulations);
paretoFigName = [directory,'/',filePrefix,'_pareto_candidates'];
for simNum = 1:nSimulations
    readOptions.nlcon.NNC.hyperPlanePt = hyperPlanePts(:,simNum);
    %Tu = readOptions.nlcon.NNC.hyperPlanePt(1)*readOptions.nlcon.NNC.normalization(1) ...
    %       + readOptions.nlcon.NNC.anchorPts(1,1)
    %Pu = readOptions.nlcon.NNC.hyperPlanePt(2)*readOptions.nlcon.NNC.normalization(2) ...
    %       + readOptions.nlcon.NNC.anchorPts(2,2)
    G_simDirectory = [directory,'/sim',num2str(simNum),'/'];
    mkdir(G_simDirectory)
    G_history.iter = [];
    G_history.fval = [];
    G_history.SD = [];
    G_history.x = [];
    G_history.cstr = [];
    G_history.inPressure = [];
    G_history.inFlow = [];
    G_history.Tave = [];
    G_history.Tpnorm = [];
    G_history.Tmax = [];
    G_history.firstorderopt = [];
    G_history.Tmax2TnRatio = [];
    G_history.volFrac = [];
    G_history.areaFrac = [];
    G_history.Tin = [];
    G_history.Tout = [];
    G_history.TaveChan = [];
    G_history.TminChan = [];
    G_history.TmaxChan = [];
    
    G_funcEvalCounter =0;
    G_constraint = [];
    G_constraint.g = []; % constraint function values
    G_constraint.history.nodalT = []; % Only applicable when the nodal temperature constraints are evaluated at different flow conditions
    
    G_out = [];
    G_out.costFuncType = costFuncType;
    G_out.gauss = femParams.gauss;
    G_out.designParams = sensitivity.designParams; % for output only
    G_out.normp = costFuncNormp;  
    G_Tmax2TnRatio = 1.2; % initial guess for ratio
    
    % IMPORTANT: must ensure that the normalization factor is positive !
    sensitivity.objScale = readOptions.nlcon.NNC.normalization(end);
    sensitivity.objOffset = anchorPts(end,end);        
    G_out.objScale = sensitivity.objScale;
    G_out.objOffset = anchorPts(end,end);
    try 
        fprintf('\n---------------------------------------------------------\n')
        fprintf('Starting simulation %i ',simNum)
        fprintf('\n---------------------------------------------------------\n')
        write_channel_file([readOptions.filePrefix,'_ini.channel'],...
                            channels,...
                            mesh.convect.Toffset, ...
                            sensitivity.designParams,...
                            'w',[])

        sensitivity.costFunction.objOpt.nDesignParams = sensitivity.designParams.nParams; 

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
                                                   sensitivity, ...
                                                   1,true,true);

        if (scaleConstraint)
            readOptions.nlcon.areaScale = mesh.domainArea;
            %readOptions.nlcon.distSqScale = lengthScale^2;
        end
        nlconfun = @ (delParams) nonlinear_constraints(delParams, ...
                                               sensitivity.designParams, ...
                                               sensitivity.restrictedParams, ...
                                               channels,...
                                               readOptions.nlcon, ...
                                               mesh, ...
                                               femParams, ...
                                               sensitivity);

        
        fminconTimer = tic;
        [paramVals, fval, exitflag, output, lambda, grad] ...
            = fmincon (costfun,...    % Pass the objective function
                       iniGuess,...    % optimization parameter initial value (x0) (multiplied by velocity). x0 can be a scalar, vector, or matrix.
                       [], [],...      % the inequalities: subject to the linear inequalities A*x <= b. Here x is param_val.
                       Aeq, beq,...      % the linear equalities Aeq*x = beq. Here x is param_val.
                       lb, ub,...      % bound for the param_val
                       nlconfun,...   % the nonlinear constriant: the nonlinear inequalities c(x) or equalities ceq(x) defined in nonlcon. 
                       sensitivity.options);       % "sensitivity.options" which was defined earlier
        fmincontime = toc(fminconTimer);
        fprintf('fmincon time = %g \n',fmincontime)
        
        % gather objectives into paretoCand
        if readOptions.nlcon.NNC.objInds(1)
            % T
            paretoCands(readOptions.nlcon.NNC.objInds(1),simNum) = G_history.Tpnorm(end);
            paretoCandsTmax(simNum) = G_history.Tmax(end);
        end
        if readOptions.nlcon.NNC.objInds(2)
            % P
            paretoCands(readOptions.nlcon.NNC.objInds(2),simNum) = G_history.inPressure(end);
        end
        if readOptions.nlcon.NNC.objInds(3)
            % A
            paretoCands(readOptions.nlcon.NNC.objInds(3),simNum) = G_history.areaFrac(end);
        end
        if nnz(readOptions.nlcon.NNC.objInds) == 2
            fig = figure('visible','off');
            plot(paretoCands(1,:),paretoCands(2,:),'ro',...
                 'markerfacecolor','r','markersize',12)
            xlabel('obj_1','fontsize',20)
            ylabel('obj_2','fontsize',20)
            set(gca,'fontsize',20)
            saveas(fig, paretoFigName, 'jpg');
        elseif nnz(readOptions.nlcon.NNC.objInds) == 3
            fig = figure('visible','off');
            plot3(paretoCands(1,:),paretoCands(2,:),paretoCands(3,:),...
                  'ro','markerfacecolor','r','markersize',12)
            xlabel('obj_1','fontsize',20)
            ylabel('obj_2','fontsize',20)
            zlabel('obj_3','fontsize',20)
            set(gca,'fontsize',20)
            saveas(fig, paretoFigName, 'jpg');
        end
        % update channel
        channels = update_channels(paramVals, ...
                                   sensitivity.designParams, ...
                                   sensitivity.restrictedParams, ...
                                   channels, ...
                                   'add', ...
                                   false);
        
        %
        filePrefix = [G_simDirectory,G_filePrefix];
        save([filePrefix,'_allVariables'])
        simStatus = 'succeed';
        fidSim = fopen(fnameSim,'a');
        fprintf(fidSim,'\nSimulation %i %s, fmincon time = %g \n',simNum,simStatus, fmincontime);
    
    catch err 
        filePrefix = [G_simDirectory,G_filePrefix];
        save([filePrefix,'_allVariables'])
        fprintf('SIMULATION %i FAILED because \n',simNum)
        warning(getReport(err))
        simStatus = 'failed';
        fidSim = fopen(fnameSim,'a');
        fprintf(fidSim,'\nSimulation %i %s \n',simNum,simStatus);
        return
    end
    
    fclose ('all');
end % loop through all simulations

      
end


