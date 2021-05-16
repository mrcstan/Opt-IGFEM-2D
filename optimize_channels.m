 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 1/17/2015
%%% Last modified date: 10/19/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   channelFiles:
%   polygonFile:
%   boundsFile_sampleFile:
%        boundsFile (*.channels)/sampleFile
%        ('*.lhs','*.rand','*.rand1','*.rand2','*.smmp')
%   jobNum_sampleNum:
%       job num (if > 0 && nSimulations > 0)/sample number (if < 0 && nSimulations == 1)
%   directory:
%   filePrefix:
%   inputFile: 
%   costFuncType:
%       i) P_NORM: costFuncNormp must be provided
%       ii) VARIANCE
%       iii) PRESSURE
%       iv) AREA
%       v) TPA: weighted combination of p-norm of temperature, pressure and
%               area
%       vi) NODAL_T_IN, NODAL_T_OUT or NODAL_T_INOUT
%   costFuncNormp:
%   costFuncWeights:
%   costFuncScale:
%   nlconType:
%       i) Pmax
%       ii) Pmin
%       iii) Tmax
%       iv) Amax
%       v) Amin
%       vi) Vmax
%       vii) Vmin
%       viii) nodalT
%   nlconMinP:
%   nlconMaxP:
%   nlconMaxTmax:
%   nlconMinAorVFrac:
%   nlconMaxAorVFrac:
%   nlconNodalTbounds: [lower Tin, upper Tin, lower Tout, upper Tout]
%   nSimulations:
%   randomizeFIrst: randomize the configuration of the first channel
%   mergeTriangles: merge the triangles of the channel network
function optimize_channels(channelFile,...
                            polygonFile,...
                            boundsFile_sampleFile,...
                            jobNum_sampleNum,...
                            directory,...
                            filePrefix,...
                            inputFile,...
                            costFuncType,...
                            costFuncNormp,...
                            costFuncWeights,...
                            costFuncScales,...
                            nlconType, ...
                            nlconMinP, ...
                            nlconMaxP, ...
                            nlconMaxTmax, ...
                            nlconMinAorVFrac, ...
                            nlconMaxAorVFrac, ...
                            nlconNodalTbounds, ...
                            nlconMassin, ...
                            nlconInBCs, ...
                            nSimulations, ...
                            randomizeFirst, ...
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

path(path, '../nurbs_toolbox')
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
minPolyAngle = 10; % in degrees
minPolyArea = 0.001 * mesh.domainArea;
minSidePolyAngle = 1; % in degrees
minSidePolyArea = 0.001* mesh.domainArea;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENSITIVITY ANALYSIS INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
senstvtyTriNpt1D = 5; % number of gauss points for line quadrature in sensitivity analysis
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
scaleObj = true;
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
        error(['Need to update objective_n_gradient.m to use channel Tave.',...
               'See PRESSURE in the script for more info'])
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

if (ischar(nlconType))
    splitStr = regexp(nlconType,',','split');
    for i = 1:numel(splitStr)       
        if (~strcmp(splitStr{i},'Pmin') ...
            && ~strcmp(splitStr{i},'Pmax') ...
            && ~strcmp(splitStr{i},'Tmax') ...
            && ~strcmp(splitStr{i},'Amin') ...
            && ~strcmp(splitStr{i},'Amax') ...
            && ~strcmp(splitStr{i},'Vmin') ...
            && ~strcmp(splitStr{i},'Vmax') ...
            && ~strcmp(splitStr{i},'nodalT') )
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
readOptions.nlcon.nodalTbounds = nlconNodalTbounds + mesh.convect.Toffset;
readOptions.nlcon.massin = nlconMassin;
readOptions.nlcon.inBCs = nlconInBCs;
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
fprintf(fidSim,'sample file %s\n',boundsFile_sampleFile);
fprintf(fidSim,'cost function type %s\n',costFuncType);
fprintf(fidSim,'cost function norm %g\n',costFuncNormp);
fprintf(fidSim,'cost function weights %g %g %g\n',costFuncWeights);
fprintf(fidSim,'cost function scales %g %g %g\n',costFuncScales);
fprintf(fidSim,'algorithm %s\n',sensitivity.options.Algorithm);
fprintf(fidSim,'TolFun %g\n',sensitivity.options.TolFun);
fprintf(fidSim,'TolX %g\n',sensitivity.options.TolX);
fprintf(fidSim,'TolCon %g\n',sensitivity.options.TolCon);
fprintf(fidSim,'obj scaled %i\n',scaleObj);
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

for simNum = 1:nSimulations
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
    G_constraint.g = []; % constraint function values
    G_constraint.history.nodalT = []; % Only applicable when the nodal temperature constraints are evaluated at different flow conditions
    G_out = [];
    G_out.costFuncType = costFuncType;
    G_out.gauss = femParams.gauss;
    G_out.normp = costFuncNormp;
    G_out.nuFlag = 1; % evaluate viscosity at average network temperature
    G_Tmax2TnRatio = 1.2; % initial guess for ratio

    sensitivity.objScale = 1; % important to reset scale !
    sensitivity.objOffset = 0.0;
    try 

        fprintf('\n---------------------------------------------------------\n')
        fprintf('Starting simulation %i ',simNum)
        fprintf('\n---------------------------------------------------------\n')
        %read channel file and design parameters
        readOptions.figBounds = [mesh.boundary.xi, mesh.boundary.xf, ...
                                 mesh.boundary.yi, mesh.boundary.yf];
        readOptions.filePrefix = [G_simDirectory,G_filePrefix];
        readOptions.polygonFile = polygonFile;
        readOptions.polygonFig = [readOptions.filePrefix,'_polygons'];
        readOptions.nlconfun = @nonlinear_constraints;
        readOptions.nlcon.sinMinPolyAngle = sin(pi/180*minPolyAngle);
        readOptions.nlcon.minPolyArea = minPolyArea;
        readOptions.nlcon.sinMinSidePolyAngle = sin(minSidePolyAngle*pi/180);
        readOptions.nlcon.minSidePolyArea = minSidePolyArea;
   
        readOptions.mergeTriangles = mergeTriangles;

        if (simNum == 1 && ~randomizeFirst)
            readOptions.boundsFile = [];
            readOptions.sampleFile = [];
            readOptions.sampleNum = [];
         elseif ~isempty(boundsFile_sampleFile)
            if (nSimulations > 0 && jobNum_sampleNum > 0)
                readOptions.sampleFile = boundsFile_sampleFile;            
                readOptions.sampleNum = (jobNum_sampleNum - 1)*nSimulations + simNum;
                fprintf('Sample number %i ',readOptions.sampleNum)
                fprintf('\n---------------------------------------------------------\n')
                readOptions.boundsFile = [];
            elseif (nSimulations == 1 && jobNum_sampleNum < 0)
                readOptions.sampleFile = boundsFile_sampleFile;            
                readOptions.sampleNum = -jobNum_sampleNum;
                fprintf('Only sample number %i is chosen !',readOptions.sampleNum)
                fprintf('\n---------------------------------------------------------\n')
                readOptions.boundsFile = [];
            else
                warning('Sample file not used')
                readOptions.boundsFile = boundsFile_sampleFile;
                readOptions.sampleFile = [];
                readOptions.sampleNum = [];
            end
        end
        readOptions.Toffset = mesh.convect.Toffset;
        [channels, sensitivity.designParams, ...
         sensitivity.restrictedParams]...
            = preprocess_channels(channelFile,readOptions);
        channels.domainVol = mesh.domainVol;
        channels.domainArea = mesh.domainArea;
        G_out.designParams = sensitivity.designParams; % for output only

        write_channel_file([readOptions.filePrefix,'_ini.channel'],...
                            channels,...
                            mesh.convect.Toffset ,...
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
                                                   1,false,false); % don't scale obj yet

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

        if (scaleObj) 
             % IMPORTANT: i) must redeclare anonymous function after any of its arguments is changed
             %            ii) must ensure that the normalization factor is positive !
            sensitivity.objScale = abs(costfun(iniGuess));
            costfun = @ (delParams) fem_n_sensitivity (delParams,  ...
                                                       mesh, ...
                                                       channels, ...
                                                       femParams, ...
                                                       sensitivity, ...
                                                       1,true,true);

            G_out.objScale = sensitivity.objScale;  
        else
            G_out.objScale = 1.0;
        end
        sensitivity.objOffset = 0.0;
        G_out.objOffset = 0.0;

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

        %%
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
    end
    
    fclose ('all');
end % loop through all simulations
  
end


