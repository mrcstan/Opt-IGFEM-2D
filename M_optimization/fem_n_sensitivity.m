%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan in 2014
%%% Modified on 9/8/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function update the channels, perform finite element analysis and
% finaly provide all the needed data for sensitivity analysis

% INPUTS
%   delParams: change in design parameters
%              if this is empty, then channels are not updated
%   mesh: All the information about mesh and sensitivity analysis
%   channels: All the information about channels
%   iniChannels: Initial channel characteristics
%   femParams: Some miscellaneous variables needed in the code
%   sensitivity:
%   autoGrad:
%       1: auto-decide whether grad(s) need to be calculated based
%          on number of outputs, i.e, grad(s) are calculated only
%          when nargout == 2 or 4
%       2: calculate grad(s)
%       3: don't calculate grads
%   scaleObj:
%       true: scale objective function
%   updateGout:
%       true: update the global variable G_out
% OUTPUT:
%   theta: Objective function
%   gradTheta: Derivative of objective function wrt design parameters
%   oneNorm: 1-norm of the temperature
%   gradOneNorm; Derivative of 1-norm of the temperature wrt design
%                  parameters
function [theta, gradTheta, oneNorm, gradOneNorm, Tmax] ...
                = fem_n_sensitivity(delParams,...
                                    iniMesh, ...
                                    iniChannels, ...
                                    femParams,...
                                    sensitivity,...
                                    autoGrad,...
                                    scaleObj,...
                                    updateGout)
global G_simDirectory
global G_filePrefix
global G_funcEvalCounter
global G_constraint
global G_out

filePrefix = [G_simDirectory,G_filePrefix];


%{
g = nonlinear_constraints(delParams,...
                          sensitivity.designParams,...
                          iniChannels,...
                          sensitivity.nlcon);
if (any(g > 0))
    warning('non-linear constraints violated, returning nan')
    theta = nan;
    if (nargout > 1)
        gradTheta = nan(sensitivity.designParams.nParams,1);
    end
    return
end
%}

%all_figs = findobj(0, 'type', 'figure');
%close(setdiff(all_figs,999)); % close all figures except the history figure
close all

% start from initial mesh
mesh = iniMesh;

try 
    switch autoGrad
        case 1
            if (nargout == 2 || nargout == 4)
                calcGrad = true;
            else
                calcGrad = false;
            end
        case 2
            calcGrad = true;
        case 3
            calcGrad = false;
        case -2 
            calcGrad = true; % called by nonlinear_constraint.m
        case -3 
            calcGrad = false; % called by nonlinear_constraint.m
        otherwise
            error('unknown auto gradient option')
    end
    invalidChannels = false;
    
    % Depending on whether delParams is (i) empty, (ii) contains inf or (iii) 
    % otherwise, the following actions are carried out:
    % (i) no geometrical, pressure and mass flow rate updates
    % (ii) only update the pressure and mass flow rate
    % (iii) update geometry, pressure and mass flow rate
    if isempty(delParams)
         channels = iniChannels;
    elseif any(isinf(delParams))    
        [channels,pressure,mass] = update_channels(delParams, ...
                                                  sensitivity.designParams, ...
                                                  sensitivity.restrictedParams, ...
                                                  iniChannels, ...
                                                  'updateMassPressure', ...
                                                  calcGrad);
        if (any(isnan(pressure)))
            warning('some nodal pressures are NaN, returning NaN ...')
            invalidChannels = true;
        end                                         
    else
        [channels,pressure,mass] = update_channels(delParams, ...
                                                  sensitivity.designParams, ...
                                                  sensitivity.restrictedParams, ...
                                                  iniChannels, ...
                                                  'add', ...
                                                  calcGrad);
        if (any(isnan(pressure)))
            warning('some nodal pressures are NaN, returning NaN ...')
            invalidChannels = true;
        end                                         
    end
    if (channels_self_intersections(channels,femParams.tol.channelSelfIntersect))
        warning('self-intersecting channels, returning NaN ...')
        invalidChannels = true;
    end
    
    if (invalidChannels)
        theta = nan;
        if (nargout > 1)
            gradTheta = nan(sensitivity.designParams.nParams+sensitivity.restrictedParams.nParams,1);
        end
        oneNorm = [];
        gradOneNorm = [];
        Tmax = [];
        return
    end
  
catch err
    fprintf('output channel error message \n')
    fcountStr = num2str(G_funcEvalCounter);
    fname = [filePrefix,'_channel_fcount_',fcountStr,'.err'];
    fid = fopen(fname,'w');
    fprintf(fid,getReport(err));
    fclose(fid);
    write_channel_file([filePrefix,'_fcount_',fcountStr,'.channel'],...
                        channels,mesh.convect.Toffset,sensitivity.designParams,'w',[]);
                    
    fig = figure;
    plot_mesh_curve(mesh.node.coords,mesh.elem.elem_node,[],channels)
    % save last channel configuration
    fname = [filePrefix,'_channel_fcount_',fcountStr];
    saveas(fig, fname, 'jpg');
    saveas(fig, fname, 'fig');
    
    theta = nan;
    if (nargout > 1)
        gradTheta = nan(sensitivity.designParams.nParams+sensitivity.restrictedParams.nParams,1);
    end
    oneNorm = [];
    gradOneNorm = [];
    Tmax = [];
    return                
end

G_funcEvalCounter = G_funcEvalCounter+1;
fprintf('\n------------------------------------------------------------\n')
fprintf('function evaluation number %i ',G_funcEvalCounter)
fprintf('\n------------------------------------------------------------\n')


if updateGout
    if ~exist('pressure','var')
        [pressure,mass] ...
            = network_pressure_mass_flow_rate(channels.contvty,...
                                              channels.nurbs,...
                                              channels.diams,...
                                              channels.heights,...
                                              channels.viscosity,...
                                              channels.inletEndPoint,...
                                              channels.massin,...
                                              channels.powerXdensity,...
                                              channels.pressureOutletEndPoint,...
                                              channels.pressureOut,...
                                              channels.crossSection);
    end
    G_out.inPressure = pressure(channels.inletEndPoint); % inlet pressure
    G_out.massin = max(mass);
    
end

%{
figure
plot_mesh_curve(mesh.node.coords, ...
                mesh.elem.elem_node, ...
                mesh.elem.heatSource, ...
                channels)

figure
plot_channel_network(channels.contvty, ...
                     channels.nurbs, ...
                     pressure/1000.0,mass*1e6*60/1065)

% save last channel configuration
fname = [filePrefix,'_last_channel'];
set(fig, 'PaperPositionMode', 'auto');
print('-djpeg ', '-r300', fname);
%
%save(fname, 'channels');

headers{1} = '-----------------------------------------------------------';
headers{2} = ['function evaluation ',num2str(G_funcEvalCounter)];
headers{3} = headers{1};
write_channel_file([filePrefix,'.allChannels'],...
                    channels,...
                    sensitivity.designParams,...
                    'a',headers)

%fname = [filePrefix,'_last_pressure_mass'];
%save(fname, 'pressure','mass');
%}

try
    % IMPORTANT: Substitute all occurences of otherFlags.calcItrsectVel with
    % calcGrad !
    if (calcGrad)
        fprintf('checking for original nodes that are intersection points \n')
        fprintf('and kinks or branch points on element edges \n')
        %nodeTimer = tic;
        [mesh.node.coords, nodesMoved] ...
                = eliminate_original_node_intersections(mesh.node.coords,...
                                                        mesh.edge.edge_node,...   
                                                        channels.nurbs,...
                                                        channels.branch_kinks.XX',...
                                                        mesh.boundary,...
                                                        femParams.tol.boundary,...
                                                        femParams.tol.node,...
                                                        femParams.moveNode.maxAttempts,...
                                                        femParams.moveNode.dist,...
                                                        femParams.moveNode.randDirection);
        %IMPORTANT: at this point, if mesh.node.coords is changed, it will
        % no longer be consistent with mesh.DT.Points. 
        % So mesh.DT should be updated before being used
        % However, if refinement is requested, mesh.elem.elem_node will be
        % ordered in a certain way for the bisection method.
        % So before the refinement, do not change mesh.elem.elem_node                                            
        %toc(nodeTimer);
    else
        nodesMoved = false;
    end
    
    fprintf('\nfinding intersection points \n')
    %intersectTimer = tic;
    [mesh.edge,...
     mesh.node,...
     mesh.elem,...
     channels.itrsectParams,...
     channels.origNodes, ...
     channels.enrichNodes, ...
     refineLevel]...
         =edges_curves_intersect(mesh.edge,mesh.node,...
                                 mesh.elem,...
                                 channels,...
                                 femParams.tol,...
                                 femParams.refine,...
                                 calcGrad);  
    if ~isempty(sensitivity.costFunction.objOpt.node)
        switch sensitivity.costFunction.objOpt.node 
            case -1
                % inlet nodes only
                sensitivity.costFunction.objOpt.node = mesh.node.inletNodes;
            case -2
                % oulet nodes only
                sensitivity.costFunction.objOpt.node = mesh.node.outletNodes;
            case -3
                % inlet and outlet nodes
                sensitivity.costFunction.objOpt.node = [mesh.node.inletNodes(:); ...
                                                        mesh.node.outletNodes(:)];
        end
    end
    if any(sensitivity.costFunction.objOpt.node > mesh.node.nOriginalNode)
        error('Nodal objective function has not been implemented for enrichment nodes')
    end
    if refineLevel || nodesMoved
        % reconstruct DT if refinement has been carried out or
        % the nodes have been moved
        [mesh.DT,mesh.elem] ...
            = update_mesh_DT_n_elem(mesh.node.coords(1:mesh.node.nOriginalNode,:), ...
                                    mesh.edge.edge_node, ...
                                    mesh.elem);
        mesh.elem.elem_edge = find_elem_edge(mesh.elem.elem_node, ...
                                             mesh.edge.edge_node, ...
                                             mesh.node.nOriginalNode);  
    elseif femParams.refine.maxRefineLevel
        % although no refinement has been carried out, elem_node has 
        % become inconsistent with DT after relabeling by the label
        % function for bisection. So, we need to 
        % revert it back to its original connectivity
        mesh.elem.elem_node = mesh.DT.ConnectivityList;
    end
    [mesh.elem.branch_kinks, ...
     mesh.node.coords, ...
     mesh.node.n_node, ...
     mesh.edge.itrsect] ...
                = elem_branching_n_kinks(mesh.DT,...
                                         mesh.elem.elem_edge,...
                                         mesh.edge.itrsect,...
                                         mesh.node.coords,...
                                         mesh.node.n_node,...
                                         channels.branch_kinks,...
                                         channels.itrsectParams,...
                                         femParams.tol.nurbsParam, ...
                                         femParams.tol.bary, ...
                                         channels.designParamNum, ...
                                         calcGrad);
    % nOriginalEnrichNode only includes enrichment nodes that arise due to 
    % the intersections
    % If NURBS-IGFEM is used, additional enrichment nodes corresponding to
    % additional control points may arise
    mesh.node.nOriginalEnrichNode = mesh.node.n_node;                                  
    %fprintf('\nintersection time = %g \n',toc(intersectTimer));

    if (mesh.node.n_node == mesh.node.nOriginalNode)
        disp('conforming mesh')
    else
        disp('nonconforming mesh')
    end
       
    fprintf('\nsetting boundary conditions \n')
    [mesh.elem, mesh.node] = set_boundary_conditions(mesh.BCs,...
                                                     mesh.elem,...
                                                     mesh.node,...
                                                     mesh.boundary,...
                                                     femParams.tol.boundary);

    % creat parent elements
    disp('constructing parent elements')
    %parentTimer = tic;
    mesh.elem.parent = element_intersections(mesh.elem.elem_node, ...
                                             mesh.elem.dualedge, ...
                                             mesh.edge.edge_node, ...
                                             mesh.edge.itrsect, ...
                                             calcGrad);
    %fprintf('\nparent time = %g \n', toc(parentTimer));
                                       

    [mesh.elem.parent, ...
     mesh.elem.cstrElems, ...
     mesh.elem.nIGFEMelems] ...
                     = parent_elements_nurbs(mesh.elem.parent,...
                                             mesh.elem.elem_node,...
                                             mesh.elem.material,...
                                             mesh.material,...
                                             mesh.elem.branch_kinks,...
                                             mesh.node.coords,...
                                             mesh.node.nOriginalNode,...
                                             mesh.node.constraint,...
                                             channels,...
                                             femParams.tol.nurbsParam,...
                                             femParams.otherFlags.polyIGFEM,...
                                             calcGrad);                                     
    %{
    figure
    elem2plot = [];
    plot_mesh_nurbs_parent(mesh.node.coords,mesh.node.nOriginalNode,...
                           mesh.elem.elem_node,mesh.elem.parent,channels, ...
                           false,false,elem2plot);  
    %}                 
    if(mesh.node.n_node==mesh.node.nOriginalNode)
        disp('conforming mesh')
    else
        disp('non-conforming mesh')
        disp('constructing child elements')
        %childTimer = tic;
        [mesh.elem.parent,mesh.node]...edge, Dirichlet, nodeCoords, boundary
            =child_elements_nurbs(mesh.elem.parent,...
                                  mesh.node,...
                                  femParams.tol.collinear, ...
                                  femParams.tol.slender,...
                                  femParams.otherFlags.polyIGFEM);
        %fprintf('\nchild time = %g \n', toc(childTimer));                      
        % the global equation number of additional equations for 
        % Lagrange multiplier method                     
        for i = 1:numel(mesh.elem.cstrElems)
            mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows...
                    =  mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows ...
                     + size(mesh.node.coords,1); 
        end                        

    end

    
    if updateGout
        G_out.elem = mesh.elem;
        G_out.node = mesh.node;
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %  Finite Element Code  %
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    % gauss quadrature schemes

    disp('performing FEM')
    % eq_num:  Equation number assigned to each node
    % n_dof: The number of degree of freedom in the model
    nLagrangeMult = numel(mesh.node.constraint.temp_node) ...
                    +size(mesh.node.constraint.nodePairs,1);
    [n_dof, eq_num] = initialize(size(mesh.node.coords,1) + nLagrangeMult, ...
                                 mesh.node.Dirichlet.n_pre_temp, ...
                                 mesh.node.Dirichlet.temp_node);

    [UUR, ~,fGloInd,~,Kff,~,Jac] = assemble_and_solve(mesh,...
                                    eq_num,...
                                    femParams.gauss,...
                                    channels,...
                                    femParams.otherFlags.supg,...
                                    femParams.otherFlags.polyIGFEM,...
                                    n_dof,...
                                    calcGrad,...
                                    []);
    %{
    fprintf('saving variables\n')
    fname = [filePrefix,'_fcount_',num2str(G_funcEvalCounter)];
    save(fname)
    write_channel_file([fname,'.channel'],...
                        channels,sensitivity.designParams,'w',[]);
    %}                
    % output paraview file
    UUR_update = update_enrichment_node_value(UUR, ...
                                              mesh.node, ...
                                              mesh.elem, ...
                                              mesh.edge);
    if nargout > 4
        Tmax = max(UUR_update);
    end
   
    if updateGout
        G_out.UUR = UUR;
        G_out.UUR_update = UUR_update;
        G_out.Toffset = mesh.convect.Toffset;
        G_out.Tunit = mesh.convect.Tunit;
        G_out.inletNodes = mesh.node.inletNodes;
        G_out.outletNodes = mesh.node.outletNodes;
        G_out.channels = channels;
    end
    fprintf('\nsensitivity_analysis\n')
    %sensitivityTimer = tic;
    
    if strcmpi(sensitivity.costFunction.type,'AREA')
        theta = sum(channels.diams.*channels.length);
        if nargout > 1
            [D2L, D2diam] = D2_channel_lengths_n_diameters(...
                                    channels,sensitivity.designParams);
            gradTheta = nan(sensitivity.designParams.nParams,1);
            for i = 1:sensitivity.designParams.nParams
                gradTheta(i) = sum(D2L(:,i).*channels.diams ...
                               +channels.length.*D2diam(:,i));
            end
            gradTheta = [gradTheta;zeros(sensitivity.restrictedParams.nParams,1)];
        else
            gradTheta = [];
        end

        oneNorm = [];
        gradOneNorm = [];
    else
        if mesh.convect.linearRad
            [theta, gradTheta,oneNorm,gradOneNorm] ...
                = objective_n_gradient(mesh, ...
                                      channels, ...
                                      sensitivity.costFunction, ...
                                      sensitivity.designParams, ...
                                      UUR, ...
                                      sensitivity.gauss, ...
                                      femParams.otherFlags.supg, ...
                                      Kff, ... 
                                      fGloInd, ...
                                      calcGrad);
        else
            [theta, gradTheta,oneNorm,gradOneNorm] ...
                = objective_n_gradient(mesh, ...
                                      channels, ...
                                      sensitivity.costFunction, ...
                                      sensitivity.designParams, ...
                                      UUR, ...
                                      sensitivity.gauss, ...
                                      femParams.otherFlags.supg, ...
                                      Jac, ... 
                                      fGloInd, ...
                                      calcGrad);
        end
        if calcGrad
            gradTheta = [gradTheta;zeros(sensitivity.restrictedParams.nParams,1)];
            if sensitivity.costFunction.objOpt.calcOneNorm
                gradOneNorm = [gradOneNorm;zeros(sensitivity.restrictedParams.nParams,1)];
            end

        end
   
    end
    %toc(sensitivityTimer)
    output_values([filePrefix,'.objGrad'], ...
                  'Counter, # obj, Obj, Gradient', ...
                  G_funcEvalCounter,[numel(theta),theta(:)',gradTheta(:)'])
    
    if scaleObj
        if (isfield(sensitivity,'objScale'))
            fprintf('objective value and gradient scaled\n')
            theta = (theta - sensitivity.objOffset)/sensitivity.objScale;
            gradTheta = gradTheta/sensitivity.objScale;
        else
            error('objective requested to be scale but scale not given')
        end
    end         
    if autoGrad == -2 || autoGrad == -3
        % Only output constraints right after they are evaluated
        output_values([filePrefix,'.cstr'], ...
                      'Counter, Constraints', G_funcEvalCounter,G_constraint.g(:)')
        output_values([filePrefix,'.inOutTcstr'], ...
                      'Counter, Applicable inlet and/or outlet temperature constraints', ...
                      G_funcEvalCounter,G_constraint.history.nodalT - mesh.convect.Toffset)
    end
catch err   
    fprintf('output error message \n')
    
    fname = [filePrefix,'_fcount_',num2str(G_funcEvalCounter),'.err'];
    fid = fopen(fname,'w');
    fprintf(fid,getReport(err));
    fclose(fid);
    fcountStr = num2str(G_funcEvalCounter);
    fprintf('output channel \n')
    write_channel_file([filePrefix,'_fcount_',fcountStr,'.channel'],...
                        channels,mesh.convect.Toffset,sensitivity.designParams,'w',[]);
                    
    fig = figure;
    options.showMesh = true;
    plot_mesh_curve(mesh.node.coords,mesh.elem.elem_node,[],channels,options)
    % save last channel configuration
    fname = [filePrefix,'_fcount_',fcountStr];
    saveas(fig, fname, 'fig');
    saveas(fig, fname, 'jpg');
    
    if (exist('mesh','var'))
        fprintf('saving mesh ...\n')
        save([filePrefix,'_last_mesh_w_err'],'mesh')
    end
    if (exist('solution','var'))
        fprintf('saving solution ...\n')
        save([filePrefix,'_last_soln_w_err'],'solution')
    end
    
    rethrow(err)
    %theta = nan;
    %gradTheta = [];
end


end

