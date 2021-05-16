%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan in 2014
%%% Modified on 2/21/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function update the channels, perform finite element analysis and
% finaly provide all the needed data for sensitivity analysis

%INPUTS
% delParams: change in design parameters
% mesh: All the information about mesh and sensitivity analysis
% channels: All the information about channels
% iniChannels: Initial channel characteristics
% femParams: Some miscellaneous variables needed in the code

%OUTPUT
% varargout{1:n} = {theta1,theta2,...,thetan}, n = number of objective functions 
% varargout{(n+1):2n} = {gradTheta1,gradTheta2,...,gradThetan}
% Note: n can be determine from the number of elements of
%       sensitivity.costFunction
function varargout = fem_n_sensitivity(delParams,...
                                    iniMesh, ...
                                    iniChannels, ...
                                    femParams,...
                                    sensitivity)
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

all_figs = findobj(0, 'type', 'figure');
close(setdiff(all_figs,999)); % close all figures except the history figure

calcItrsectVel = true;

% start from initial mesh
mesh = iniMesh;

% determine number of objective functions
nObjFuncs = numel(sensitivity.costFunction);

try 
    [channels,pressure] = update_channels(delParams, ...
                                          sensitivity.designParams, ...
                                          sensitivity.restrictedParams, ...
                                          iniChannels, ...
                                          'add');

    if (channels_self_intersections(channels,femParams.tol.channelSelfIntersect))
        warning('self-intersecting channels, returning nan')
        [varargout{1:nObjFuncs}] = deal(nan);
        [varargout{(nObjFuncs+1):2*nObjFuncs}] ...
            = deal(nan(sensitivity.designParams.nParams+sensitivity.restrictedParams.nParams,1));        
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
                        channels,sensitivity.designParams,'w',[]);
                    
    fig = figure('visible','off');
    plot_mesh_curve(mesh.node.coords,mesh.elem.elem_node,channels,false,false)
    % save last channel configuration
    fname = [filePrefix,'_channel_fcount_',fcountStr];
    set(fig, 'PaperPositionMode', 'auto');
    print('-djpeg ', '-r300', fname); 
    saveas(fig, fname, 'fig');
    
    [varargout{1:nObjFuncs}] = deal(nan);
    [varargout{(nObjFuncs+1):2*nObjFuncs}] ...
        = deal(nan(sensitivity.designParams.nParams+sensitivity.restrictedParams.nParams,1));          
    return                
end

G_funcEvalCounter = G_funcEvalCounter+1;
counter = G_funcEvalCounter;
fprintf('\n------------------------------------------------------------\n')
fprintf('function evaluation number %i ',G_funcEvalCounter)
fprintf('\n------------------------------------------------------------\n')


G_out.inPressure = pressure(1); % inlet pressure
G_out.channels = channels;
%{
fig = figure('visible','off');
plot_mesh_curve(mesh.node.coords,mesh.elem.elem_node,channels,false,false)

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
    if (calcItrsectVel)
        fprintf('checking for original nodes that are intersection points \n')
        fprintf('and kinks or branch points on element edges \n')
        nodeTimer = tic;
        mesh.node.coords ...
                = eliminate_original_node_intersections(mesh.node.coords,...
                                                        mesh.edge.edge_node,...   
                                                        channels.nurbs,...
                                                        [channels.kinks.x;
                                                        channels.kinks.y]',...
                                                        channels.pts(channels.junc,:),...
                                                        mesh.boundary,...
                                                        femParams.tol.boundary,...
                                                        femParams.tol.node,...
                                                        femParams.moveNode.maxAttemps,...
                                                        femParams.moveNode.dist);
        toc(nodeTimer);
    end

    
    fprintf('\nfinding intersection points \n')
    intersectTimer = tic;
    [mesh.edge,mesh.node,mesh.elem]...
         =edges_curves_intersect(mesh.edge,mesh.node,...
                                 mesh.elem,...
                                 channels,...
                                 femParams.tol,...
                                 femParams.opt,...
                                 calcItrsectVel);                       
    fprintf('\nintersection time = %g \n',toc(intersectTimer));

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
    parentTimer = tic;
    mesh.elem.parent = element_intersections(mesh.elem.elem_node, ...
                                             mesh.elem.dualedge, ...
                                             mesh.edge.edge_node, ...
                                             mesh.edge.itrsect, ...
                                             calcItrsectVel);
    fprintf('\nparent time = %g \n', toc(parentTimer));
                                       

    [mesh.elem.parent, ...
     mesh.elem.cstrElems, ...
     mesh.elem.nIGFEMelems] ...
                     = parent_elements_nurbs(mesh.elem.parent,...
                                             mesh.elem.elem_node,...
                                             mesh.elem.elem_edge,...
                                             mesh.elem.material,...
                                             mesh.material,...
                                             mesh.elem.junc,...
                                             mesh.elem.kinks,...
                                             mesh.node.coords,...
                                             mesh.node.nOriginalNode,...
                                             mesh.node.constraint,...
                                             channels,...
                                             femParams.tol.nurbsParam,...
                                             femParams.polyIGFEM,...
                                             calcItrsectVel);                                     

    %plot_mesh_nurbs_parent(mesh.node.coords,mesh.node.nOriginalNode,...
    %                     mesh.elem.elem_node,mesh.elem.parent,channels,false,false);  
                     
    if(mesh.node.n_node==mesh.node.nOriginalNode)
        disp('conforming mesh')
    else
        disp('non-conforming mesh')
        disp('constructing child elements')
        childTimer = tic;
        [mesh.elem.parent,mesh.node]...edge, Dirichlet, nodeCoords, boundary
            =child_elements_nurbs(mesh.elem.parent,...
                                  mesh.node,...
                                  femParams.slenderTol,...
                                  femParams.polyIGFEM);
        fprintf('\nchild time = %g \n', toc(childTimer));                      
        % the global equation number of additional equations for 
        % Lagrange multiplier method                     
        for i = 1:numel(mesh.elem.cstrElems)
            mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows...
                    =  mesh.elem.parent(mesh.elem.cstrElems(i)).cstrRows ...
                     + size(mesh.node.coords,1); 
        end                        

    end
    %{
    fprintf('saving mesh ...\n')
    savetimer = tic;
    save([filePrefix,'_last_mesh'],'mesh')
    toc(savetimer)
    %}
    G_out.elem = mesh.elem;
    G_out.node = mesh.node;
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %  Finite Element Code  %
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    % gauss quadrature schemes

    disp('performing FEM')
    % eq_num:  Equation number assigned to each node
    % n_dof: The number of degree of freedom in the model
    nLagrangeMult = numel(mesh.node.constraint.temp_node);
    [n_dof, eq_num] = initialize(size(mesh.node.coords,1) + nLagrangeMult, ...
                                 mesh.node.Dirichlet.n_pre_temp, ...
                                 mesh.node.Dirichlet.temp_node);

    % Assemble the stiffness matrix
    if (femParams.polyIGFEM)
          %
          Up =  assemble_prescribed_node(mesh.node.Dirichlet,eq_num);
          [Kff,Kfp,Kpf,Kpp,Pf,Pp] ...
                        = mx_assemble_sparse(mesh.node.coords',...
                                             mesh.elem.elem_node',...
                                             mesh.elem.heatSource,...
                                             mesh.convect,...
                                             eq_num,...
                                             femParams.gauss,...
                                             mesh.elem.parent,...
                                             channels,...
                                             mesh.node.Dirichlet,...
                                             mesh.elem.Neumann,...
                                             femParams.supg);

        %{
        [Kpp,Kpf,Kfp,Kff,Pp,Pf,Up] = ...
            assemble (mesh, eq_num, n_dof,channels,femParams.gauss,femParams.supg);
        %}
    else
        [Kpp,Kpf,Kfp,Kff,Pp,Pf,Up] = ...
            assemble (mesh, eq_num, n_dof,channels,femParams.gauss,femParams.supg);
    end
    fprintf('\nsolving equation\n')
    [UUR, ~,fGloInd,~] ...
        = solve_matrix_eqn(Kpp,Kpf,Kfp,Kff,Pp,Pf,Up, eq_num); 




    % output paraview file
    UUR_update = update_enrichment_node_value(UUR, ...
                                              mesh.node, ...
                                              mesh.elem, ...
                                              mesh.edge);
    %{
    fprintf('\nsaving solution ...\n')
    savetimer = tic;
    save([filePrefix,'_last_soln'],'solution') 
    toc(savetimer)
    %}
    G_out.UUR = UUR;
    G_out.UUR_update = UUR_update;

    fprintf('\nsensitivity_analysis\n')
    sensitivityTimer = tic;
    
    if (nargout == 2*nObjFuncs)
        calcGrad = true;
    else
        calcGrad = false;
    end
    varargout = cell(2*nObjFuncs,1);
    for i = 1:nObjFuncs
        i2 = nObjFuncs + i;
        [varargout{i},varargout{i2}] ...
            = objective_n_gradient(mesh, ...
                                 channels, ...
                                 sensitivity.costFunction(i), ...                                            
                                 UUR, ...
                                 UUR_update, ...
                                 nLagrangeMult, ...
                                 sensitivity.gauss, ...
                                 femParams.supg, ...
                                 Kff, ... 
                                 fGloInd, ...
                                 calcGrad);
        if (sensitivity.restrictedParams.nParams && calcGrad)
            varargout{i2} = [varargout{i2};zeros(sensitivity.restrictedParams.nParams,1)];
        end
        if (isfield(sensitivity,'objScale'))
            fprintf('objective value and gradient scaled\n')
            varargout{i} = varargout{i}/sensitivity.objScale;

            varargout{i2} = varargout{i2}/sensitivity.objScale;
        end
    end
    toc(sensitivityTimer)
    % just output the first objective function and its gradient
    output_opt (delParams, [varargout{1};G_constraint], ...
                varargout{2},[filePrefix,'.opt'],counter);

catch err   
    fprintf('output error message \n')
    
    fname = [filePrefix,'_fcount_',num2str(G_funcEvalCounter),'.err'];
    fid = fopen(fname,'w');
    fprintf(fid,getReport(err));
    fclose(fid);
    fcountStr = num2str(G_funcEvalCounter);
    fprintf('output channel \n')
    write_channel_file([filePrefix,'_fcount_',fcountStr,'.channel'],...
                        channels,sensitivity.designParams,'w',[]);
                    
    fig = figure('visible','off');
    plot_mesh_curve(mesh.node.coords,mesh.elem.elem_node,channels,false,false)
    % save last channel configuration
    fname = [filePrefix,'_fcount_',fcountStr];
    set(fig, 'PaperPositionMode', 'auto');
    print('-djpeg ', '-r300', fname); 
    saveas(fig, fname, 'fig');
    
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

