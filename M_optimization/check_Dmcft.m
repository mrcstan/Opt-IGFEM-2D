%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Ahmad Raeisi Najafi on 10/08/2012
%%% Last modified date: 03/28/2014
%%% Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Function and sensitivity evaluation module for 
%                 Optimization analysis

% This function update the channels, perform finite element analysis and
% finaly provide all the needed data for sensitivity analysis

%INPUTS
% delParams: change in design parameters
% mesh: All the information about mesh and sensitivity analysis
% channels: All the information about channels
% iniChannels: Initial channel characteristics
% femParams: Some miscellaneous variables needed in the code

%OUTPUT
% theta: Objective and constraint functions
% gradTheta: Derivative of objective and constraint functions
% n_funct: Number of functions (including all the objective and
%          constraints function)

function [mcft,Dmcft] = check_Dmcft(delParams,...
                                    iniMesh, ...
                                    iniChannels, ...
                                    femParams,...
                                    sensitivity)



calcItrsectVel = true;

% start from initial mesh
mesh = iniMesh;

[channels,pressure] = update_channels(sensitivity.designParams, ...
                                      delParams, ...
                                      iniChannels, ...
                                      'add', ...
                                      false);
                                  

plot_mesh_curve(mesh.node.coords,mesh.elem.elem_node,channels,false,false)


fprintf('\nfinding intersection points \n')
% DEBUG DEBUG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%{
if (calcItrsectVel)
    tic
    mesh.node.coords ...
            = eliminate_original_node_intersections(mesh.node.coords,...
                                                    channels.nurbs,...
                                                    mesh.boundary,...
                                                    femParams.tol.boundary,...
                                                    femParams.tol.node,...
                                                    femParams.moveNode.maxAttemps,...
                                                    femParams.moveNode.del);
    toc
end
%}
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

%plot_mesh_curve_itrsect_junc(mesh.node,mesh.elem.elem_node,channels,...
%               mesh.edge.itrsect,mesh.elem.junc,mesh.elem.kinks,false,false);

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
%                     mesh.elem.elem_node,mesh.elem.parent,channels,true,true);                                      

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

nChannels = numel(channels.nurbs);
mcft = cell(nChannels,1);
for i = 1:numel(nChannels)
    nLineSegs = channels.nurbs(i).number - 1; 
    mcft{i} = zeros(2,nLineSegs);
    for j = 1:nLineSegs
        mcft{i}(:,j) = channels.nurbs(i).coefs(1:2,j+1) - channels.nurbs(i).coefs(1:2,j);
        mcft{i}(:,j) = mcft{i}(:,j)/norm(mcft{i}(:,j)); 
    end
    mcft{i} = channels.mcf(i)*mcft{i};
end

if (nargout > 1) 
    el = 2;
    ch = 1;
    lineSegNum = 1;
    Dmcft = cell(nChannels,1);
    for i = 1:nChannels
        nLineSegs = channels.nurbs(i).number - 1;
        Dmcft{i} = zeros(2,sensitivity.designParams.nParams,nLineSegs);
    end
    for i = 1:numel(mesh.elem.parent(el).child(ch).channelNum)
        pts = mesh.node.coords(mesh.elem.parent(el).child(ch).channelNodes{i},:)';
        [DtDpts,channelLocPaEnNodesOffset] ...
        = channel_unit_vec_der_wrt_end_pts(mesh.elem.parent(el).child(ch).locPaNodes,...
                                           mesh.elem.parent(el).child(ch).channelLocNodes(:,i),...
                                           pts)
        chanNum = mesh.elem.parent(el).child(ch).channelNum(i);
        for j = 1:numel(mesh.elem.parent(el).designParamNum)
            paramNum = mesh.elem.parent(el).designParamNum(j);
            for k = 1:numel(channelLocPaEnNodesOffset)
                Dmcft{chanNum}(:,paramNum,lineSegNum) = Dmcft{chanNum}(:,paramNum,lineSegNum) ...
                    + DtDpts(1:2,(2*k-1):2*k)*mesh.elem.parent(el).vel(:,channelLocPaEnNodesOffset(k),j);
            end
        end
        Dmcft{chanNum} = channels.mcf(chanNum)*Dmcft{chanNum}; 
    end
end  

end

function [DtDpts,channelLocPaEnNodesOffset] ...
            = channel_unit_vec_der_wrt_end_pts(locPaNodes, ...
                                               channelLocNodes , ...
                                               pts)
chanLocNodesOffset = locPaNodes(channelLocNodes) - 3;
vecLen = norm(pts(:,2) - pts(:,1));
vecLenCube = vecLen^3;
if (all(chanLocNodesOffset > 0))
    DtDpts = zeros(2,4);
    DtDpts(1,1) = -(pts(2,1) - pts(2,2))^2;
    DtDpts(2,1) = (pts(2,1) - pts(2,2))*(pts(1,1) - pts(1,2));
    DtDpts(1,2) = DtDpts(2,1);
    DtDpts(2,2) = -(pts(1,1) - pts(1,2))^2;
    DtDpts(1,3) = -DtDpts(1,1);
    DtDpts(2,3) = -DtDpts(2,1);
    DtDpts(1,4) = -DtDpts(2,1);
    DtDpts(2,4) = -DtDpts(2,2);
    DtDpts = DtDpts/vecLenCube;
    channelLocPaEnNodesOffset = zeros(2,1);
    channelLocPaEnNodesOffset(1) = chanLocNodesOffset(1);
    channelLocPaEnNodesOffset(2) = chanLocNodesOffset(2);
elseif (chanLocNodesOffset(1) > 0)
    DtDpts = zeros(2,2);
    DtDpts(1,1) = -(pts(2,1) - pts(2,2))^2;
    DtDpts(2,1) = (pts(2,1) - pts(2,2))*(pts(1,1) - pts(1,2));
    DtDpts(1,2) = DtDpts(2,1);
    DtDpts(2,2) = -(pts(1,1) - pts(1,2))^2;
    channelLocPaEnNodesOffset = chanLocNodesOffset(1);
elseif (chanLocNodesOffset(2) > 0)
    DtDpts = zeros(2,2);
    DtDpts(1,1) = (pts(2,1) - pts(2,2))^2;
    DtDpts(2,1) = (pts(2,1) - pts(2,2))*(pts(1,2) - pts(1,1));
    DtDpts(1,2) = DtDpts(2,1);
    DtDpts(2,2) = (pts(1,1) - pts(1,2))^2;
    channelLocPaEnNodesOffset = chanLocNodesOffset(2);
else
    DtDpts = [];
    channelLocPaEnNodesOffset = [];
end
    
end