%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan in 2014
%%% Last modified date: 9/23/2015
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear constraint function for parallel channels
%INPUT:


%OUTPUT
% g and gradg: Nonlinear inequalitis
% geq and gradgeq: Nonlinear equalitis

function [ g, geq, gradg, gradgeq] ...
                = nonlinear_constraints(delParams, ...
                                        designParams,...
                                        restrictedParams, ...
                                        iniChannels, ...
                                        nlcon, ...
                                        mesh, ...
                                        femParams, ...
                                        sensitivity)
% disp('NONLINEAR CONSTRAINTS !!!!!!!!!!!!!!!!!!!!!!!!!')
% VERY IMPORTANT TO REINITIALIZE THESE GLOBAL VARIABLES OUTSIDE OF THIS 
% FUNCTION AFTER EVERY SIMULATION !!!!!!!!!!!!!! 
global G_constraint
global G_Tmax2TnRatio

[channels,~,~] = update_channels(delParams, ...
                                 designParams, ...
                                 restrictedParams, ...
                                 iniChannels, ...
                                 'add', ...
                                 (nargout > 2));
                           
if ~isfield(channels,'polygons') || isempty(channels.polygons) || isempty(channels.vertexCoords)
    channels.polygons = [];
    nPolygons = 0;
    nV = 0;
    ng = 0;
    %warning('no polygons specified for nonlinear constraint')
else
    nPolygons = numel(channels.polygons);
    nV = sum(cat(1,channels.polygons.nVertices));
    %nVertexPairs = size(channels.vertexPairs,1);
    ng = nV + nPolygons;
end

if (~isfield(nlcon,'type'))
    nlcon.type = [];
end

PminConOn = false;
PmaxConOn = false;
TmaxConOn= false;
AminConOn = false;
AmaxConOn = false;
VminConOn = false;
VmaxConOn = false;
% flag for NODAL_T_OUT_max, NODAL_T_OUT_min,
%          NODAL_T_IN_max or NODAL_T_IN_min
TnodalConOn = false; 
if (ischar(nlcon.type))
    splitStr = regexp(nlcon.type,',','split');
    for i = 1:numel(splitStr)
        switch splitStr{i}
             case 'Pmin'
                if (~isfield(nlcon,'minP'))
                    error('nlcon must have maxP field when nlcon.type=="P"')
                end
                PminConOn = true;
            case 'Pmax'
                if (~isfield(nlcon,'maxP'))
                    error('nlcon must have maxP field when nlcon.type=="P"')
                end
                PmaxConOn = true;
            case 'Tmax'
                if (~isfield(nlcon,'maxTmax'))
                    error('nlcon must have maxTmax field when nlcon.type=="Tmax"')
                end
                TmaxConOn = true;
            case 'Amin'
                if (~isfield(nlcon,'area'))
                    error('nlcon must have area field when nlcon.type=="Amin"')
                end
                if (~isfield(nlcon,'minAorVFrac'))
                    error('nlcon must have minAorVFrac field when nlcon.type=="Amin"')
                end
                AminConOn = true;
            case 'Amax'
                if (~isfield(nlcon,'area'))
                    error('nlcon must have area field when nlcon.type=="Amax"')
                end
                if (~isfield(nlcon,'maxAorVFrac'))
                    error('nlcon must have maxAorVFrac field when nlcon.type=="Amax"')
                end
                AmaxConOn = true;
            case 'Vmin'
                if (~isfield(channels,'domainVol'))
                    error('channels must have domainVol field when nlcon.type=="Vmin"')
                end
                if (~isfield(nlcon,'minAorVFrac'))
                    error('nlcon must have minAorVFrac field when nlcon.type=="Vmin"')
                end
                VminConOn = true;    
            case 'Vmax'
                if (~isfield(channels,'domainVol'))
                    error('channels must have domainVol field when nlcon.type=="Vmax"')
                end
                if (~isfield(nlcon,'maxAorVFrac'))
                    error('nlcon must have maxAorVFrac field when nlcon.type=="Vmax"')
                end
                VmaxConOn = true;
            case 'nodalT'
                if ~isfield(nlcon,'nodalTbounds')
                    error('nlcon must have nodalTbounds when nlcon.type=="NODAL_T"')
                end
                if numel(nlcon.nodalTbounds) ~= 4
                    error('nlcon.nodalTbounds must have 4 entries (set irrelevant entries to nan')
                end
                TnodalConOn = true;
                % nodalTbounds is a 4-entry vector with the first to fourth
                % respectively denoting the lower and upper bounds of the
                % inlet temperature and that of the outlet temperature,
                % respectively           
                % If the bound is not applicable set to nan
                if ~isfield(nlcon,'massin')
                    error('nlcon must have massin when nlcon.type=="NODAL_T"')
                end
                if numel(nlcon.massin) ~= 4
                    error('nlcon.massin must have 4 entries (set irrelevant entries to nan')
                end
                if ~isfield(nlcon,'inBCs')
                    error('nlcon must have inBCs when nlcon.type=="NODAL_T"')
                end
                if numel(nlcon.inBCs) ~= 4
                    error('nlcon.inBCs must have 4 entries (set irrelevant entries to nan')
                end
            otherwise
                error('unknown nonlinear constraint type')
        end
    end
end

%{
if ~isempty(sensitivity) && ~isempty(sensitivity.costFunction.type)
    if strcmpi(sensitivity.costFunction.type,'AREA')
        % update Gout in fem_n_sensitivity when it is called from this function
        % only when the objective function is area, since the calculation the
        % area objective function does not go through any FEM calculation
        updateGout = true;
    else
        updateGout = false;
    end
else
    updateGout = false;
end
%}

sensitivity.costFunction.objOpt.node = [];
if (PmaxConOn && TmaxConOn)
    Pind = ng + 1;
    ng = ng + 2;
    Tind = ng;
    sensitivity.costFunction.type = 'P_NORM';
    sensitivity.costFunction.objOpt.normp = 8;
    sensitivity.costFunction.objOpt.calcOneNorm = false;
    sensitivity.costFunction.objOpt.intDomainType = 0;
    if TnodalConOn 
        error('Pressure, Tmax and nodal temperatures constraints not yet implemented')
    end
elseif PmaxConOn || PminConOn  
    Pind = ng + 1;
    ng = ng + 1;
    sensitivity.costFunction.type = 'PRESSURE';
    sensitivity.costFunction.objOpt.normp = 1;
    sensitivity.costFunction.objOpt.calcOneNorm = false;
    sensitivity.costFunction.objOpt.intDomainType = 1;
    if TnodalConOn 
        error('Pressure and nodal temperatures constraints not yet implemented')
    end
elseif TmaxConOn
    ng = ng + 1;
    Tind = ng;
    sensitivity.costFunction.type = 'P_NORM';
    sensitivity.costFunction.objOpt.normp = 8;
    sensitivity.costFunction.objOpt.calcOneNorm = false;
    sensitivity.costFunction.objOpt.intDomainType = 0;
    if TnodalConOn 
        error('Tmax and nodal temperatures constraints not yet implemented')
    end
end
if  TnodalConOn 
    TnodalInd = ng + 1;
    ng = ng + nnz(~isnan(nlcon.nodalTbounds));
    if all(isnan(nlcon.massin)) && all(isnan(nlcon.inBCs))
        evalNodeByNode = false;
    else
        evalNodeByNode = true;
    end
end

if AminConOn
    ng = ng + 1;
    AminInd = ng;
elseif VminConOn
    ng = ng + 1;
    VminInd = ng;
end

if AmaxConOn
    ng = ng + 1;
    AmaxInd = ng;
elseif VmaxConOn
    ng = ng + 1;
    VmaxInd = ng;
end


% normalized normal constraint for generating Pareto front
if isfield(nlcon,'NNC')
    nObjs = size(nlcon.NNC.anchorPts,1); 
    if (nObjs ~= size(nlcon.NNC.anchorPts,2) ...
        || nObjs ~= nnz(nlcon.NNC.objInds))
        error('number of non-zero entries of objInds must equal to each dimension of the anchorPts matrix')
    end
    NNCind = ng;
    ng = ng + nObjs - 1;    
    normalizedObj = nan(nObjs,1);
    if (nargout > 2)
        dNormalizedObj = zeros(designParams.nParams+restrictedParams.nParams,nObjs);
    end
end

g = zeros(ng, 1);
geq = [];

% angle constraints
if (nargout > 2)
    gradg = zeros(designParams.nParams+restrictedParams.nParams,ng);
    gradgeq = [];
    i1 = 0;
    for i = 1:nPolygons
        designParamNum = designParams.vertices2params(channels.polygons(i).connectivity,:)'; 
        ind = ~isnan(designParamNum);
        designParamNum(~ind) = [];
        
        XX = channels.vertexCoords(channels.polygons(i).connectivity,:)';
        i2 = i1 + channels.polygons(i).nVertices;        
        [sineAngles,dSineAngles] = polygon_sine_angles(XX(1,:),...
                                                       XX(2,:)); 
        if (channels.polygons(i).isSideTriangle)
            g((i1+1):i2) = nlcon.sinMinSidePolyAngle - sineAngles; 
        else
            g((i1+1):i2) = nlcon.sinMinPolyAngle - sineAngles; 
        end
        gradg(designParamNum,(i1+1):i2) = -dSineAngles(ind(:),:);    
        i1 = i2;
        
        [area,dArea] = polygon_area(XX(1,:),...
                                    XX(2,:));
        if (channels.polygons(i).isSideTriangle)
            g(nV+i) = nlcon.minSidePolyArea - area;
        else
            g(nV+i) = nlcon.minPolyArea - area;
        end
        gradg(designParamNum,nV+i) = -dArea(ind(:));

    end
    if (PmaxConOn && TmaxConOn)
        % It is more efficient to output both the n-norm of the temperature
        % over the whole domain and the 1-norm of the temperature over the
        % channels in one pass to the function fem_n_sensitivity
        % However, this requires modification of the
        % mx_assemble_pseudo_adjoint_forces.cpp.
        %
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn,D2Tn,~, ~, Tmax] = fem_n_sensitivity ([],...
                                                          mesh, ...
                                                          channels, ...
                                                          femParams,...
                                                          sensitivity,...
                                                          -2, ...
                                                          false, ... % don't scale objective
                                                          false); % Tave is not scaled by area
                                                      
        if (isnan(Tn))
            warning('pressure and max temperature constraints returning NaN') 
            g(:) = NaN;
            gradg(:) = NaN;
            return;
        end               
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = 'PRESSURE';
        sensitivity.costFunction.objOpt.intDomainType = 1;
        [delP, DdelP] = fem_n_sensitivity([], ...
                                          mesh, ...
                                          channels, ...
                                          femParams,...
                                          sensitivity, ...
                                          -2, ...
                                          false, ...
                                          false);
       
        g(Pind) = delP/nlcon.maxP - 1.0;
        gradg(:,Pind) = DdelP/nlcon.maxP;
        % For debugging
        %g(Tind) = Tn - 1.0;
        %gradg(:,Tind) = D2Tn;
        %
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        gradg(:,Tind) = G_Tmax2TnRatio*D2Tn/nlcon.maxTmax;
        G_Tmax2TnRatio = Tmax/Tn;
        %
    elseif PmaxConOn || PminConOn   
        [delP, DdelP] = fem_n_sensitivity([], ...
                                          mesh, ...
                                          channels, ...
                                          femParams,...
                                          sensitivity, ...
                                          -2, ...
                                          false, ...
                                          false);
         if PminConOn
              g(Pind) = 1 - delP/nlcon.minP;
              gradg(:,Pind) = -DdelP/nlcon.minP;
         else
              g(Pind) = delP/nlcon.maxP - 1.0;
              gradg(:,Pind) = DdelP/nlcon.maxP;
         end
    elseif TmaxConOn
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn, D2Tn, ~, ~, Tmax] = fem_n_sensitivity([],...
                                                   mesh, ...
                                                   channels, ...
                                                   femParams,...
                                                   sensitivity, ...
                                                   -2, ...
                                                   false, ...
                                                   false); 
        if (isnan(Tn))
            warning('max temperature constraint returning NaN') 
            g(:) = NaN;
            gradg(:) = NaN;
            return;
        end     
        % For debugging
        %{
        g(Tind) = Tn/nlcon.maxTmax - 1.0;
        gradg(:,Tind) = D2Tn/nlcon.maxTmax;
        %}
        %
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        gradg(:,Tind) = G_Tmax2TnRatio*D2Tn/nlcon.maxTmax;        
        G_Tmax2TnRatio = Tmax/Tn;
        %
    end
    if TnodalConOn 
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = '';
        sensitivity.costFunction.objOpt.intDomainType = 0;
        %historyNodalTrow = size(G_constraint.history.nodalT,1) + 1;
        historyNodalTrow = 1;
        historyNodalTcol = 0;
        if evalNodeByNode
            gInd = TnodalInd;
            originalMassIn = channels.massin;
            originalPairIndHeat = channels.pair_ind_heat;
            originalPtTemp = channels.pt_temp;
            % Evaluate node by node since each constraint is subject to a
            % different flow rate and channel BC
            for i = 1:numel(nlcon.nodalTbounds)
                if ~isnan(nlcon.nodalTbounds(i))
                    if i > 2 
                        sensitivity.costFunction.objOpt.node = -2;
                    else
                        sensitivity.costFunction.objOpt.node = -1;
                    end
                     % Assume that all the inlets take the same value
                    if ~isnan(nlcon.massin(i))
                        % change the mass flow rate if applicable
                        channels.massin = nlcon.massin(i);
                    end
                    if ~isnan(nlcon.inBCs(i))
                        % change the temperature difference or the inlet
                        % temperature depending which one is specified
                        % originally
                        if size(channels.pair_ind_heat,2) == 3
                            channels.pair_ind_heat(:,3) = nlcon.inBCs(i);
                        elseif size(channels.pt_temp,2) == 2
                            channels.pt_temp(1,2) = nlcon.inBCs(i);
                        end
                    end
                    [nodalT, gradNodalT] = fem_n_sensitivity(inf, ...
                                                  mesh, ...
                                                  channels, ...
                                                  femParams,...
                                                  sensitivity, ...
                                                  -2, ...
                                                  false, ...
                                                  false);
                    historyNodalTcol = historyNodalTcol+1;                          
                    G_constraint.history.nodalT(historyNodalTrow,historyNodalTcol) = nodalT;
                    if rem(i,2) == 1
                       % lower bound
                       g(gInd) = nlcon.nodalTbounds(i) - nodalT;
                       gradg(:,gInd) = -gradNodalT;
                    else
                       % upper bound
                       g(gInd) = nodalT - nlcon.nodalTbounds(i);
                       gradg(:,gInd) = gradNodalT;
                   end                                
                                              
                    gInd = gInd + 1;    
                    channels.massin = originalMassIn;
                    channels.pair_ind_heat = originalPairIndHeat;
                    channels.pt_temp = originalPtTemp;
                end
            end
        else
            inletNodeFlag = any(~isnan(nlcon.nodalTbounds(1:2)));
            outletNodeFlag = any(~isnan(nlcon.nodalTbounds(3:4)));
            if inletNodeFlag && outletNodeFlag 
                % Inlet and outlet nodes required
                sensitivity.costFunction.objOpt.node = -3;
            elseif inletNodeFlag
                % Inlet nodes only
                sensitivity.costFunction.objOpt.node = -1;
            elseif outletNodeFlag
                % Outlet nodes only
                sensitivity.costFunction.objOpt.node = -2;
            end
            [nodalT, gradNodalT] = fem_n_sensitivity([], ...
                                                  mesh, ...
                                                  channels, ...
                                                  femParams,...
                                                  sensitivity, ...
                                                  -2, ...
                                                  false, ...
                                                  false);
            if sensitivity.costFunction.objOpt.node == -3 && numel(nodalT) < 2
                % Handle the case where at least two nodes cannot be found
                nodalT = nan(2,1);
                gradNodalT = nan(size(gradNodalT,1),2);
            end                                  
            gInd = TnodalInd;
            for i = 1:numel(nlcon.nodalTbounds)
               if ~isnan(nlcon.nodalTbounds(i))
                   if i > 2 && sensitivity.costFunction.objOpt.node == -3
                       nodalTind = 2;
                   else
                       nodalTind = 1;
                   end

                   if rem(i,2) == 1
                       % lower bound
                       g(gInd) = nlcon.nodalTbounds(i) - nodalT(nodalTind);
                       gradg(:,gInd) = -gradNodalT(:,nodalTind);
                   else
                       % upper bound
                       g(gInd) = nodalT(nodalTind) - nlcon.nodalTbounds(i);
                       gradg(:,gInd) = gradNodalT(:,nodalTind);
                   end      
                   gInd = gInd + 1;
               end
            end
        end
    end
    % area or volume fraction constraint
    if AminConOn || AmaxConOn 
        areaFrac = sum(channels.diams.*channels.length)/nlcon.area;
        [D2L, D2diam] = D2_channel_lengths_n_diameters(channels,designParams);
        dAreaFrac = nan(designParams.nParams,1);
        for i = 1:designParams.nParams
            dAreaFrac(i) = sum(D2L(:,i).*channels.diams ...
                           +channels.length.*D2diam(:,i))/nlcon.area;
        end
        if AminConOn
            g(AminInd) = nlcon.minAorVFrac - areaFrac;
            gradg(1:designParams.nParams,AminInd) = -dAreaFrac;
        end
        if AmaxConOn
            g(AmaxInd) = areaFrac - nlcon.maxAorVFrac;
            gradg(1:designParams.nParams,AmaxInd) = dAreaFrac;
        end
    elseif VminConOn || VmaxConOn
        volFrac = channels.vol/channels.domainVol;
        [D2L, D2diam] = D2_channel_lengths_n_diameters(channels,designParams);
        if (strcmpi(channels.crossSection,'circular'))
            factor = pi/(4*channels.domainVol);
        elseif (strcmpi(channels.crossSection,'square'))
            factor = 1/channels.domainVol;
        else
            error('unknown cross section')
        end
        dVolFrac = nan(designParams.nParams,1);
        for i = 1:designParams.nParams
            dVolFrac(i) = factor*sum(D2L(:,i).*channels.diams.^2 ...
                          +2*channels.length.*channels.diams.*D2diam(:,i));
        end
        if VminConOn
            g(VminInd) = nlcon.minAorVFrac - volFrac;
            gradg(1:designParams.nParams,VminInd) = -dVolFrac;
        end
        if VmaxConOn            
            g(VmaxInd) = volFrac - nlcon.maxAorVFrac;
            gradg(1:designParams.nParams,VmaxInd) = dVolFrac;
        end
    end
    
    if isfield(nlcon,'NNC')
        if (nlcon.NNC.objInds(1) && nlcon.NNC.objInds(2))
            % T and P
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = false;
            sensitivity.costFunction.objOpt.intDomainType = 0;
            [Tn,D2Tn] = fem_n_sensitivity ([],...
                                            mesh, ...
                                            channels, ...
                                            femParams,...
                                            sensitivity,...
                                            -2, ...
                                            false, ...
                                            false); % Tave is not scaled by area
            if (isnan(Tn))
                warning('pressure and max temperature constraints returning NaN') 
                g(:) = NaN;
                gradg(:) = NaN;
                return;
            end   
            
            sensitivity.costFunction.objOpt.normp = 1;
            sensitivity.costFunction.objOpt.calcOneNorm = false;
            sensitivity.costFunction.type = nlcon.NNC.obj2type;
            sensitivity.costFunction.objOpt.intDomainType = nlcon.NNC.obj2intDomainType;
            [delP, DdelP] = fem_n_sensitivity([], ...
                                          mesh, ...
                                          channels, ...
                                          femParams,...
                                          sensitivity, ...
                                          -2, ...
                                          false, ...
                                          false);
            % P
            normalizedObj(nlcon.NNC.objInds(2)) = delP;       
            % scale the gradient of pressure channels.D2P, 
            % which is evaluated at a reference viscosity and does not include
            % the effect of the change in viscosity with temperature
            dNormalizedObj(:,nlcon.NNC.objInds(2)) = DdelP; 
            % T       
            normalizedObj(nlcon.NNC.objInds(1)) = Tn;                                                
            dNormalizedObj(:,nlcon.NNC.objInds(1)) = D2Tn;      
        elseif nlcon.NNC.objInds(1)
            % T
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = false;   
            sensitivity.costFunction.objOpt.intDomainType = 0;
            [Tn, D2Tn] = fem_n_sensitivity([],...
                                           mesh, ...
                                           channels, ...
                                           femParams,...
                                           sensitivity, ...
                                           -2, ...
                                           false, ...
                                           false); 
            if (isnan(Tn))
                warning('max temperature constraint returning NaN') 
                g(:) = NaN;
                gradg(:) = NaN;
                return;
            end     
            normalizedObj(nlcon.NNC.objInds(1)) = Tn;
            dNormalizedObj(:,nlcon.NNC.objInds(1)) = D2Tn;
        elseif nlcon.NNC.objInds(2)
            % P
            sensitivity.costFunction.type = 'PRESSURE';
            sensitivity.costFunction.objOpt.normp = 1;
            sensitivity.costFunction.objOpt.calcOneNorm = false;
            sensitivity.costFunction.type = nlcon.NNC.obj2type;
            sensitivity.costFunction.objOpt.intDomainType = nlcon.NNC.obj2intDomainType;
            [delP, DdelP] = fem_n_sensitivity([], ...
                                              mesh, ...
                                              channels, ...
                                              femParams,...
                                              sensitivity, ...
                                              -2, ...
                                              false, ...
                                              false);
             normalizedObj(nlcon.NNC.objInds(2)) = delP;                      
             dNormalizedObj(:,nlcon.NNC.objInds(2)) = DdelP;
        end
        if nlcon.NNC.objInds(3)
            areaFrac = sum(channels.diams.*channels.length)/nlcon.area;
            [D2L, D2diam] = D2_channel_lengths_n_diameters(channels,designParams);
            dAreaFrac = nan(designParams.nParams,1);
            for i = 1:designParams.nParams
                dAreaFrac(i) = sum(D2L(:,i).*channels.diams ...
                               +channels.length.*D2diam(:,i))/nlcon.area;
            end
            normalizedObj(nlcon.NNC.objInds(3)) = areaFrac;
            dNormalizedObj(:,nlcon.NNC.objInds(3)) = dAreaFrac;
        end
        % normalized the objective functions
        for i = 1:numel(nlcon.NNC.objInds)
            if nlcon.NNC.objInds(i) > 0
                normalizedObj(nlcon.NNC.objInds(i)) ...
                    = (normalizedObj(nlcon.NNC.objInds(i)) ...
                        - nlcon.NNC.anchorPts(nlcon.NNC.objInds(i),nlcon.NNC.objInds(i)))...
                        /nlcon.NNC.normalization(nlcon.NNC.objInds(i));
                dNormalizedObj(1:designParams.nParams,nlcon.NNC.objInds(i)) ...
                    = dNormalizedObj(1:designParams.nParams,nlcon.NNC.objInds(i)) ...
                        /nlcon.NNC.normalization(nlcon.NNC.objInds(i));
            end
        end
        for i = 1:(nObjs-1)
            g(NNCind+i) = dot(nlcon.NNC.UtopiaLineVec(:,i),normalizedObj-nlcon.NNC.hyperPlanePt);
            gradg(:,NNCind+i) = dNormalizedObj*nlcon.NNC.UtopiaLineVec(:,i);
        end
    end
    

else
    
    i1 = 0;
    for i = 1:nPolygons
        XX = channels.vertexCoords(channels.polygons(i).connectivity,:)';
        i2 = i1 + channels.polygons(i).nVertices;
        sineAngles = polygon_sine_angles(XX(1,:),...
                                         XX(2,:)); 
        if (channels.polygons(i).isSideTriangle)
            g((i1+1):i2) = nlcon.sinMinSidePolyAngle - sineAngles; 
        else
            g((i1+1):i2) = nlcon.sinMinPolyAngle - sineAngles; 
        end
       
            
        i1 = i2;
        
        area = polygon_area(XX(1,:),...
                            XX(2,:));
        if (channels.polygons(i).isSideTriangle)
            g(nV+i) = nlcon.minSidePolyArea - area;
        else
            g(nV+i) = nlcon.minPolyArea - area;
        end

    end
   
   
    if PmaxConOn && TmaxConOn
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn,~,~,~,Tmax] = fem_n_sensitivity([],...
                                               mesh, ...
                                               channels, ...
                                               femParams, ...
                                               sensitivity, ...
                                               -3, ...
                                               false, ...
                                               false); % Tave is not scaled by area
        if isnan(Tn)
            warning('pressure and max temperature constraints returning NaN') 
            g(:) = NaN;            
            return;
        end 
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = 'PRESSURE';
        sensitivity.costFunction.objOpt.intDomainType = 1;
        delP = fem_n_sensitivity([], ...
                                mesh, ...
                                channels, ...
                                femParams,...
                                sensitivity, ...
                                -3, ...
                                false, ...
                                false);
        g(Pind) = delP/nlcon.maxP - 1.0;
        % For debugging
        %g(Tind) = Tn;
        %
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        G_Tmax2TnRatio = Tmax/Tn;
        %
    elseif PmaxConOn || PminConOn 
       delP = fem_n_sensitivity([], ...
                                mesh, ...
                                channels, ...
                                femParams,...
                                sensitivity, ...
                                -3, ...
                                false, ...
                                false);
        if PminConOn
             g(Pind) = 1.0 - delP/nlcon.minP;
         else
             g(Pind) = delP/nlcon.maxP - 1.0;
        end            
    elseif TmaxConOn
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn,~,~,~,Tmax] = fem_n_sensitivity([],...
                                            mesh, ...
                                            channels, ...
                                            femParams,...
                                            sensitivity, ...
                                            -3, ...
                                            false, ...
                                            false); 
        if isnan(Tn)
            warning('max temperature constraint returning NaN') 
            g(:) = NaN;         
            return;
        end  
        % For debugging
        %g(Tind) = Tn/nlcon.maxTmax - 1.0;
        %
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        G_Tmax2TnRatio = Tmax/Tn;
        %
    end
    if TnodalConOn 
        sensitivity.costFunction.objOpt.normp = 1;
        sensitivity.costFunction.objOpt.calcOneNorm = false;
        sensitivity.costFunction.type = '';
        sensitivity.costFunction.objOpt.intDomainType = 0;
        if evalNodeByNode
            originalMassIn = channels.massin;
            originalPairIndHeat = channels.pair_ind_heat;
            originalPtTemp = channels.pt_temp;
            gInd = TnodalInd;
            % Evaluate node by node since each constraint is subject to a
            % different flow rate and channel BC
            for i = 1:numel(nlcon.nodalTbounds)
                if ~isnan(nlcon.nodalTbounds(i))
                    if i > 2 
                        sensitivity.costFunction.objOpt.node = -2;
                    else
                        sensitivity.costFunction.objOpt.node = -1;
                    end
                    % Assume that all the inlets take the same value
                    if ~isnan(nlcon.massin(i))
                        % change the mass flow rate if applicable
                        channels.massin = nlcon.massin(i);
                    end
                    if ~isnan(nlcon.inBCs(i))
                        % change the temperature difference or the inlet
                        % temperature depending which one is specified
                        % originally
                        if size(channels.pair_ind_heat,2) == 3
                            channels.pair_ind_heat(:,3) = nlcon.inBCs(i);
                        elseif size(channels.pt_temp,2) == 2
                            channels.pt_temp(1,2) = nlcon.inBCs(i);
                        end
                    end
                    nodalT = fem_n_sensitivity(inf, ...
                                              mesh, ...
                                              channels, ...
                                              femParams,...
                                              sensitivity, ...
                                              -3, ...
                                              false, ...
                                              false);
                    if rem(i,2) == 1
                       % lower bound
                       g(gInd) = nlcon.nodalTbounds(i) - nodalT;
                    else
                       % upper bound
                       g(gInd) = nodalT - nlcon.nodalTbounds(i);
                    end                                                                              
                    gInd = gInd + 1;  
                    channels.massin = originalMassIn;
                    channels.pair_ind_heat = originalPairIndHeat;
                    channels.pt_temp = originalPtTemp;
                end
            end
        else
            inletNodeFlag = any(~isnan(nlcon.nodalTbounds(1:2)));
            outletNodeFlag = any(~isnan(nlcon.nodalTbounds(3:4)));
            if inletNodeFlag && outletNodeFlag 
                % Inlet and outlet nodes required
                sensitivity.costFunction.objOpt.node = -3;
            elseif inletNodeFlag
                % Inlet nodes only
                sensitivity.costFunction.objOpt.node = -1;
            elseif outletNodeFlag
                % Outlet nodes only
                sensitivity.costFunction.objOpt.node = -2;
            end
            nodalT = fem_n_sensitivity([], ...
                                      mesh, ...
                                      channels, ...
                                      femParams,...
                                      sensitivity, ...
                                      -3, ...
                                      false, ...
                                      false);
            if sensitivity.costFunction.objOpt.node == -3 && numel(nodalT) < 2
                % Handle the case where at least two nodes cannot be found
                nodalT = nan(2,1);
            end                                  
            gInd = TnodalInd;
            for i = 1:numel(nlcon.nodalTbounds)
               if ~isnan(nlcon.nodalTbounds(i))
                   if i > 2 && sensitivity.costFunction.objOpt.node == -3
                       nodalTind = 2;
                   else
                       nodalTind = 1;
                   end

                   if rem(i,2) == 1
                       % lower bound
                       g(gInd) = nlcon.nodalTbounds(i) - nodalT(nodalTind);
                   else
                       % upper bound
                       g(gInd) = nodalT(nodalTind) - nlcon.nodalTbounds(i);
                   end      
                   gInd = gInd + 1;
               end
            end
        end
    end
    % area or volume fraction constraint
    if AminConOn || AmaxConOn
         areaFrac = sum(channels.diams.*channels.length)/nlcon.area;
        if AminConOn
            g(AminInd) = nlcon.minAorVFrac - areaFrac;
        end
        if AmaxConOn
            g(AmaxInd) = areaFrac - nlcon.maxAorVFrac;
        end
    elseif VminConOn || VmaxConOn
        volFrac = channels.vol/channels.domainVol;
        if VminConOn
            g(VminInd) = nlcon.minAorVFrac - volFrac;
        end
        if VmaxConOn            
            g(VmaxInd) = volFrac - nlcon.maxAorVFrac;
        end
    end
    
    if isfield(nlcon,'NNC')
        if (nlcon.NNC.objInds(1) && nlcon.NNC.objInds(2))
            % T and P
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = true;
            sensitivity.costFunction.objOpt.intDomainType = 0;
            Tn = fem_n_sensitivity ([],...
                                    mesh, ...
                                    channels, ...
                                    femParams,...
                                    sensitivity,...
                                    -3, ...
                                    false, ...
                                    false); % Tave is not scaled by area

            if (isnan(Tn))
                warning('pressure and max temperature constraints returning NaN') 
                g(:) = NaN;
                return
            end                                         
            sensitivity.costFunction.objOpt.normp = 1;
            sensitivity.costFunction.objOpt.calcOneNorm = false;
            sensitivity.costFunction.type = nlcon.NNC.obj2type;
            sensitivity.costFunction.objOpt.intDomainType = nlcon.NNC.obj2intDomainType;
            delP = fem_n_sensitivity([], ...
                                      mesh, ...
                                      channels, ...
                                      femParams,...
                                      sensitivity, ...
                                      -3, ...
                                      false, ...
                                      false);
            % P
            normalizedObj(nlcon.NNC.objInds(2)) = delP;   
            % T       
            normalizedObj(nlcon.NNC.objInds(1)) = Tn;                                                    
        elseif nlcon.NNC.objInds(1)
            % T
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = false; 
            sensitivity.costFunction.objOpt.intDomainType = 0;
            Tn = fem_n_sensitivity([],...
                                   mesh, ...
                                   channels, ...
                                   femParams,...
                                   sensitivity, ...
                                   -3, ...
                                   false, ...
                                   false); 
            if (isnan(Tn))
                warning('max temperature constraint returning NaN') 
                g(:) = NaN;
                return;
            end     
            normalizedObj(nlcon.NNC.objInds(1)) = Tn;
        elseif nlcon.NNC.objInds(2)
            % P
            sensitivity.costFunction.type = 'PRESSURE';
            sensitivity.costFunction.objOpt.normp = 1;
            sensitivity.costFunction.objOpt.calcOneNorm = false;
            sensitivity.costFunction.type = nlcon.NNC.obj2type;
            sensitivity.costFunction.objOpt.intDomainType = nlcon.NNC.obj2intDomainType;
            delP = fem_n_sensitivity([], ...
                                      mesh, ...
                                      channels, ...
                                      femParams,...
                                      sensitivity, ...
                                      -3, ...
                                      false, ...
                                      false);
             normalizedObj(nlcon.NNC.objInds(2)) = delP;                      
        end
        if nlcon.NNC.objInds(3)
            areaFrac = sum(channels.diams.*channels.length)/nlcon.area;
            normalizedObj(nlcon.NNC.objInds(3)) = areaFrac;
        end
        % normalized the objective functions
        for i = 1:numel(nlcon.NNC.objInds)
            if nlcon.NNC.objInds(i) > 0
                normalizedObj(nlcon.NNC.objInds(i)) ...
                    = (normalizedObj(nlcon.NNC.objInds(i)) ...
                        - nlcon.NNC.anchorPts(nlcon.NNC.objInds(i),nlcon.NNC.objInds(i)))...
                        /nlcon.NNC.normalization(nlcon.NNC.objInds(i));
            end
        end
        for i = 1:(nObjs-1)
            g(NNCind+i) = dot(nlcon.NNC.UtopiaLineVec(:,i),normalizedObj-nlcon.NNC.hyperPlanePt);
        end
    end
end

if (isfield(nlcon,'areaScale'))
    fprintf('area constraint scaled \n')
    ind = (nV+1):(nV+nPolygons); 
    g(ind) = g(ind)./nlcon.areaScale; 
    if (nargout > 2)
        gradg(:,ind) = gradg(:,ind)./nlcon.areaScale;
    end
end


G_constraint.g = g;
end

