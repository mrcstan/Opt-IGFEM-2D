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
                = nonlinear_constraints_blocked_channels(delParams, ...
                                                         blockedSets, ...
                                                         designParams,...
                                                         restrictedParams, ...
                                                         iniChannels, ...
                                                         nlcon, ...
                                                         mesh, ...
                                                         femParams, ...
                                                         sensitivity)
% VERY IMPORTANT TO REINITIALIZE THESE GLOBAL VARIABLES OUTSIDE OF THIS 
% FUNCTION AFTER EVERY SIMULATION !!!!!!!!!!!!!! 
global G_constraint
global G_Tmax2TnRatio
% the last parameter is the extra parameter introduced for the minimax
% problem
[channels,~,~] = update_channels(delParams(1:end-1), ...
                                 designParams, ...
                                 restrictedParams, ...
                                 iniChannels, ...
                                 'add', ...
                                 (nargout > 2));
if (~isfield(channels,'polygons'))
    channels.polygons = [];
    nPolygons = 0;
    nV = 0;
    ng = 0;
    warning('no polygons specified for nonlinear constraint')
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

if (ischar(nlcon.type))
    splitStr = regexp(nlcon.type,',','split');
    sensitivityPTAV = sensitivity;
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
            otherwise
                error('unknown nonlinear constraint type')
        end
    end
end

%{
if strcmpi(sensitivity.costFunction.type,'AREA')
    % update Gout in fem_n_sensitivity when it is called from this function
    % only when the objective function is area, since the calculation the
    % area objective function does not go through any FEM calculation
    updateGout = true;
else
    costTypeSplit = regexp(sensitivity.costFunction.type,',','split');
    if ~isempty(cell2mat(strfind(costTypeSplit,'A')))
        updateGout = true;
    else
        updateGout = false;
    end
end
%}

if (PmaxConOn && TmaxConOn)
    Pind = ng + 1;
    ng = ng + 2;
    Tind = ng;
    sensitivityPTAV.costFunction.type = 'P_NORM';
    sensitivityPTAV.costFunction.objOpt.normp = 8;
    sensitivityPTAV.costFunction.objOpt.calcOneNorm = false;
    sensitivityPTAV.costFunction.objOpt.intDomainType = 0;
elseif PmaxConOn || PminConOn  
    Pind = ng + 1;
    ng = ng + 1;
    sensitivityPTAV.costFunction.type = 'PRESSURE';
    sensitivityPTAV.costFunction.objOpt.normp = 1;
    sensitivityPTAV.costFunction.objOpt.calcOneNorm = false;
    sensitivityPTAV.costFunction.objOpt.intDomainType = 1;
elseif TmaxConOn
    ng = ng + 1;
    Tind = ng;
    sensitivityPTAV.costFunction.type = 'P_NORM';
    sensitivityPTAV.costFunction.objOpt.normp = 8;
    sensitivityPTAV.costFunction.objOpt.calcOneNorm = false;
    sensitivityPTAV.costFunction.objOpt.intDomainType = 0;
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

% nonlinear constraints for minimax optimization
nMaxObjs = numel(blockedSets);
maxObjInd = ng + 1;
ng = ng + nMaxObjs;


g = zeros(ng, 1);
geq = [];

% angle constraints
if (nargout > 2)
    totParams = designParams.nParams+restrictedParams.nParams+1;
    gradg = zeros(totParams,ng); % last row corresponds to gradient wrt minimax parameter
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
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn,D2Tn,~, ~, Tmax] = fem_n_sensitivity ([],...
                                                  mesh, ...
                                                  channels, ...
                                                  femParams,...
                                                  sensitivityPTAV,...
                                                  -2, ...
                                                  false, ... % don't scale objective
                                                  false); % Tave is not scaled by area
                                                      
        if (isnan(Tn))
            warning('pressure and max temperature constraints returning NaN') 
            g(:) = NaN;
            gradg(:) = NaN;
            return;
        end                                         
        sensitivityPTAV.costFunction.objOpt.normp = 1;
        sensitivityPTAV.costFunction.objOpt.calcOneNorm = false;
        sensitivityPTAV.costFunction.type = 'PRESSURE';
        sensitivityPTAV.costFunction.objOpt.intDomainType = 1;
        [delP, DdelP] = fem_n_sensitivity([], ...
                                          mesh, ...
                                          channels, ...
                                          femParams,...
                                          sensitivityPTAV, ...
                                          -2, ...
                                          false, ...
                                          false);
       
        g(Pind) = delP/nlcon.maxP-1.0;
        gradg(:,Pind) = DdelP/nlcon.maxP;
                 
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        gradg(:,Tind) = G_Tmax2TnRatio*D2Tn/nlcon.maxTmax;
        G_Tmax2TnRatio = Tmax/Tn;  
    elseif PmaxConOn || PminConOn   
        [delP, DdelP] = fem_n_sensitivity([], ...
                                          mesh, ...
                                          channels, ...
                                          femParams,...
                                          sensitivityPTAV, ...
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
                                                   sensitivityPTAV, ...
                                                   -2, ...
                                                   false, ...
                                                   false); 
        if (isnan(Tn))
            warning('max temperature constraint returning NaN') 
            g(:) = NaN;
            gradg(:) = NaN;
            return;
        end     
        %g(Tind) = Tn/nlcon.maxTmax - 1.0;
        %gradg(:,Tind) = D2Tn/nlcon.maxTmax;
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        gradg(:,Tind) = G_Tmax2TnRatio*D2Tn/nlcon.maxTmax;        
        G_Tmax2TnRatio = Tmax/Tn;
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
   
    [maxObjs,gradMaxObjs] = fem_n_sensitivity_blocked_channels([], ...
                                                               blockedSets, ...
                                                               mesh, ...
                                                               channels, ...
                                                               femParams, ...
                                                               sensitivity, ...
                                                               false,true,true);
    if isfield(sensitivity,'objScale')
        g(maxObjInd:end) = (maxObjs - delParams(end))/sensitivity.objScale;
        gradg(:,maxObjInd:end) = [gradMaxObjs;-ones(1,nMaxObjs)]/sensitivity.objScale;
    else
        g(maxObjInd:end) = maxObjs - delParams(end);
        gradg(:,maxObjInd:end) = [gradMaxObjs;-ones(1,nMaxObjs)];
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
                                               sensitivityPTAV, ...
                                               -3, ...
                                               false, ...
                                               false); % Tave is not scaled by area
        if isnan(Tn)
            warning('pressure and max temperature constraints returning NaN') 
            g(:) = NaN;            
            return;
        end  
        sensitivityPTAV.costFunction.objOpt.normp = 1;
        sensitivityPTAV.costFunction.objOpt.calcOneNorm = false;
        sensitivityPTAV.costFunction.type = 'PRESSURE';
        sensitivityPTAV.costFunction.objOpt.intDomainType = 1;
        delP = fem_n_sensitivity([], ...
                                  mesh, ...
                                  channels, ...
                                  femParams,...
                                  sensitivityPTAV, ...
                                  -3, ...
                                  false, ...
                                  false);
        g(Pind) = delP/nlcon.maxP-1.0;
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        G_Tmax2TnRatio = Tmax/Tn; 
    elseif PmaxConOn || PminConOn 
       delP = fem_n_sensitivity([], ...
                                mesh, ...
                                channels, ...
                                femParams,...
                                sensitivityPTAV, ...
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
                                            sensitivityPTAV, ...
                                            -3, ...
                                            false, ...
                                            false); 
        if isnan(Tn)
            warning('max temperature constraint returning NaN') 
            g(:) = NaN;         
            return;
        end  
        %g(Tind) = Tn/nlcon.maxTmax - 1.0;
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        G_Tmax2TnRatio = Tmax/Tn;       
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
    
   maxObjs = fem_n_sensitivity_blocked_channels([], ...
                                                blockedSets, ...
                                                mesh, ...
                                                channels, ...
                                                femParams, ...
                                                sensitivity, ...
                                                false,true,true);
    if isfield(sensitivity,'objScale')
        g(maxObjInd:end) = (maxObjs - delParams(end))/sensitivity.objScale;
    else
        g(maxObjInd:end) = maxObjs - delParams(end);
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
G_constraint = g;
end

