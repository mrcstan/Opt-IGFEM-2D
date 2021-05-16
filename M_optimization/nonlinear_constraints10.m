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
% VERY IMPORTANT TO REINITIALIZE THESE GLOBAL VARIABLES OUTSIDE OF THIS 
% FUNCTION AFTER EVERY SIMULATION !!!!!!!!!!!!!!                                    
global G_constraint
global G_Tmax2TnRatio
[channels,~,~] = update_channels(delParams, ...
                                 designParams, ...
                                 restrictedParams, ...
                                 iniChannels, ...
                                 'add');
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

PConOn = false;
TmaxConOn= false;
AminConOn = false;
AmaxConOn = false;
VminConOn = false;
VmaxConOn = false;

if (ischar(nlcon.type))
    splitStr = regexp(nlcon.type,',','split');
    for i = 1:numel(splitStr)
        switch splitStr{i}
            case 'P'
                if (~isfield(nlcon,'area'))
                    error('nlcon must have area field when nlcon.type=="P"')
                end
                if (~isfield(nlcon,'maxP'))
                    error('nlcon must have maxP field when nlcon.type=="P"')
                end
                if (~isfield(nlcon,'PScale'))
                    error('nlcon must have PScale field when nlcon.type=="P"')
                end
                PConOn = true; 
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
    


if (PConOn && TmaxConOn)
    if (isempty(nlcon.minP))
        PMinInd = [];
        PMaxInd = ng + 1;
        ng = ng + 2;
        Tind = ng;
    else
        PMinInd = ng + 1;
        PMaxInd = ng + 2;
        ng = ng + 3;
        Tind = ng;
    end
    sensitivity.costFunction.type = 'P_NORM';
    sensitivity.costFunction.objOpt.normp = 8;
    sensitivity.costFunction.objOpt.calcOneNorm = true;
elseif PConOn    
    if (isempty(nlcon.minP))
        PMinInd = [];
        PMaxInd = ng + 1;
        ng = ng + 1;
    else               
        PMinInd = ng + 1;
        PMaxInd = ng + 2;
        ng = ng + 2;
    end
    %sensitivity.costFunction.type = 'P_NORM';
    %sensitivity.costFunction.objOpt.normp = 1;
    %sensitivity.costFunction.objOpt.calcOneNorm = false;
    sensitivity.costFunction.type = 'PRESSURE';
    sensitivity.costFunction.objOpt.normp = 1;
    sensitivity.costFunction.objOpt.calcOneNorm = false;
    sensitivity.costFunction.area = nlcon.area;
elseif (TmaxConOn)
    ng = ng + 1;
    Tind = ng;
    sensitivity.costFunction.type = 'P_NORM';
    sensitivity.costFunction.objOpt.normp = 8;
    sensitivity.costFunction.objOpt.calcOneNorm = false;
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
    
    if (PConOn && TmaxConOn)
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn,D2Tn,Tave, D2Tave, Tmax] = fem_n_sensitivity ([],...
                                                          mesh, ...
                                                          channels, ...
                                                          femParams,...
                                                          sensitivity,...
                                                          2, ...
                                                          false, ...
                                                          false); % Tave is not scaled by area
                                                      
        if (isnan(Tn))
            warning('pressure and max temperature constraints returning NaN') 
            g(:) = NaN;
            gradg(:) = NaN;
            return;
        end                                         
        Tave = Tave/nlcon.area;
        [nu,DTnu] = kinematic_viscosity(channels.density,Tave,1);
        [pressure,~,~,~,~,~,DnuPressure] = network_pressure_mass_flow_rate(channels.contvty,...
                                                      channels.nurbs,...
                                                      channels.diams,...
                                                      nu,...
                                                      channels.inletEndPoint,...
                                                      channels.massin,...
                                                      channels.pressureOutletEndPoint,...
                                                      channels.pressureOut,...
                                                      channels.crossSection);
       
        g(PMaxInd) = (pressure(channels.inletEndPoint) ...
                    -pressure(channels.pressureOutletEndPoint) ...
                    - nlcon.maxP)/nlcon.PScale;
        D2Tave = D2Tave/nlcon.area;
        % scale the gradient of pressure channels.D2P, 
        % which is evaluated at a reference viscosity and does not include
        % the effect of the change in viscosity with temperature
        gradg(1:designParams.nParams,PMaxInd) ...
                = (DnuPressure(channels.inletEndPoint)*DTnu*D2Tave(1:designParams.nParams)...
                 + nu/channels.viscosity*channels.D2P(channels.inletEndPoint,:)')...
                 /nlcon.PScale;  
        if (~isempty(nlcon.minP))                                          
            g(PMinInd) = (nlcon.minP - pressure(channels.inletEndPoint) ...
                       + pressure(channels.pressureOutletEndPoint))/nlcon.PScale;
            gradg(1:designParams.nParams,PMinInd) = -gradg(1:designParams.nParams,PMaxInd);        
        end     
            
        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        gradg(:,Tind) = G_Tmax2TnRatio*D2Tn/nlcon.maxTmax;
        G_Tmax2TnRatio = Tmax/Tn;  
    elseif PConOn   
        %{
        [Tave, D2Tave] = fem_n_sensitivity([], ...
                                           mesh, ...
                                           channels, ...
                                           femParams,...
                                           sensitivity, ...
                                           2, ...
                                           false, ...
                                           false); % Tave is not scaled by area
        if (isnan(Tave))
            warning('pressure constraint returning NaN') 
            g(:) = NaN;
            gradg(:) = NaN;
            return;
        end
        Tave = Tave/nlcon.area;
        [nu,DTnu] = kinematic_viscosity(channels.density,Tave,1);
        [pressure,~,~,~,~,~,DnuPressure] = network_pressure_mass_flow_rate(channels.contvty,...
                                                      channels.nurbs,...
                                                      channels.diams,...
                                                      nu,...
                                                      channels.inletEndPoint,...
                                                      channels.massin,...
                                                      channels.pressureOutletEndPoint,...
                                                      channels.pressureOut,...
                                                      channels.crossSection);
         
         g(PMaxInd) = (pressure(channels.inletEndPoint) ...
                    -pressure(channels.pressureOutletEndPoint) ...
                    - nlcon.maxP)/nlcon.PScale;
        D2Tave = D2Tave/nlcon.area;
        % scale the gradient of pressure channels.D2P, 
        % which is evaluated at a reference viscosity and does not include
        % the effect of the change in viscosity with temperature
        gradg(1:designParams.nParams,PMaxInd) ...
                = (DnuPressure(channels.inletEndPoint)*DTnu*D2Tave(1:designParams.nParams)...
                 + nu/channels.viscosity*channels.D2P(channels.inletEndPoint,:)')...
                 /nlcon.PScale;
        if (~isempty(nlcon.minP))
            g(PMinInd) = (nlcon.minP - pressure(channels.inletEndPoint) ...
                       + pressure(channels.pressureOutletEndPoint))/nlcon.PScale;
            gradg(1:designParams.nParams,PMinInd) = -gradg(1:designParams.nParams,PMaxInd);       
        end
        %}
        [delP, DdelP] = fem_n_sensitivity([], ...
                                          mesh, ...
                                          channels, ...
                                          femParams,...
                                          sensitivity, ...
                                          2, ...
                                          false, ...
                                          false);
         g(PMaxInd) = (delP - nlcon.maxP)/nlcon.PScale;
         gradg(:,PMaxInd) = DdelP/nlcon.PScale;  
         if ~isempty(nlcon.minP)
            g(PMinInd) = (nlcon.minP - delP)/nlcon.PScale;
            gradg(1:designParams.nParams,PMinInd) = -gradg(1:designParams.nParams,PMaxInd);    
         end
    elseif TmaxConOn
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn, D2Tn, ~, ~, Tmax] = fem_n_sensitivity([],...
                                                   mesh, ...
                                                   channels, ...
                                                   femParams,...
                                                   sensitivity, ...
                                                   2, ...
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
    if isfield(nlcon,'NNC')
        if (nlcon.NNC.objInds(1) && nlcon.NNC.objInds(2))
            % T and P
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = true;
            [Tn,D2Tn,Tave, D2Tave] = fem_n_sensitivity ([],...
                                                        mesh, ...
                                                        channels, ...
                                                        femParams,...
                                                        sensitivity,...
                                                        2, ...
                                                        false, ...
                                                        false); % Tave is not scaled by area

            if (isnan(Tn))
                warning('pressure and max temperature constraints returning NaN') 
                g(:) = NaN;
                gradg(:) = NaN;
                return;
            end                                         
            Tave = Tave/nlcon.area;
            [nu,DTnu] = kinematic_viscosity(channels.density,Tave,1);
            [pressure,~,~,~,~,~,DnuPressure] ...
                    = network_pressure_mass_flow_rate(channels.contvty,...
                                                      channels.nurbs,...
                                                      channels.diams,...
                                                      nu,...
                                                      channels.inletEndPoint,...
                                                      channels.massin,...
                                                      channels.pressureOutletEndPoint,...
                                                      channels.pressureOut,...
                                                      channels.crossSection);
            % P
            normalizedObj(nlcon.NNC.objInds(2)) ...
                = pressure(channels.inletEndPoint) ...
                   -pressure(channels.pressureOutletEndPoint);       
            D2Tave = D2Tave/nlcon.area;
            % scale the gradient of pressure channels.D2P, 
            % which is evaluated at a reference viscosity and does not include
            % the effect of the change in viscosity with temperature
            dNormalizedObj(1:designParams.nParams,nlcon.NNC.objInds(2)) ...
                    = DnuPressure(channels.inletEndPoint)*DTnu*D2Tave(1:designParams.nParams)...
                     + nu/channels.viscosity*channels.D2P(channels.inletEndPoint,:)'; 
            % T       
            normalizedObj(nlcon.NNC.objInds(1)) = Tn;                                                
            dNormalizedObj(:,nlcon.NNC.objInds(1)) = D2Tn;      
        elseif nlcon.NNC.objInds(1)
            % T
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = false;   
            [Tn, D2Tn] = fem_n_sensitivity([],...
                                           mesh, ...
                                           channels, ...
                                           femParams,...
                                           sensitivity, ...
                                           2, ...
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
            sensitivity.costFunction.area = nlcon.area;
            [delP, DdelP] = fem_n_sensitivity([], ...
                                              mesh, ...
                                              channels, ...
                                              femParams,...
                                              sensitivity, ...
                                              2, ...
                                              false, ...
                                              false);
             normalizedObj(nlcon.NNC.objInds(2)) = delP;                      
             dNormalizedObj(1:designParams.nParams,nlcon.NNC.objInds(2)) = DdelP;
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
            dNormalizedObj(1:designParams.nParams,nlcon.NNC.objInds(3)) = dAreaFrac;
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
   
   
    if (PConOn && TmaxConOn)
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn,~,Tave,~,Tmax] = fem_n_sensitivity([],...
                                               mesh, ...
                                               channels, ...
                                               femParams, ...
                                               sensitivity, ...
                                               3, ...
                                               false, ...
                                               false); % Tave is not scaled by area
        if isnan(Tn)
            warning('pressure and max temperature constraints returning NaN') 
            g(:) = NaN;            
            return;
        end  
        Tave = Tave/nlcon.area;
        nu = kinematic_viscosity(channels.density,Tave,1);
        pressure = network_pressure_mass_flow_rate(channels.contvty,...
                                                   channels.nurbs,...
                                                   channels.diams,...
                                                   nu,...
                                                   channels.inletEndPoint,...
                                                   channels.massin,...
                                                   channels.pressureOutletEndPoint,...
                                                   channels.pressureOut,...
                                                   channels.crossSection);
        if (~isempty(nlcon.minP))
            g(PMinInd) = (nlcon.minP - pressure(channels.inletEndPoint) ...
                       + pressure(channels.pressureOutletEndPoint))/nlcon.PScale;
        end
        g(PMaxInd) = (pressure(channels.inletEndPoint) ...
                    -pressure(channels.pressureOutletEndPoint) ...
                    - nlcon.maxP)/nlcon.PScale;

        g(Tind) = G_Tmax2TnRatio*Tn/nlcon.maxTmax - 1.0;
        G_Tmax2TnRatio = Tmax/Tn; 
    elseif PConOn
        %{
        Tave = fem_n_sensitivity([],...
                                 mesh, ...
                                 channels, ...
                                 femParams,...
                                 sensitivity, ...
                                 3, ...
                                 false, ...
                                 false); % Tave is not scaled by area
        if (isnan(Tave))          
            warning('pressure constraint returning NaN')
            g(:) = NaN;
            return;
        end                        
        Tave = Tave/nlcon.area;
        nu = kinematic_viscosity(channels.density,Tave,1);
        pressure = network_pressure_mass_flow_rate(channels.contvty,...
                                                  channels.nurbs,...
                                                  channels.diams,...
                                                  nu,...
                                                  channels.inletEndPoint,...
                                                  channels.massin,...
                                                  channels.pressureOutletEndPoint,...
                                                  channels.pressureOut,...
                                                  channels.crossSection);                                     
       if (~isempty(nlcon.minP))
            g(PMinInd) = (nlcon.minP - pressure(channels.inletEndPoint) ...
                       + pressure(channels.pressureOutletEndPoint))/nlcon.PScale;
       end
       g(PMaxInd) = (pressure(channels.inletEndPoint) ...
                    -pressure(channels.pressureOutletEndPoint) ...
                    - nlcon.maxP)/nlcon.PScale;
       %}
       delP = fem_n_sensitivity([], ...
                                mesh, ...
                                channels, ...
                                femParams,...
                                sensitivity, ...
                                3, ...
                                false, ...
                                false);
        if ~isempty(nlcon.minP)
            g(PMinInd) = (nlcon.minP - delP)/nlcon.PScale;
        end                       
        g(PMaxInd) = (delP - nlcon.maxP)/nlcon.PScale;             
    elseif TmaxConOn
        % NOTE: Tn in this case is the n-norm of temperature not the n-mean
        [Tn,~,~,~,Tmax] = fem_n_sensitivity([],...
                                            mesh, ...
                                            channels, ...
                                            femParams,...
                                            sensitivity, ...
                                            3, ...
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
    
    if isfield(nlcon,'NNC')
        if (nlcon.NNC.objInds(1) && nlcon.NNC.objInds(2))
            % T and P
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = true;
            [Tn,~,Tave] = fem_n_sensitivity ([],...
                                             mesh, ...
                                             channels, ...
                                             femParams,...
                                             sensitivity,...
                                             3, ...
                                             false, ...
                                             false); % Tave is not scaled by area

            if (isnan(Tn))
                warning('pressure and max temperature constraints returning NaN') 
                g(:) = NaN;
                return
            end                                         
            Tave = Tave/nlcon.area;
            nu = kinematic_viscosity(channels.density,Tave,1);
            pressure = network_pressure_mass_flow_rate(channels.contvty,...
                                                      channels.nurbs,...
                                                      channels.diams,...
                                                      nu,...
                                                      channels.inletEndPoint,...
                                                      channels.massin,...
                                                      channels.pressureOutletEndPoint,...
                                                      channels.pressureOut,...
                                                      channels.crossSection);
            % P
            normalizedObj(nlcon.NNC.objInds(2)) ...
                = pressure(channels.inletEndPoint) ...
                   -pressure(channels.pressureOutletEndPoint);   
            % T       
            normalizedObj(nlcon.NNC.objInds(1)) = Tn;                                                    
        elseif nlcon.NNC.objInds(1)
            % T
            sensitivity.costFunction.type = 'P_NORM';
            sensitivity.costFunction.objOpt.normp = 8;
            sensitivity.costFunction.objOpt.calcOneNorm = false;   
            Tn = fem_n_sensitivity([],...
                                   mesh, ...
                                   channels, ...
                                   femParams,...
                                   sensitivity, ...
                                   3, ...
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
            sensitivity.costFunction.area = nlcon.area;
            delP = fem_n_sensitivity([], ...
                                      mesh, ...
                                      channels, ...
                                      femParams,...
                                      sensitivity, ...
                                      3, ...
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
G_constraint = g;

end

