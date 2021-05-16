close all
clear all
path(path,'../M_channels')
path(path,'../M_geom_toolbox')
path(path,'../channelFiles')
path(path,'../../NURBS/nurbs_toolbox')
path(path,'../../Opt-IGFEM-Curves-2D/M_optimization')
channelFile = 'parallel2_start.channel';
%readOptions.boundsFile = 'parallel2_start_bounds_w_diams.channel';
readOptions.boundsFile = [];
readOptions.lhsFile = 'parallelTwo_NE_a.lhs';
readOptions.sampleNum = 16;
readOptions.margeTriangles = false;
readOptions.polygonFile = 'parallelTwo_NE.polygon';
readOptions.nlconfun = @nonlinear_constraints;
readOptions.nlcon.minPolyArea = 0.001*0.15*0.2;
readOptions.nlcon.sinMinPolyAngle = sin(5*pi/180);
readOptions.nlcon.minSidePolyArea = 0.001*0.15*0.2;
readOptions.nlcon.sinMinSidePolyAngle = sin(1*pi/180);
readOptions.nlcon.minDistSq = (0.001*0.15)^2;
readOptions.nlcon.areaScale = 0.15*0.2;
readOptions.nlcon.distSqScale = 0.15^2;

[channels,~,designParams,restrictedParams] = preprocess_channels(channelFile,readOptions);


del = 1e-7;
method = 'central';

func = @(x) nonlinear_constraints(x,designParams,restrictedParams,...
                                  channels,readOptions.nlcon);
nParams = designParams.nParams+restrictedParams.nParams;
xo = zeros(nParams,1);
[gOld,~,gradg] = func(xo);
FD_g = nan(nParams,numel(gOld));
for i = 1:nParams
    if (strcmpi(method,'forward'))
        xdel = xo;
        xdel(i) = xdel(i) + del; 
        gNew1 = func(xdel);
        FD_g(i,:) = (gNew1 - gOld)'/del;
        
    elseif (strcmpi(method,'central'))
        xdel = xo;
        xdel(i) = xdel(i) + del; 
        gNew2 = func(xdel);

        xdel = xo;
        xdel(i) = xdel(i) - del; 
        gNew1 = func(xdel);

        FD_g(i,:) = 0.5*(gNew2 - gNew1)'/del;

    else
        error('unrecognized method')
    end
end
abs_gradg_diff = abs(FD_g - gradg);
fprintf('maximum abs diff = %g \n',max(abs_gradg_diff(:)));

rel_diff = abs(1 - FD_g./gradg);
fprintf('maximum rel diff = %g \n',max(rel_diff(:)));
