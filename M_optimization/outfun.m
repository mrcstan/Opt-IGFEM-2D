function stop = outfun(x,optimValues,state)
    global G_simDirectory
    global G_history
    global G_filePrefix
    global G_constraint
    global G_out
    global G_Tmax2TnRatio
    stop = false;
    
   
    switch state
         %case 'init'    
          
         case 'iter'
            fprintf('Output results at iteration %i\n',optimValues.iteration)
            % outputTimer = tic;
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
            % NOTE: when minimizing the Tpnorm, the fval in G_history is
            %       used to be Tpmean. Now (Oct 5, 2015), it just Tpmean
            G_history.fval = [G_history.fval, optimValues.fval*G_out.objScale+G_out.objOffset];
            G_history.iter = [G_history.iter, optimValues.iteration];
            G_history.x = [G_history.x, x];
            %G_history.cstrViolation = [G_history.cstrViolation; optimValues.constrviolation];
            G_history.cstr = [G_history.cstr,G_constraint.g];
            %G_history.inPressure = [G_history.inPressure,G_out.inPressure];
            G_history.firstorderopt = [G_history.firstorderopt,optimValues.firstorderopt]; 
            G_history.Tmax2TnRatio = [G_history.Tmax2TnRatio,G_Tmax2TnRatio];
            
            filePrefix = [G_simDirectory,G_filePrefix];
            
            % output vtk file of this iteration
            strIter = num2str(optimValues.iteration);
            %{
            load([filePrefixLast,'_last_mesh'])
            load([filePrefixLast,'_last_soln'])
            %}
  
            matlab2vtk_scalar([filePrefix,'_',strIter,'.vtk'], ...
                              'T', ...
                              G_out.elem, ...
                              G_out.node, ...
                              G_out.UUR_update - G_out.Toffset);
            % NOTE: Tave and Tmax are converted to the temperature unit
            %       specified in the input file
            Tave = average_temp(G_out.elem,G_out.node.coords,G_out.UUR_update) - G_out.Toffset;
            Tmax = max(G_out.UUR_update) - G_out.Toffset;
            G_history.Tave = [G_history.Tave,Tave];
            G_history.Tmax = [G_history.Tmax,Tmax];
            % The nodal temperature is converted to the temperature unit
            % specified in the input file
            if ~isempty(regexpi(G_out.costFuncType,'NODAL_T'))
                G_history.fval(end) = G_history.fval(end) - G_out.Toffset;
            end
            
            if ~strcmpi(G_out.costFuncType,'P_NORM')
                 G_history.Tpnorm = [G_history.Tpnorm,...
                                 field_p_norm(G_out.elem,G_out.node.coords,G_out.UUR,...
                                              G_out.normp,Tmax,G_out.gauss)];
            end
            % only take the first node even if there are multiple nodes 
            Tin = G_out.UUR_update(G_out.inletNodes(1)) - G_out.Toffset;                    
            Tout = G_out.UUR_update(G_out.outletNodes(1)) - G_out.Toffset;
            G_history.Tin(end+1)= Tin;
            G_history.Tout(end+1)= Tout;
            
            TchanStats = channel_T_stats(G_out.UUR_update,G_out.channels.origNodes, ...
                                         G_out.channels.enrichNodes) - G_out.Toffset;
            G_history.TaveChan = [G_history.TaveChan,TchanStats(1)];
            G_history.TminChan = [G_history.TminChan,TchanStats(2)];
            G_history.TmaxChan = [G_history.TmaxChan,TchanStats(3)];
            if G_out.nuFlag == 0
                % evaluate viscosity at panel average temperature
                nuTave = Tave; 
            elseif G_out.nuFlag == 1
                % evaluate viscosity at network average temperature
                nuTave = TchanStats(1);
            end
            nu =  kinematic_viscosity(G_out.channels.viscosity, ...
                                      G_out.channels.density,nuTave, ...
                                      G_out.channels.viscosityModel,...
                                      G_out.Tunit);
            G_history.inPressure = [G_history.inPressure,...
                                    G_out.inPressure*nu/G_out.channels.viscosity];
            if isempty(G_out.channels.massin)
                G_out.channels.massin = max(G_out.channels.mcf)/G_out.channels.heatCapacity;
                G_history.inFlow = [G_history.inFlow, ...
                                    G_out.channels.massin...
                                    /G_out.channels.density*6e7]; % ml/min
            end
                       
            G_history.volFrac = [G_history.volFrac,G_out.channels.vol/G_out.channels.domainVol];
            G_history.areaFrac = [G_history.areaFrac,...
                                  sum(G_out.channels.diams.*G_out.channels.length)...
                                  /G_out.channels.domainArea];
            save([filePrefix,'_history'],'G_history');
            
            % output history figure
            %fig = figure('visible','off');
            fig = figure;
            %plot(optimValues.iteration,optimValues.fval,'ro','markerfacecolor','r');
            switch G_out.costFuncType
                case {'P_NORM','P_NORM_CHANNEL','NODAL_T_OUT','NODAL_T_IN'}
                    plot_opt_history(G_history.iter, ...
                                     G_history.fval, ...
                                     G_history.Tmax, ...
                                     'doubleY',false, ...
                                     'ylabel1','T');
                case 'P_NORM_CHANNEL_W_OFFSET'
                    % NOTE: fval is the p-mean temperature
                    plot_opt_history(G_history.iter, ...
                                     G_history.fval, ...
                                     G_history.TminChan, ...
                                     'doubleY',true, ...
                                     'ylabel1','Obj','ylabel2','T');
                case 'VARIANCE'
                    G_history.SD = [G_history.SD;sqrt(G_history.fval(end))];
                    plot_opt_history(G_history.iter, ...
                                     G_history.SD, ...
                                     G_history.Tmax, ...
                                     'doubleY',true, ...
                                     'ylabel1','\sigma(T)','ylabel2','T_{max}');              
                case 'PRESSURE' 
                    plot_opt_history(G_history.iter, ...
                                     G_history.inPressure, ...
                                     G_history.Tmax, ...
                                     'doubleY',true, ...
                                     'ylabel1','\Delta P','ylabel2','T_{max}');
                otherwise
                    plot_opt_history(G_history.iter, ...
                                     G_history.fval, ...
                                     G_history.Tmax, ...
                                     'doubleY',true, ...
                                     'ylabel1','obj','ylabel2','T_{max}');      
            end
            figName = [filePrefix,'_history'];
            %set(fig, 'PaperPositionMode', 'auto');
            %print('-djpeg ', '-r300', figName);
            %set(gcf,'units','normalized','outerposition',[0 0 1 1])
            saveas(fig, figName, 'jpg');
            
            
            %fig = figure('visible','off');
            fig = figure;
            plot(G_history.iter, G_history.inPressure/1000.0,...
                'linestyle','-','linewidth',2,...
                'marker','p','color','r', ...
                'markerfacecolor','r','markersize',12)
            xlabel('Iteration','fontsize',20)
            ylabel('\Delta P/1000','fontsize',20)
            set(gca,'fontsize',20)
            figName = [filePrefix,'_pressure_history'];
            saveas(fig, figName, 'jpg');
            
            if ~isempty(G_history.inFlow)
                 %fig = figure('visible','off');
                 fig = figure;
                 plot(G_history.iter, G_history.inFlow,...
                    'linestyle','-','linewidth',2,...
                    'marker','p','color','r', ...
                    'markerfacecolor','r','markersize',12)
                 xlabel('Iteration','fontsize',20)
                 ylabel('Flow rate (ml/min)','fontsize',20)
                 set(gca,'fontsize',20)
                 figName = [filePrefix,'_flow_rate_history'];
                 saveas(fig, figName, 'jpg');
            end
            
            %{
            %fig = figure('visible','off');
            fig = figure;
            plot(G_history.iter, G_history.areaFrac,...
                'linestyle','-','linewidth',2,...
                'marker','p','color','r', ...
                'markerfacecolor','r','markersize',12)
            xlabel('Iteration','fontsize',20)
            ylabel('Area fraction','fontsize',20)
            set(gca,'fontsize',20)
            figName = [filePrefix,'_areaFrac_history'];
            saveas(fig, figName, 'jpg');
            %}
            
            % load([filePrefixLast,'_last_channel'])     
            write_channel_file([filePrefix,'_',strIter,'.channel'],...
                                G_out.channels,...
                                G_out.Toffset, ...
                                G_out.designParams,'w',[]);
            %{
            if (isfield(G_out.channels,'polygons'))
                figName = [filePrefix,'_channel_polygons_',strIter];
                %fig = figure('visible','off');
                fig = figure;
                plot_channel_polygons(G_out.channels.polygons,...
                                      G_out.channels.vertexCoords,...
                                      G_out.designParams.vertices2params,...
                                      false)
                saveas(fig, figName, 'jpg');
            end
            %}
            % output channel configuration
            %fig = figure('visible','off');
            fig = figure;
            meshCurveOption.showMesh = false;
            plot_mesh_curve(G_out.node.coords, ...
                            G_out.elem.elem_node, ...
                            [], ...
                            G_out.channels,meshCurveOption)
            %set(fig,'ResizeFcn','set(gcf,''visible'',''on'')'); 
            set(fig,'ResizeFcn')
            switch G_out.costFuncType
                case 'P_NORM'
                    title(['i=',strIter,', T_{ave}=',num2str(Tave,4),...
                           ', T_{max}=',num2str(Tmax,4)])
                case 'P_NORM_CHANNEL'
                    title(['i=',strIter,', T_{chan,ave}=',num2str(TchanStats(1),4),...
                           ', T_{chan,max}=',num2str(TchanStats(3),4)])
                case 'P_NORM_CHANNEL_W_OFFSET'
                    title(['i=',strIter,', T_{chan,ave}=',num2str(TchanStats(1),4),...
                           ', T_{chan,min}=',num2str(TchanStats(2),4)])
                case 'VARIANCE'
                    title(['i=',strIter,...
                           ', \sigma(T)=',num2str(sqrt(G_history.fval(end)),4),...
                           ', T_{max}=',num2str(Tmax,4)])
                case 'PRESSURE'
                    title(['i=',strIter,', \DeltaP =',...
                           num2str(G_history.inPressure(end)/101325.0,4),...
                           ', T_{chan,ave}=',num2str(TchanStats(1),4),...
                           ', T_{chan,min}=',num2str(TchanStats(2),4)])
                case {'NODAL_T_IN','NODAL_T_OUT'}
                    title(['i=',strIter,', T_{in}=',num2str(Tin,4),...
                           ', T_{out}=',num2str(Tout,4)])        
                otherwise
                    title(['i=',strIter,', \DeltaP/1e3=',...
                           num2str(G_history.inPressure(end)/1000.0,3),...
                           ', T_{max}=',num2str(Tmax,3),...
                           ', A_f=',num2str(G_history.areaFrac(end),3)])
                       
                
            end
            figName = [filePrefix,'_channel_',strIter];
            %set(fig, 'PaperPositionMode', 'auto');
            %print('-djpeg ', '-r300', figName);
            saveas(fig, figName, 'jpg');
            saveas(fig, figName, 'fig');
            %toc(outputTimer);
 
         case 'done'
             hold off
         otherwise
     end           
end   