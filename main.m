%diary('debug.out')
optimize_channels('parallel2_start.channel',...% channelFile
                  'parallelTwo_NE.polygon',... % polygonFile
                  'parallelTwo_NE_a.rand',... % boundsFile/sampleFile
                   -5,... % job num (if > 0 && nSimulations > 0)/sample number (if < 0 && nSimulations == 1)
                  'test',... % directory
                  'insulated_composite',...% filePrefix
                  'battery_composite_panel.in',... % input file
                  'PRESSURE_AVE_PANEL',... % costFuncType
                   8, ... % costFuncNormp
                   [1,1,1]/3.0, ... % costFuncWeights
                   [30,1e5,0.1], ... % costFuncScales
                   [],... % nlcon.Type
                   1000, ... % nlconMinP
                   20000, ... % nlconMaxP
                   44, ... % nlconMaxTmax
                   0.018, ... % nlconMinAorVFrac
                   0.04, ... % nlconMaxAorVFrac
                   [nan,80,nan,nan], ... % nlconNodalTbounds [lower Tin, upper Tin, lower Tout, upper Tout]
                   [nan,0.001418,nan,nan], ... % nlconMassIn [massin or nan, , ,]
                   [nan,80,nan,nan], ...% nlconInBCs [Tin or Q or nan, , ,]
                   1,... % nSimulations
                   true,... % randomizeFirst
                   false); % mergeTriangles
%diary off
