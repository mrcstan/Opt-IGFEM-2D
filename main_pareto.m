generate_pareto_front('parallel2_start.channel',...% channelFile
                      'parallelTwo_NE.polygon',... % polygonFile
                      'parallelTwo_NE_a.rand',... % boundsFile/sampleFile
                       5,... % job num (if > 0 && nSimulations > 0)/sample number (if < 0 && nSimulations == 1)
                       3,... % nSimulations
                       true,... % randomizeFirst
                      'test',... % directory
                      'insulated_composite',...% filePrefix
                      'battery_composite_panel.in',... % input file
                      'T,Ppa',... % costFuncType
                       8, ... % costFuncNormp
                       [25.49, 246.5;
                        46994, 261], ... % anchor points
                       22, ... % number of points in first objective direction on hyperplane (includes both end points)
                       21, ... % 1 <= hyperplane point number <= nHyperPlanePts1-2
                       [],... % nlconType
                       [], ... % nlconMinP
                       0, ... % nlconMaxP
                       0, ... % nlconMaxTmax
                       0.018, ... % nlconMinAorVFrac
                       0.2, ... % nlconMaxAorVFrac
                       false); % mergeTriangles

% [261, 46994;
% 246.5, 25.49], ... % anchor points                   

                   