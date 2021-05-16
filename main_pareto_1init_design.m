generate_pareto_front('parallelTwoNEoptimal.channel',...% channelFile
                      'parallelTwo_NE.polygon',... % polygonFile
                      'test2',... % directory
                      'insulated_composite',...% filePrefix
                      'battery_composite_panel.in',... % input file
                      'T,P',... % costFuncType
                       8, ... % costFuncNormp
                       [25.49,246.5;
                        46994,261], ... % anchor points
                       6, ... % number of points in first objective direction on hyperplane
                       [],... % nlconType
                       [], ... % nlconMinP
                       0, ... % nlconMaxP
                       0, ... % nlconMaxTmax
                       0.018, ... % nlconMinAorVFrac
                       0.2, ... % nlconMaxAorVFrac
                       false); % mergeTriangles
