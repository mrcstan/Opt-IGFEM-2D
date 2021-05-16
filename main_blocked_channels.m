%diary('debug.out')
optimize_blocked_channels('3x3localRefLp075deg6_p68xp53mm.channel',...% channelFile
                          '3x3Lp075NE.polygon',... % polygonFile
                          '3x3Lp075_a.rand',... % boundsFile_sampleFile
                          '3x3deg6nBlk2.blk',... % file indicating sets of channels that are blocked
                           -5,... % job num (if > 0 && nSimulations > 0)/sample number (if < 0 && nSimulations == 1)
                          'test1',... % directory
                          'insulated_composite',...% filePrefix
                          'square_PDMS.in',...% input file
                          'P_NORM',... % costFuncType
                           8, ... % costFuncNormp
                           [],... % nlconType
                           1000, ... % nlconMinP
                           20000, ... % nlconMaxP
                           44, ... % nlconMaxTmax
                           0.018, ... % nlconMinAorVFrac
                           0.2, ... % nlconMaxAorVFrac
                           [nan,nan,nan,nan], ... % nlconNodalTbounds [lower Tin, upper Tin, lower Tout, upper Tout]
                           1,... % nSimulations
                           true,... % randomizeFirst
                           false); % mergeTriangles
                           
%diary off  

