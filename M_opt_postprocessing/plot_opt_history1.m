%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/5/2014
%%% Last modified date: 11/5/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   X
%   Y1
%   Y2
function [] = plot_opt_history(X,Y1,Y2)
ftsize = 16;
[ax,p1,p2] = plotyy(X, Y1, ...
                    X, Y2);
set(p1,'Marker','o','linewidth',2,'linestyle','-')
set(p2,'Marker','^','linewidth',2,'linestyle','--')
xlabel(ax(1),'iteration',  'FontSize', ftsize)
ylabel(ax(1),'objective',  'FontSize',  ftsize)
xlabel(ax(2),'iteration',  'FontSize', ftsize)
ylabel(ax(2),'constraint',  'FontSize',  ftsize)
set(ax(1),'fontsize',ftsize)
set(ax(2),'fontsize',ftsize)

end