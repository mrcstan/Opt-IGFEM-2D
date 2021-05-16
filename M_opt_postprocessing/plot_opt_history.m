%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 11/5/2014
%%% Last modified date: 11/5/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   1) X
%   2) Y1
%   3) Y2
%   4) 6) 8) ... 'property name'
%   5) 7) 9) ... property value
    
function [ax,p1,p2] = plot_opt_history(varargin)
if (nargin < 3)
    error('must have more than 2 arguments')
end

if (rem(nargin-3,2) ~= 0)
    error('each property name must be followed by a property value')
end

xlab = 'Iteration';
ylab1 = 'Objective';
ylab2 = 'Constraint';
ftsize = 20;
doubleY = true;

lnwidth = 2;
marker1 = 'o';
marker2 = '^';
style1 = '-';
style2 = '--';
color1 = 'b'; 
color2 = 'k';
for i = 4:2:nargin-1
    if (~ischar(varargin{i}))
        error('property name must be a string')
    end
    
    if (strcmpi(varargin{i},'xlabel'))
        if (~ischar(varargin{i+1}))
            error('xlabel must be a string')
        end
        xlab = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'ylabel1'))
        if (~ischar(varargin{i+1}))
            error('ylab1 must be a string')
        end
        ylab1 = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'ylabel2'))
        if (~ischar(varargin{i+1}))
            error('ylab2 must be a string')
        end
        ylab2 = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'marker1'))
        if (~ischar(varargin{i+1}))
            error('marker1 must be a string')
        end
        marker1 = varargin{i+1};
    end
    

    
    if (strcmpi(varargin{i},'marker2'))
        if (~ischar(varargin{i+1}))
            error('marker2 must be a string')
        end
        marker2 = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'style1'))
        if (~ischar(varargin{i+1}))
            error('style1 must be a string')
        end
        style1 = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'style2'))
        if (~ischar(varargin{i+1}))
            error('style2 must be a string')
        end
        style2 = varargin{i+1};
    end
    
     if (strcmpi(varargin{i},'color1'))
        if (~ischar(varargin{i+1}))
            error('color1 must be a string')
        end
        color1 = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'color2'))
        if (~ischar(varargin{i+1}))
            error('color2 must be a string')
        end
        color2 = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'fontsize'))
        if (~isnumeric(varargin{i+1}))
            error('fontsize must be a number')
        end
        ftsize = varargin{i+1};
    end
    
    if (strcmpi(varargin{i},'doubleY'))
        if (ischar(varargin{i+1}))
            error('doubleY must be a logical value or a number')
        end
        if (varargin{i+1})
            doubleY = true;
        else
            doubleY = false;
        end
      
    end
end

if (doubleY)
    [ax,p1,p2] = plotyy(varargin{1}, varargin{2}, ...
                        varargin{1}, varargin{3});
    set(p1,'Marker',marker1,'linewidth',lnwidth,'linestyle',style1)
    set(p2,'Marker',marker2,'linewidth',lnwidth,'linestyle',style2)
    xlabel(ax(1),xlab,  'FontSize', ftsize)
    ylabel(ax(1),ylab1,  'FontSize',  ftsize)
    xlabel(ax(2),xlab,  'FontSize', ftsize)
    ylabel(ax(2),ylab2,  'FontSize',  ftsize)
    set(ax(1),'fontsize',ftsize)
    set(ax(2),'fontsize',ftsize)
    set(ax(1),'units','normalized','position',[0.12,0.15,0.75,0.75])
    set(ax(2),'units','normalized','position',[0.12,0.15,0.75,0.75])
    %{
    xlim(ax(1),[0,60])
    ylim(ax(1),[40,80])
    xlim(ax(2),[0,60])
    ylim(ax(2),[-2e-3,0])
    %}
else
     ax = [];
     p1 = plot(varargin{1}, varargin{2}, 'Marker',marker1,'linewidth',...
          lnwidth,'linestyle',style1,'color',color1);
     hold on
     p2 = plot(varargin{1}, varargin{3}, 'Marker',marker2,'linewidth',...
          lnwidth,'linestyle',style2,'color',color2);
     xlabel(xlab,  'FontSize', ftsize)
     ylabel(ylab1,  'FontSize',  ftsize)
     set(gca(),'fontsize',ftsize)
     hold off
end

%set(gcf,'units','normalized','outerposition',[0 0 1 1])
end