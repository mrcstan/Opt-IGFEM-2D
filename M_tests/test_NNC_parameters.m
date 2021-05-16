clear
close all
path(path,'../M_optimization')
%anchorPts = [25.49,246.5,137.2;
%            46994,261,1958;
%            0.0274,0.018,7.17e-3];
%anchorPts = [0,5,0;
%             0,0,1;
%             1,0,0];
anchorPts = [25.49,246.5;
            46994,261];
[NNC.normalization,...
 NNC.scaledAnchorPts,...
 NNC.UtopiaLineVec,...
 NNC.planePts]= NNC_parameters(anchorPts,6);

tol = 1e-8;
NNC.planePts = pareto_filter(NNC.planePts,tol);
nObjs = size(anchorPts,1);
switch nObjs
    case 2
        plot(NNC.scaledAnchorPts(1,:),NNC.scaledAnchorPts(2,:),'bs',...
             'markerfacecolor','b','markersize',12)
        hold on
        for i = 1:size(NNC.planePts,2)
            plot(NNC.planePts(1,i),NNC.planePts(2,i),'ro',...
                'markerfacecolor','r','markersize',12)
            text(NNC.planePts(1,i),NNC.planePts(2,i),num2str(i),...
                'fontsize',30)
        end
        set(gca,'fontsize',40)
        xlabel('T','fontsize',40)
        ylabel('P','fontsize',40)
    case 3
        plot3(NNC.scaledAnchorPts(1,:),NNC.scaledAnchorPts(2,:),...
              NNC.scaledAnchorPts(3,:),'bs',...
             'markerfacecolor','b','markersize',12)
        
        hold on
        for i = 1:size(NNC.scaledAnchorPts,2)
            text(NNC.scaledAnchorPts(1,i),NNC.scaledAnchorPts(2,i),...
                 NNC.scaledAnchorPts(3,i),num2str(i),'fontsize',30)
        end
        for i = 1:size(NNC.planePts,2)
            plot3(NNC.planePts(1,i),NNC.planePts(2,i),...
                  NNC.planePts(3,i),'ro',...
                'markerfacecolor','r','markersize',12)
            text(NNC.planePts(1,i),NNC.planePts(2,i),...
                 NNC.planePts(3,i),num2str(i),...
                'fontsize',30)
        end
        set(gca,'fontsize',40)
        xlabel('T','fontsize',40)
        ylabel('P','fontsize',40)
        zlabel('A','fontsize',40)
        view(45,45)
    otherwise
        error('number of objectives not considered')
end