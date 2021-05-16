%close all
path(path,'../M_optimization')
tolVert = tan(89.9*pi/180);
tolItrsect = 1e-13;
%{
edgeA = [-1,1;
         -1,1];
edgeB = [-1,1;
          1,-1];
%}

% both edges are random
%{
edgeA = rand(2,2);
edgeB = rand(2,2);
intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
while (isempty(intPtsOld))
    edgeA = rand(2,2);
    edgeB = rand(2,2);
    intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
end
%}

% edge A is vertical
%{
edgeA = [0,0; 
        -1,1];
%{
edgeB = [-1,1;
         -1,1];
%}
edgeB = 2*rand(2,2)-1;
intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
while (isempty(intPtsOld))
    edgeB = 2*rand(2,2)-1;
    intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
end
%}

% edge B is vertical
%{
edgeB = [0,0; 
        -1,1];
%}
%{
edgeA = [-0.5,0.5;
         -1,1];
%}
%{
edgeA = 2*rand(2,2)-1;
intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
%
while (isempty(intPtsOld))
    edgeA = 2*rand(2,2)-1;
    intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
end
%
%}


% edge A is horizontal
%
edgeA = [0,0; 
        -1,1];
%{
edgeB = [0,0;
         -1,1];
%}
edgeB = 2*rand(2,2)-1;
intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
while (isempty(intPtsOld))
    edgeB = 2*rand(2,2)-1;
    intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
end
%

% edge B is horizontal
edgeB = [0,0; 
        -1,1];
edgeA = 2*rand(2,2)-1;
intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
while (isempty(intPtsOld))
    edgeA = 2*rand(2,2)-1;
    intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);
end

%edgeA = [0.1,0.11;
%         0.17,0.17];
%edgeA = [0.1,0.11;
%         0.17,0.18];
%edgeA = [0,0.11;
%         0.11,0.11];     
%edgeA = [0.1,0.11;
%         0.11,0.12];

%edgeB = [0.105,0.105;
%         0.172,0.112];
%edgeB = [0 0.18 
%        0.105 0.172]';  
%edgeB = [0.105 0.112 
%        0.047 0.085]';
%intPtsOld = intersect_edges(edgeA(:)',edgeB(:)',tolItrsect);     


[vx,vy] = intersection_velocities(edgeA,edgeB,tolVert);


figure
plot(edgeA(1,:),edgeA(2,:),'k-','linewidth',2)
hold on
plot(edgeB(1,:),edgeB(2,:),'r-','linewidth',2)

plot(intPtsOld(1),intPtsOld(2),'ko','linewidth',2,'markersize',8)

del = 1e-7;

FDvx = nan(1,4);
FDvy = nan(1,4);
for i = 1:4
    newEdgeB = edgeB;
    newEdgeB(i) = newEdgeB(i) + del;
    intPtsNew = intersect_edges(edgeA(:)',newEdgeB(:)',tolItrsect);
    FDvx(i) = (intPtsNew(1) - intPtsOld(1))/del;
    FDvy(i) = (intPtsNew(2) - intPtsOld(2))/del;
end

disp('vx = ')
disp(vx)
disp('FDvx = ')
disp(FDvx)

disp('vy = ')
disp(vy)
disp('FDvy = ')
disp(FDvy)

fprintf('vel x max diff = %g \n',max(abs(vx - FDvx)))
fprintf('vel y max diff = %g \n',max(abs(vy - FDvy)))

fprintf('vel x max rel diff = %g \n',max(abs(1 - FDvx./vx)))
fprintf('vel y max rel diff = %g \n',max(abs(1 - FDvy./vy)))