%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 12/21/2014
%%% Copyright 2014 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates the intersection point velocities wrt the coordinates of B
% INPUT:
%   edgeA: [x1,x2;y1,y2] assumed fixed
%   edgeB: [x3,x4;y3,y4] 
%   tolVert: tan(positive angle) or the maximum positive slope beyond which
%            an edge is considered vertical
% OUTPUT: 
%   vx: [dxI/dx3,dxI/dy3,dxI/dx4,dxI/dy4]
%   vy: [dyI/dx3,dyI/dy3,dyI/dx4,dyI/dy4],
%       (xI,yI) is the intersection point
% ASSUMPTIONS:
%   i) mA ~= mB; if not they wouldn't have any intersection
% REMARK: Currently, it works for straight segments only
function [vx,vy] = intersection_velocities(edgeA,edgeB,tolVert)
mA = (edgeA(2,2) - edgeA(2,1))/(edgeA(1,2) - edgeA(1,1));
cA = edgeA(2,1) - mA*edgeA(1,1);

mB = (edgeB(2,2) - edgeB(2,1))/(edgeB(1,2) - edgeB(1,1));
cB = edgeB(2,1) - mB*edgeB(1,1);

if (abs(mB) <= tolVert)
    denom = edgeB(1,2) - edgeB(1,1);
    denomsq = (edgeB(1,2) - edgeB(1,1))^2;
    dcB_dx3 = edgeB(1,2)*(edgeB(2,1)-edgeB(2,2))/denomsq; %x4(y3-y4)/(x4-x3)^2
    dmB_dx3 = (edgeB(2,2)-edgeB(2,1))/denomsq; % (y4-y3)/(x4-x3)^2
    dcB_dy3 = edgeB(1,2)/denom; % x4/(x4-x3);
    dmB_dy3 = -1/denom; % 1/(x3-x4);
    dcB_dx4 = edgeB(1,1)*(edgeB(2,2)-edgeB(2,1))/denomsq; % x3(y4-y3)/(x4-x3)^2
    dcB_dy4 = -edgeB(1,1)/denom; % x3/(x3-x4)
    % DCM = [dcB/dx3,dcB/dy3,dcB/dx4,dcB/dy4;
    %        dmB/dx3,dmB/dy3,dmB/dx4,dmB/dy4];
    DCM = [dcB_dx3,dcB_dy3,dcB_dx4,dcB_dy4;
           dmB_dx3,dmB_dy3,-dmB_dx3,-dmB_dy3];
       
    if (abs(mA) <= tolVert)
        vx = [1/(mA-mB),(cB-cA)/(mA-mB)^2]*DCM;
        vy = mA*vx;
    else
        vx = zeros(1,4);
        vy = [1,edgeA(1,1)]*DCM;
    end   
else
    if (abs(mA) <= tolVert)
        dxIdx3 = (edgeB(1,2)*(edgeA(2,2)-edgeA(2,1))...
                 +edgeA(1,2)*(edgeA(2,1)-edgeB(2,2))...
                 +edgeA(1,1)*(edgeB(2,2)-edgeA(2,2)))...
                 /((edgeA(1,2)-edgeA(1,1))*(edgeB(2,1)-edgeB(2,2)));
        dxIdx4 = (edgeB(1,2)*(edgeA(2,2)-edgeA(2,1))...
                 +edgeA(1,2)*(edgeA(2,1)-edgeB(2,1))...
                 +edgeA(1,1)*(edgeB(2,1)-edgeA(2,2)))...
                 /((edgeA(1,1)-edgeA(1,2))*(edgeB(2,1)-edgeB(2,2)));
        vx = [dxIdx3,0,dxIdx4,0];
        vy = [mA*dxIdx3,0,mA*dxIdx4,0];
    else
        % case of two line segments being vertical
        %warning('two line segments are both almost vertical: slope A = %g, slope B = %g \n',mA,mB)
        vx = zeros(1,4);
        vy = zeros(1,4);
    end
end
end
