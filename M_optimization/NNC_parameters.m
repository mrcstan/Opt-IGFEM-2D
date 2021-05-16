%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Marcus Tan on 9/23/2015
%%% Copyright 2015 University of Illinois at Urbana-Champaign. All rights reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the parameters for the normalized normal
% constraint method for generating Pareto frontier
% Currently the construction of hyperplane points only works for nObjs == 2
% INPUT:
%   anchorPts: A matrix of unnormalized anchor points with the i-th column 
%              being the objective values corresponding to the minimization
%               of the i-th objective alone
%              Example: say the objective values corresponding to the 
%              minimization of the p-norm alone are respectively T^{(T)}
%              and P^{(T)} and those corresponding to the minimzation of 
%              the pressure drop alone are respectively T^{(P)} and P^{(P)}.
%               Then the matrix should be [T^{(T)},T^{(P)};P^{(T)},P^{(P)}].
%   nHyperPlanePts: Number of points in first objective direction on hyperplane 
%                   including the end points. 
%                   If the Pareto front of three objectives are needed, 
%                   the number of points in the other direction would be
%                   calculated by such that the spacings of the points are equal 
%                   in each direction.
%   hyperPlanePts: For two objective functions, see the following:
%                  i) Before 1/28/2017, the first nHyperPlanePts1-2 columns of
%                     hyperPlanePts are the interior hyperplane points while
%                     the last two are NaN columns
%                  ii) From 1/28/2017, the first and last columns of hyperPlanePts
%                      are the scaled anchor points, while the the remaining columns
%                      are the interior hyperplane points
% NOTE: Cannot handle more than three objectives
function [normalization,...
          scaledAnchorPts,...
          UtopiaLineVec,...
          hyperPlanePts] = NNC_parameters(anchorPts,nHyperPlanePts1)
    nObjs = size(anchorPts,1);
    if (nObjs ~= size(anchorPts,2))
        error('anchorPts must be a square matrix')
    end

    % check that each diagonal entry is the minimum in its row
    for i = 1:nObjs
        [~,ind] = min(anchorPts(i,:));
        if ind ~= i
            warning('each diagonal entry point of the anchorPts matrix should be the minimum of its row')
            break
        end
    end
    % calculate normalization factors
    normalization = nan(nObjs,1);
    for i = 1:nObjs
        normalization(i) = max(anchorPts(i,:))-anchorPts(i,i);
    end
    if any(normalization == 0)
        disp(normalization)
        error('Some normalization entries vanish') 
    end
    % translate and normalize anchor points
    scaledAnchorPts = nan(nObjs,nObjs);
    for i = 1:nObjs
        scaledAnchorPts(i,:) ...
            = (anchorPts(i,:)-anchorPts(i,i))/normalization(i);
    end
    if isnan(scaledAnchorPts)
        disp(scaledAnchorPts)
        error('Some scaled anchor points are NaN')
    end
    % calculate Utopian line vector
    UtopiaLineVec = nan(nObjs,nObjs-1);
    for i = 1:(nObjs-1)
        UtopiaLineVec(:,i) = scaledAnchorPts(:,end) ...
                             - scaledAnchorPts(:,i);
    end
    
    % calculate evenly distributed points on Utopian hyperplane
    delw = 1/(nHyperPlanePts1-1);
    switch nObjs
        case 2
            hyperPlanePts = nan(nObjs,nHyperPlanePts1);
            for i = 1:nHyperPlanePts1
                im1 = i-1;
                hyperPlanePts(:,i) = (1-im1*delw)*scaledAnchorPts(:,1) ...
                                     + im1*delw*scaledAnchorPts(:,2);
            end
        case 3
            nHyperPlanePts = nHyperPlanePts1*(nHyperPlanePts1+1)/2-3; % subtract the 3 corners
            hyperPlanePts = nan(nObjs,nHyperPlanePts);
            count = 0;
            
            % for arranging the points in a serpentine manner
            % starting from the first anchor point
            sign = -1; 
            
            for i = 0:(nHyperPlanePts1-2)
                if i == 0
                    jstart = 1;
                    jend = nHyperPlanePts1 - 2;
                else
                    jstart = 0;
                    jend = nHyperPlanePts1 - 1 - i;
                end
                alpha3 = i*delw;
                delw2 = (1-alpha3)/(nHyperPlanePts1-1-i);
                
                for j = jstart:jend
                    count = count + 1;
                    alpha1 = j*delw2;
                    % for arranging the points in a serpentine manner
                    if sign == 1
                        hyperPlanePts(:,count) = alpha1*scaledAnchorPts(:,1) ...
                                                +(1-alpha1-alpha3)*scaledAnchorPts(:,2) ...
                                                +alpha3*scaledAnchorPts(:,3);
                    else
                        hyperPlanePts(:,count) = (1-alpha1-alpha3)*scaledAnchorPts(:,1) ...
                                                +alpha1*scaledAnchorPts(:,2) ...
                                                +alpha3*scaledAnchorPts(:,3);
                    end
                end
                sign = sign*-1;
            end
        otherwise
            error('number of objectives not considered')
            
    end
end