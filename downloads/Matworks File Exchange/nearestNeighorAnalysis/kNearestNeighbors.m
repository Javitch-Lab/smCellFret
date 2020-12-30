function [neighborIds neighborDistances] = kNearestNeighbors(dataMatrix, queryMatrix, k)
%--------------------------------------------------------------------------
% Program to find the k - nearest neighbors (kNN) within a set of points. 
% Distance metric used: Euclidean distance
% 
% Usage:
% [neighbors distances] = kNearestNeighbors(dataMatrix, queryMatrix, k);
% dataMatrix  (N x D) - N vectors with dimensionality D (within which we 
%                         search for the nearest neighbors)
% queryMatrix (M x D) - M query vectors with dimensionality D
% k           (1 x 1) - Number of nearest neighbors desired
% 
% Example:
% a = [1 1; 2 2; 3 2; 4 4; 5 6]; a=rand(100,2) (dataMatrix)
% b = [1 1; 2 1; 6 2];           b=rand(100,2) (queryMatrix, Abfrage Matrix)
% [neighbors distances] = kNearestNeighbors(a,b,2);
% 
% Output (pos. in vector a):
% neighbors =
%     1     2
%     1     2
%     4     3
% Interpretation of results
%     k=1          k=2
%    a=1|b=1     a=2|b=1
%    a=1|b=2     a=2|b=2                                                          
%    a=4|b=3     a=3|b=3
% 
% distances =
%       k=1       k=2
%          0    1.4142
%     1.0000    1.0000
%     2.8284    3.0000
% Plot distribution of distances
% figure; 
% hist(distances)
%--------------------------------------------------------------------------

neighborIds = zeros(size(queryMatrix,1),k);
neighborDistances = neighborIds;

numDataVectors = size(dataMatrix,1);
numQueryVectors = size(queryMatrix,1);
% Plot data set1 
%x=dataMatrix(:,1);
%y=dataMatrix(:,2);
%qx=queryMatrix(:,1);
%qy=queryMatrix(:,2);
%figure
%plot(x,y,'.','Marker','s','MarkerSize',3,'color','r'); hold on
%Plot data set2 which will be used to find nearest neigbors to data set 1
%plot(qx,qy,'.','Marker','s','MarkerSize',3,'color','g');

for i=1:numQueryVectors,
    dist = sum((repmat(queryMatrix(i,:),numDataVectors,1)-dataMatrix).^2,2);
    [sortval sortpos] = sort(dist,'ascend');
    neighborIds(i,:) = sortpos(1:k);
    neighborDistances(i,:) = sqrt(sortval(1:k));
end