function [ tform_lwm21, FRE, TRE ] = imageRegistrationTrackMate( calfile, calcTRE,...
    rows,cols, roiData)
%--------------------------------------------------------------------------
%
% imageRegistrationTrackMate:
%   Calculates the lwm geometric transformation function. 
% 
% Description:
%   This function is used to calculate the mapping function from a stack of 
%   Tiff images [1]. Each image shows a pair of control points. Together,  
%   the control points form a grid of 26 rows and 24 columns in each  
%   emission channel (acceptor emission: left control point and donor  
%   emission: right control point).The position of each control point was
%   measured with the ImageJ plugin TrackMate [2] and saved under the name 
%   'Spots in tracks statistics.txt'.
% 
% Syntax:  
%   [ tform_lwm21, FRE, TRE ] = imageRegistrationTrackMate(calfile, ...
%                               calcTRE, rows, cols, roiData)
% 
% Inputs:
%   1. calfile - ‘Spots in tracks statistics.txt’ 
%   2. calcTRE - Binary variable that enables (true) or disables (false)
%      the calculation of the target registration error. 
%   3. row - Number of lines of the calibration grid
%   4. col - Number of columns of the calibration grid
%   5. roiData - Specifies the size of the calibration grid used 
%      Vector with entries:  x origin of image coordinates (0),
%      y origin of image coordinates (0), image width, image height
%
%   NOTE: gridSize  : unit um
%         pixCal    : unit um/pixel
% 
% Outputs:
%   1. tform_lwm21 - Spatial transformation derived from control point 
%      pairs, which is returned as TFORM structure. 
%   2. FRE - fiducial registration error.
%   3. TRE - target registration error.
% 
% See also: 
%   scriptGetFretTraces.m, trackMate2traces.m ,txtReadTrackMate.m 
%
% References:
%   [1] Churchman, L. S., et al. (2005). "Single molecule high-resolution 
%   colocalization of Cy3 and Cy5 attached to macromolecules measures 
%   intramolecular distances through time." Proc Natl Acad Sci U S A 102(5)
%   : 1419-1423.
%   [2] Tinevez, J. Y., et al. (2017). "TrackMate: An open and extensible 
%   platform for single-particle tracking." Methods 115: 80-90.
%
% Authors: 
%   -PG 06/12
%   -PG 12/13
%   -PG 06/14
%   -PG 08/14
%   -PG 10/14
%   -PG 01/16 fitgeotrans instaed of cp2tform
%   -PG 02/17 corrected bug in tform calc
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% START 
disp('Calculating mapping function...');

% calfile = 'D:\#smCellFRET_v1.1\examples\171025 Nanostage Calibration\Grid #01\Spots in tracks statistics.txt';
% calcTRE = 1;
% rows = 26;
% cols = 24;
% roiData=[256 256 256];
%%
%---
%--- Initialize parameters
%---
constants     = smCellConstants('');
dt            = constants.dt;
translCols    = constants.xTranslation; % correction of a xCoordinate shift 
translRows    = constants.yTranslation; % correction of a yCoordinate shift 
xOriginUTrack = 0;
yOriginUTrack = 0;
imWidth       = roiData(3); 

%---
%--- Get controlpoints for both channels
%---

% Load SPT coordinates
nFrames=rows*cols;
xInputPts = NaN(nFrames,1);
yInputPts = NaN(nFrames,1);
xBasePts  = NaN(nFrames,1);
yBasePts  = NaN(nFrames,1);
[~,data]  = trackMate2traces(calfile,dt,nFrames,[]);
idxRow=zeros(1,2);
nRows=0;
for j=1:nFrames
    idxCol  = data.x(:,j)>0;
    x  = data.x(idxCol,j);
    y  = data.y(idxCol,j);
    %if sum(idxCol)== 2 % only use entries with 2 detection points
        x1 = x(1);y1 = y(1);
        x2 = x(2);y2 = y(2);
        
        % Calculate nRows 
        tmpIdx  = find(data.x(:,j));
        idxRow(1,1) = max(tmpIdx);
        if idxRow(1,1) > idxRow(1,2)
            nRows=nRows+1;
        end
        
        %save previous row
        idxRow(1,2) = idxRow(1,1);
        
        % save grid coordinates 
        if min(x1, x2) == x1
            % Red Channel (Ch1)
            idx=j;
            xInputPts(idx)= x1;
            yInputPts(idx)= y1;
            xBasePts(idx) = x2;
            yBasePts(idx) = y2;
        else
            % Green Channel (Ch2)
            idx=j;
            xInputPts(idx)= x2;
            yInputPts(idx)= y2;
            xBasePts(idx) = x1;
            yBasePts(idx) = y1;
        end
    %end
end
% Combine data in a new matrix 'inputPts' and basePts
inputPts = zeros(nFrames,2);
basePts  = zeros(nFrames,2);
                                                       
for k=1:nFrames                                        % Linear Coordinate Transformation (trackMate -> Utrack)
    inputPts(k,1) = xInputPts(k)+xOriginUTrack;        % add 1 if grid is tracked with trackMate 
    inputPts(k,2) = yInputPts(k)+yOriginUTrack;        % (image origin is: 0,0) and the  
    basePts(k,1)  = xBasePts(k)-imWidth+xOriginUTrack; % receptors are tracked with 
    basePts(k,2)  = yBasePts(k)+yOriginUTrack;         % uTrack (image origin is 1,1 -> Matlab)
end

%%
% Translate grid in Ch2 (basePts,green) relative to the grid in Ch1 
% (inputPts, red) if the grid in Ch1 and Ch2 are shifted with resprect to 
% each other because of bad initial grid alignmet.
if translCols > 0 || translRows > 0
    basePts(:,1)=basePts(:,1)-translCols; %xCoordinate
    basePts(:,2)=basePts(:,2)-translRows; %yCoordinate
end

%%
%---
%--- Calculate tform
%---
cpChannel1=inputPts; % red
cpChannel2=basePts;  % green

%cpChannel1=basePts; moving points
%cpChannel2=inputPts; fixed points

% Use the lwm transformation (local weighted mean), when the distortion
% varies locally and piecewise linear is not sufficient.
% movingPoints and fixedPoints, which define matched control points in the
% moving and fixed images, respectively. 
% moving points: Ch2 (Donor)
% fixed Points:  Ch1 (Acceptor)
% The n closest points are used to infer a second degree polynomial 
% transformation for each control point pair.
tform_lwm21=fitgeotrans(cpChannel2,cpChannel1,'lwm',12);

%---Calculate fiducial registration error
trans_cpChannel1=transformPointsInverse(tform_lwm21,cpChannel1);

% Calc the fiducial registration error
FRE21=sqrt(sum(sum((cpChannel2-trans_cpChannel1).^2,2))/length(cpChannel2));
FRE  =FRE21*160; %160nm per pixel, unit nm

%---Calculate target registration error
if calcTRE == 1
    %Find number of control point
    number_cp=length(cpChannel1);
    %loop through the control points
    for i=1:number_cp
        %take out that control point
        remain_cp = [1:i-1 i+1:number_cp];
        
        %calc the transformation without the ith control point
        tform_lwm21=fitgeotrans(cpChannel2(remain_cp,:),cpChannel1(remain_cp,:),'lwm',12);
        
        %transform the leftout control point with the above found transformation
        trans_cpChannel1(i,:)=transformPointsInverse(tform_lwm21,cpChannel1(i,:));
    end
    
    TRE21=sqrt(sum(sum((cpChannel2-trans_cpChannel1).^2,2))/length(cpChannel2));
    TRE  =TRE21*160; %160nm per pixel, unit nm
    
else
    TRE=[];
end
%%




