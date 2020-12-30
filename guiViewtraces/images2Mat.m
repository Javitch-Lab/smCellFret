function [ imMat ] = images2Mat( fname, roiData, intRange)

%--------------------------------------------------------------------------
%
% images2Mat:
%   Subfunction of cellFretViewtraces that converts an image series into a  
%   3D pixel matrix.
%
% Description:
%   The function processes a series of images (fov1_0001.tif, 
%   fov1_0002.tif, fov1_0003.tif, ...) which are saved in the 'ImageData' 
%   folder.
% 
% Syntax:  
%   [ imMat ] = images2Mat( fname, roiData, intRange)
% 
% Inputs:
%   1. fname - Path name of the folder 'ImageData'
%   2. roiData - (optional) Vector with entries:  x origin of image 
%   coordinates, y origin of image coordinates, image width, image height
%   3. intRange - (optional) Specifies the contrast limits of the function
%   imadjust
% 
% Outputs:
%   1. imMat(x,y,frame) - the image series is saved in a 3D array with the 
%   three subscripts x (x-coordinate), y (y-coordinate), frame (frame number)  
% 
% Authors: 
%   - P.G. Aug. 2012
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

% Get number of Frames
nFrames=getNFrames([fname '\ImageData\']);
if nargin<2; roiData=[];end
if nargin<3; intRange=[];end

% Get image size
info=imfinfo([fname 'ImageData\fov1_' num2str(1,'%05i') '.tif']);
originalWidth=info.Width;
oriognialHeight=info.Height;

% Allocate memory for images
if ~isempty(roiData)
    xZero    = roiData(1);
    yZero    = roiData(2);
    imWidth  = roiData(3);
    imHeight = roiData(4);
    imMat=zeros(imHeight+1, imWidth+1, nFrames);
else
    yZero=0;
    xZero=0;
    imWidth=originalWidth;
    imHeight=oriognialHeight;
    imMat=zeros(oriognialHeight,originalWidth, nFrames);
end

% Save images in matrix
disp('Loading Images...');
parfor i=1:nFrames
    imFile=[fname 'ImageData\fov1_' num2str(i,'%05i') '.tif'];
    tmpIm = imread(imFile);
    % Calculate ROI submatrix
    if ~isempty(roiData)
        idx = zeros(oriognialHeight,originalWidth);
        idx(yZero:yZero+imHeight, xZero:xZero+imWidth) = 1;
        tmpIm = reshape(tmpIm(find(idx)), max(sum(idx)), max(sum(idx,2)));
    end
    % Adjust image intensity
    if ~isempty(intRange)
        tmpIm = imadjust(tmpIm,[intRange(1) intRange(2)],[0 1]);
    end
    imMat(:,:,i) = tmpIm;
end
