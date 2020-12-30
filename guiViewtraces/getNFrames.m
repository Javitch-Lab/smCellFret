function nFrames = getNFrames(imageDir) 
%--------------------------------------------------------------------------
%
% getNFrames:
%   Subfunction of images2Mat.m which returns the number of frames in the 
%   folder 'ImageData'
% 
% Syntax:  
%   nFrames = getNFrames(imageDir) 
% 
% Inputs:
%   imageDir - pathname of the 'ImageData' folder
% 
% Outputs:
%   nFrames - variable with the assigned value of the number of images in
%   the 'ImageData' folder.
%
% See also: 
%   images2Mat.m
%
% Authors: 
%   - P.G. Oct 2011
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

% Get number of images
listing=ls([imageDir '*.tif']);
nFrames=size(listing,1);

end

