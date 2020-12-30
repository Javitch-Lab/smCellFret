function [ nFrames, dataDirCh1, dataDirCh2 ] = stk2tiff( filename, ...
    nChannels,startFrame, endFrame)
%--------------------------------------------------------------------------
%
% stk2tiff:
%   Function that converts and spectrally separates a stk file movie into a 
%   donor and a acceptor image series.
% 
% Description:
%   This script is a batch processing script that first converts a 
%   metamorph stk file into a sequence of single tiff images. Second, the
%   script separates the two spectrally different regions of the donor and 
%   acceptor in each image into two independent numbered image series.  The
%   acceptor image series is saved in the #Ch1/ImageData folder while the 
%   donor image series is saved in the #Ch2/ImageData folder.
% 
% Syntax:  
%   [ nFrames, dataDirCh1, dataDirCh2 ] = stk2tiff( filename, ...
%    nChannels,startFrame, endFrame)
% 
% Inputs:
%   The function requires the manual selection of the current working  
%   folder with the Matlab file explorer.
%   1. filename - File name of the metamorph stk file, consisting of a base
%      name and a 2-digit number extension (e.g. 40msG3mW80Pd100G_01). 
%   2. nChannels - Either 1 (for one-color images, no image splitting) or 2 
%      (for two-color FRET images). 
%   3. startFrame - First image of the image sequence to be converted and 
%      split.      
%   4. endFrame - Last image of the image sequence to be converted and 
%      split.    
% 
% Outputs:
%   1. nFrames - Number of total frames in a stk movie file.
%   2. dataDirCh1 - Directory where the acceptor data is saved. 
%   3. dataDirCh2 - Directory where the donor data is saved. 
% 
% Other m-files required: 
%   Subfunctions: The script uses the function tiffread.m [1].
% 
% See also: 
%   scriptStk2tiff.m, scriptGetFrettraces.m
%
% References:
%   [1] Nedelec, F., et al. (2001). "Dynamic concentration of motors in
%   microtubule arrays." Phys Rev Lett 86(14): 3192-3195.
%
% Permission: 
%   tiffread is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, version 3 of the License.
%   Copyright (C) 1999-2014 Francois Nedelec
%
% Authors: 
%   - P.G.
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

% Create directories
currentPath  = what;
basepath     = char(currentPath.path);
filenameSegs = textscan(filename,'%s %02s %s','delimiter','_');
switch nChannels
    case 1
        dataDirCh1   = [basepath '\#' char(filenameSegs{2}) 'Ch1']; 
        mkdir(dataDirCh1);
        imageDirCh1  = [dataDirCh1 '\ImageData'];
        mkdir(imageDirCh1);
    case 2
        dataDirCh1   = [basepath '\#' char(filenameSegs{2}) 'Ch1']; 
        mkdir(dataDirCh1);
        dataDirCh2   = [basepath '\#' char(filenameSegs{2}) 'Ch2']; 
        mkdir(dataDirCh2);
        imageDirCh1  = [dataDirCh1 '\ImageData']; 
        mkdir(imageDirCh1);
        imageDirCh2  = [dataDirCh2 '\ImageData'];
        mkdir(imageDirCh2);
    otherwise
        error('number of channels >2 is not supported');
end

% Allocate digits for filename 
enumString = repmat('0',99999,5);
for i = 1 : 9
    enumString(i,:) = ['0000' num2str(i)];
end
for i = 10 : 99
    enumString(i,:) = ['000' num2str(i)];
end
for i = 100 : 999
    enumString(i,:) = ['00' num2str(i)];
end
for i = 1000 : 9999
    enumString(i,:) = ['0' num2str(i)];
end
for i = 10000 : 99999
    enumString(i,:) = num2str(i);
end

% Allocate digits for filename 
%enumString = repmat('0',9999,4);
%for i = 1 : 9
%    enumString(i,:) = ['000' num2str(i)];
%end
%for i = 10 : 99
%    enumString(i,:) = ['00' num2str(i)];
%end
%for i = 100 : 999
%    enumString(i,:) = ['0' num2str(i)];
%end
%for i = 1000 : 9999
%    enumString(i,:) = num2str(i);
%end


% Load stk file
stack=tiffreadCell(filename);
nFrames    = length(stack);
imageSizeX = stack(1,1).width;
imageSizeY = stack(1,1).height;


%Set start frame
%startFrame=1+flagFirstFrameDark;
frameTag=0;
for i = startFrame:endFrame
    frameTag=frameTag+1;
    filename=(['fov1_' enumString(frameTag,:) '.tif']);
    imTmp=stack(i).data;
    switch nChannels
       case 1
           imCh1=imTmp;
           imwrite(imCh1,[imageDirCh1 '\' filename],'tif');
       case 2
           imCh1 = imTmp(1:imageSizeY,1:imageSizeX/2);
           imwrite(imCh1,[imageDirCh1 '\' filename],'tif');
           imCh2 = imTmp(1:imageSizeY,1+(imageSizeX/2):imageSizeX);
           imwrite(imCh2,[imageDirCh2 '\' filename],'tif');
    end
end

