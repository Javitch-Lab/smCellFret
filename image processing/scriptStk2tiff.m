%--------------------------------------------------------------------------
%
% scriptStk2tiff:
%   Script that converts and spectrally separates a stk file movie into a 
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
%   scriptStk2tiff
% 
% Inputs:
%   The script will prompt the user to select one or multiple Metamorph stk
%   files.
% 
% Outputs:
%   1. acceptor image series saved in the #Ch1/ImageData folder
%   2. donor image series saved in the #Ch2/ImageData folder
%
% Other m-files required: 
%   Subfunctions: The script uses the function tiffread.m [1]. 
% 
% See also: 
%   stk2tiff.m, scriptGetFrettraces.m
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

%% Initialze
clear;
% Number of emission detection channels
% For single color emission : 1
% For dual color emission (acceptor left / donor right): 2 
nEmissions  = 2;

% frame of the stk file with which to begin
startFrame = 2; % first frame is dark.

%frame of the stk file with which to end
endFrame   = 4001;

%% Select multiple stk files
files = getFiles('*.stk','Highlight and Open *.stk files (HIT CANCEL TO END SELECTION)');

%% Allocate digits for filename 
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

%% Convert stk file to tiff sequences 
for i=1:length(files)
    
    [basepath,stk_file,ext] = fileparts(files{i});
    
    chId = textscan(stk_file,'%s %s','delimiter','_');
    
    switch nEmissions
        case 1
            dataDirCh1   = [basepath '\#' char(chId{2}) 'Ch1']; 
            mkdir(dataDirCh1);
            imageDirCh1  = [dataDirCh1 '\ImageData'];
            mkdir(imageDirCh1);
        case 2
            dataDirCh1   = [basepath '\#' char(chId{2}) 'Ch1']; 
            mkdir(dataDirCh1);
            dataDirCh2   = [basepath '\#' char(chId{2}) 'Ch2']; 
            mkdir(dataDirCh2);
            imageDirCh1  = [dataDirCh1 '\ImageData']; 
            mkdir(imageDirCh1);
            imageDirCh2  = [dataDirCh2 '\ImageData'];
            mkdir(imageDirCh2);
        otherwise
            error('number of channels >2 is not supported');
    end

    % Load Metamorph stk file and convert into a series of individual tif images
    stack   = tiffread([basepath filesep stk_file ext]);
    nFrames = length(stack);
    width   = stack(1,1).width;
    height  = stack(1,1).height;

    %Split images vertically into donor and acceptor images   
    frameTag=0;
    for j = startFrame:endFrame
        frameTag = frameTag+1;
        tif_file = (['fov1_' enumString(frameTag,:) '.tif']);
        imTmp    = stack(j).data;
        switch nEmissions
           case 1
               imCh1=imTmp;
               imwrite(imCh1,[imageDirCh1 '\' tif_file],'tif');
           case 2
               imCh1 = imTmp(1:height,1:width/2);
               imwrite(imCh1,[imageDirCh1 '\' tif_file],'tif');
               imCh2 = imTmp(1:height,1+(width/2):width);
               imwrite(imCh2,[imageDirCh2 '\' tif_file],'tif');
        end
    end
end

