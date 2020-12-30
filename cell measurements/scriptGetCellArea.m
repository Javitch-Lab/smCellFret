%--------------------------------------------------------------------------
%
% scriptGetCellArea:
%   Calculates the cell area from the cell boundary.
% 
% Syntax:  
%   scriptGetCellArea.m
% 
% Inputs:
%   The File Explorer prompts the user to specify directories with file
%   names of type STD*.txt  
% 
% Outputs:
%   1. Figure with the cell boundary
%   2. workspace variable pArea:
%      pArea(i,1) = cell area in um^2
%      pArea(i,2) = cell area in pixel^2
% 
% Author: 
%   P.G. Oct 2018
%
% Copyright:
%   2011-2020 Javitch Lab - Columbia Univiversity / RFMH
%
%
% -------------------------------------------------------------------------


%% Initialize 
clear;
pixCal = 0.160; % um/pixel

%% Get directories from user
dirInfo  = getDirs;
nFolders = numel(dirInfo(1,:));
pArea = zeros(nFolders,2);

%% Calculate Polygon Area in um2
for i=1:nFolders
    % Load cell ROI
    path    = dirInfo{2,i};
    dataSet = dirInfo{1,i};
    fName = 'STD_40msG3mWR5080Pd100G_';
    
    tmpNum=textscan(dataSet,'%c %u %s');
    N=tmpNum{1,2};
    if N<10
        enumString = ['0' num2str(N)];
    else
        enumString = num2str(N);
    end
    
    fExt  = '.txt';
    
    file=[path filesep dataSet filesep fName enumString fExt ];
    xy=dlmread(file);
    
    % Draw Polygon
    figure
    plot(xy(:,1),xy(:,2));
    
    % Calc Polygon Area
    % in um^2 (1st Column)
    pArea(i,1)=polyarea(xy(:,1),xy(:,2))*pixCal^2;
    % in pixel^2 (2nd Column)
    pArea(i,2)=polyarea(xy(:,1),xy(:,2));
    
end

%% End
