function [rgbColor,Name]=colorSelection(number)
%--------------------------------------------------------------------------
%
% colorSelection.m:
%   This function is a sub-function of lscriptFretTracesStat.mlx.
% 
% Description:
%   The function is used to assign a color from a palette of 12 colors to a
%   graphics object. 
% 
% Syntax:  
%   [rgbColor,Name] = colorSelection(number)
% 
% Inputs:
%   number      - Integer number between 1 and 12
% 
% Outputs:
%   1. rgbColor - rgb vector 
%   2. Name     - color name
% 
% See also: 
%   lscriptFretTracesStat.mlx, 
%
% Authors: 
%   - P.G.
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


% Define colormap
cMap                = [0,   0,   0      % black
                       1,   0,   0      % red
                       0,   1,   0      % green
                       0,   0,   1      % blue
                       0,   1,   0.75   % turquise 
                       1    0.75 0      % gold
                       0,   1,   1      % cyan
                       0,   0.5, 1      % light blue
                       0.5, 0,   1      % purple
                       1,   0,   1      % magenta
                       1,   0,   0.5    % pink
                       1,   0.3, 0]';   % orange

colorNames          = {'Black',...
                       'Red',...
                       'Green',...   
                       'Blue',...
                       'Turquise',...
                       'Gold',...
                       'Cyan',...
                       'L-Blue',...
                       'Purple',...
                       'Magenta',...
                       'Pink',...
                       'Orange'};
rgbColor = cMap(:,number);
Name  = colorNames{number}; 