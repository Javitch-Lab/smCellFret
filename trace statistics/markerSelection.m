function [ markerType,markerSize ] = markerSelection( value )
%--------------------------------------------------------------------------
%
% markerSelection.m:
%   A subfunction of lscriptFretTraceStat.mlx that assigns a marker symbol
%   to a data set
% 
% Description:
%   The function expects as input an integer number and returns as response
%   a marker symbol from a list of 13 different symbols and the marker
%   size.
% 
% Syntax:  
%   [ markerType, markerSize ] = markerSelection( value )
% 
% Inputs:
%   value         - Integer number specifying e.g. a data set.
% 
% Outputs:
%   1. markerType - String variable holding the assigned marker symbol.  
%   2. markerSize - Integer number that specifies the size of the marker
%                   symbol. The marker size increases if the input value 
%                   is > 13 and then >= 26, 39,...,n*13  
% 
% See also: 
%   lscriptFretTracesStat.mlx,  colorSelection.m
%
% Authors: 
%   - P.G.
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

markers = {'o',...  % Circle
           '+',...  % Plus Sign 
           '*',...  % Asterix
           '.',...  % Point
           'x',...  % Cross
           's',...  % Square
           'd',...  % Diamond
           '^',...  % Upward Triangle
           'v',...  % Downward Triangle
           '>',...  % Right Triangle
           '<',...  % Left Triangle
           'p',...  % Pentagram
           'h'};    % hexagram

nMarkers = size(markers,2);
markerSize = 30;
offset     = 10;
if value > nMarkers 
    b = (1+(nMarkers.*floor(value./nMarkers))/nMarkers);
    markerSize = b*markerSize-offset;
    value = mod(value,nMarkers);
    if value == 0
        value=nMarkers;
    end
else
    b=1;
    markerSize = b*markerSize;
end

markerType = markers{value};

end

