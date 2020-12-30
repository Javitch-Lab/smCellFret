function deleteTraces( directory )
%--------------------------------------------------------------------------
%
% deleteTraces.m:
%   Subfunction of scriptDeleteTraces.m
% 
% Description:
%   Deletes all traces in the directory (but not subdirectories)
% 
% Syntax:  
%   deleteTraces( directory )
% 
% Inputs:
%   directory - Path of the folder where all *.traces files should be
%   deleted.
% 
% Outputs:
%   Deletes all trace files in the selected directory
%   
% See also: 
%   scriptDeleteTraces.m
%
% Authors: 
%   - Jozsef Meszaros 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------
    tracesFiles = dir( directory )
    boolean_check = arrayfun( @(x) and(~x.isdir,~isempty(strfind(x.name,'.traces'))), tracesFiles )
    ch_folders = arrayfun( @(x) fullfile(tracesFiles(x).folder,tracesFiles(x).name), find(boolean_check==1) , 'UniformOutput', false);

    % Deletion Step %
    cellfun( @(x) delete(x), ch_folders );

end