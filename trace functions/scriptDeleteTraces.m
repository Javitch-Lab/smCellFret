%--------------------------------------------------------------------------
%
% scriptDeleteTraces.m:
%   The script finds all of the .traces files in an experiment directory
%   and deletes them from all subdirectories in order to conserve storage 
%   space (this will reduce data storage consumption by approximately 80%)
% 
% Syntax:  
%   scriptDeleteTraces.m
% 
% Inputs:
%   The script prompts the user to select a directory in which to remove 
%   .traces files. 
% 
% Other m-files required: 
%   Subfunctions: deleteTraces.m
% 
% Authors: 
%   - Jozsef Meszaros 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


% Uncomment if you wish to delete the traces files from the examples data
%
% dir_to_clean = 'C:\#smCellFRET_v1.1\examples\180109 mGlu2 Data'
if ~exist('dir_to_clean'); dir_to_clean = uigetdir(); end

myfolders = dir( dir_to_clean );

% Check myfolders is both a directory and it contains 'Ch'
boolean_check = arrayfun( @(x) and(x.isdir,~isempty(strfind(x.name,'Ch'))), myfolders );
% Extract only those directory addresses that contain Ch
ch_folders = arrayfun( @(x) fullfile(myfolders(x).folder,myfolders(x).name), find(boolean_check==1) , 'UniformOutput', false);
% Go into those directories and delete all traces files
deleteTraces( ch_folders{1} );
clearvars('dir_to_clean');