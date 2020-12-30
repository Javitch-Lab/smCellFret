%--------------------------------------------------------------------------
%
% getIds.m:
%   The script returns trace identifiers as numbers to a workspace
%   structure variable.
% 
% Syntax:  
%   getIds
% 
% Inputs:
%   The function opens the file explorer and prompts the user to select one 
%   or more fretTraces files that have either the .mat or .traces 
%   extension.
% 
% Outputs:
%   The script saves the output in the workspace structure variable 
%   'results'. The variable has the following fields 
%    .path    - Contains the path of all selected folders.
%    .file    - Contains the name of all selected files.
%    .ids     - Column vector containing only the identifier number of the
%               trace. 
%    .nTraces - Column vector containing the number of traces of each 
%               selected dataset.
%
% Authors: 
%   - P.G. Dec 2011
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Set initial Parameters
clear;

%% -------------------------------------------------------------------------
% Get files from user
% -------------------------------------------------------------------------

% Get FILES
tfiles = getFiles('*.*','Select (*.mat) or (*.traces) files');



% Check if files array is empty
if isempty(tfiles)== true; return;end

% Show selection
for z=1:numel(tfiles)
    [path,name,ext]=fileparts(tfiles{z});
    disp(path);
    disp([name ext]);
end

nCells  = numel(tfiles);
results = struct([]);

%% Load data
f = waitbar(0,'Load Data ...');

% Loop over selected experiment folders
for j=1:nCells
    
    clear tmp;
    waitbar(j / nCells,f)
    
    % save filename
    tracesFile = tfiles{j};
    [path,name,ext]  = fileparts(tracesFile);
    results(j).path = path;
    results(j).file = name;
    
    switch ext
        case '.traces'
            fretTraces=loadTraces(tracesFile);
            ids = {fretTraces.traceMetadata.ids}.';
            
        case '.mat'
            load(tracesFile);
            ids = {fretTraces.Ch1.traceMetadata.ids}.';
    end
    
    % read Ids as real numbers
    nTraces = size(ids,1);
    idNumber=zeros(1,nTraces);
    for k=1:nTraces
        tmpId=textscan(ids{k}, '%*s %d','delimiter','_');
        idNumber(k)=cell2mat(tmpId);
    end
    results(j).ids = idNumber';
    results(j).nTraces = nTraces;
    

end
close(f)

%% Stat

nMol=sum([results.nTraces].');
disp(['nTraces = ' num2str(nMol)]);

%% End
