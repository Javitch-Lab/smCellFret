%--------------------------------------------------------------------------
%
% getCrtValue:
%   Script that retrieves the crosstalk values saved in the files 
%   fretTracesDif*.mat or fretTracesDif*.traces. 
% 
% Description:
%   The script is used to retrieve measured crosstalk values from multiple 
%   files and to save all values in a new workspace variable. The combined 
%   data can be imported into a graphics software and displayed e.g. as a 
%   histogram.  
% 
% Syntax:  
%   getCrtValue
% 
% Inputs:
%   The script prompts the user to select fretTracesDif* files with either
%   the .mat or .traces extension
% 
% Outputs:
%   crtValues - workspace variable containing the combined crosstalk values
%
% Author: 
%   - P.G. Aug 2020
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
if any(cellfun('isempty',tfiles))== true; return;end

% Show selection
for z=1:numel(tfiles)
    [path,name,ext]=fileparts(tfiles{z});
    disp(path);
    disp([name ext]);
end


%% Get Crt

nCells = numel(tfiles);
data    =struct([]);

% Initialize Intensity Histogram
edges = 0:0.01:0.5;
nBins=numel(edges)-1;
histSnrDon=zeros(1,nCells+1);
histSnrAcc=zeros(1,nCells+1);

%% Analyze Data 

% Initialize figures and waitbar
%figure(1); % crt val
f = waitbar(0,'Load Data ...');
crtValues=[];
% Loop over selected experiment folders
for j=1:nCells
    
    clear tmp;
    waitbar(j / nCells,f)
    
    % save filename
    tracesFile = tfiles{j};
    [path,name]  = fileparts(tracesFile);
    data(j).path = path;
    data(j).file = name;
    
    switch ext
        case '.traces'
            fretTraces=loadTraces(tracesFile);
            crosstalk = [fretTraces.traceMetadata.crosstalk]';
        case '.mat'
            load(tracesFile);
            crosstalk = [fretTraces.Ch1.traceMetadata.crosstalk]';
    end
    
    
    % Get Crt Values
    crosstalk = crosstalk(crosstalk ~=0);
    data(j).crosstalk = crosstalk;
    
    % Combine 
    crtValues=[crtValues;crosstalk]; %#ok<AGROW>
    
    
end
close(f)

%% combine data
