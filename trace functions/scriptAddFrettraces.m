%--------------------------------------------------------------------------
%
% scriptAddFretTraces.m:
%   The script combines data files of type 'fretTraces*.mat' from several
%   experiment folders into a single data file to perform an analysis of 
%   all files together.
%
% Description:
%   The script is used when generating FRET histograms and analyzing track 
%   lifetimes.
% 
% Syntax:  
%   scriptAddFretTraces
% 
% Inputs:
%   The file explorer will propmt the user to select multiple
%   'fretTraces*.mat'.   
% 
% Outputs:
%   The script prompts the user to enter a file name or confirm the default 
%   name (the same name as the input file with the extension _add.mat'). 
%   The file with all combined fretTraces data is then saved under the 
%   specified file name in the last selected experiment folder.
%
% Authors: 
%   - P.G.  2018
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% Load files

% Get FILES
clear;
files = getFiles('*.mat','Select fretTraces.mat files to add');

% Check if files array is empty
if any(cellfun('isempty',files))== true; return;end

% Show selection in Cmd Window
for z=1:numel(files)
    [path,name,ext]=fileparts(files{z});
    disp(path);
    disp([name ext]);
end

% Create an output filename
[path,name,ext] = fileparts(files{end});
fName=[name ext];
fName=strrep(fName,'.','_add.');

%% Combine Arrays 
h = waitbar(0,'Concatenating Data Arrays...');

% Filed Names 
fretTracesFields  = {'x','y','xCorr','yCorr','int','snr',...
                        'idlTotal','total','pstFltVal'};

% Concatenate data arrays
tmp.Ch1=struct ('x',[],...
                'y',[],...
                'xCorr',[],...
                'yCorr',[],...
                'int',[],...
                'snr',[],...
                'traceMetadata',[]);
            
tmp.Ch2=struct ('x',[],...
                'y',[],...
                'int',[],...
                'snr',[],...
                'traceMetadata',[]);
            
tmp.Fret=struct ('idlTotal',[],...
                'total',[],...
                'traceStatSel',[],...
                'pstFltVal',[]);
            
metaFN = {'ids',...
          'lenBackset',...
          'startOfTrace',...
          'endOfTrace',...
          'lenBaseline',...
          'traceLen',...
          'meanInt',...
          'meanXIntAcc',...
          'meanXIntDon',...
          'meanXSnrAcc',...
          'meanXSnrDon',...
          'meanXCrt',...
          'crosstalk'}; 

chName = {'Ch1','Ch2'};
      
            
nExperiments = numel(files);
allTraces=0; 
for j=1:nExperiments
    waitbar(j / nExperiments)

    fretFile = files{j};
    load(fretFile);
   
    % nTraces
    nTraces = size(fretTraces.Ch1.x,1);
        
    % Load data
    for i=1:numel(fretTracesFields)
              
        % Donor 
        fn=fretTracesFields{i};
        if isfield(fretTraces.Ch1,fn)
            tmp.Ch1.(fn)=vertcat(tmp.Ch1.(fn),fretTraces.Ch1.(fn));
        end
        
        % Acceptor
        if isfield(fretTraces.Ch2,fn)
            tmp.Ch2.(fn)=vertcat(tmp.Ch2.(fn),fretTraces.Ch2.(fn));
        end
        
        % FRET 
        if isfield(fretTraces,'Fret') && isfield(fretTraces.Fret,fn)
            tmp.Fret.(fn)=vertcat(tmp.Fret.(fn),fretTraces.Fret.(fn));
        end
    end
    
    % Update Metadata
    % Traces Loop
    for i=1:nTraces
        allTraces=allTraces+1;
        %Channels Loop
        for c=1:2
            Ch=chName{c};
            % Fileds Loop
            for k=1:numel(metaFN)
                fn=metaFN{k};
                if isfield(fretTraces.(Ch).traceMetadata,fn)
                    tmp.(Ch).traceMetadata(allTraces).(fn) =...
                    fretTraces.(Ch).traceMetadata(i).(fn);
                else
%                     tmp.(Ch).traceMetadata(allTraces).(fn) = [];
                end
            end 
            % End Fields
        end
        % End Channles
    end
    % End Traces
end
% End Experiments

close(h)

% Set time
tmp.Ch1.time = fretTraces.Ch1.time;
tmp.Ch2.time = fretTraces.Ch2.time;

% Reorder Fields 
tmp.Ch1 = orderfields(tmp.Ch1,fretTraces.Ch1);
tmp.Ch2 = orderfields(tmp.Ch2,fretTraces.Ch2); 

% Check the dimensions
try
    assert(size(tmp.Ch1.x,1)==size(tmp.Ch2.x,1))
catch
    warning('Acceptor and Donor file have different sizes');
end

clear fretTraces;
fretTraces = tmp;
clear tmp;
%% Save Combined Data
prompt = {'Enter Filename'};
title = 'Save File As...';
definput = {fName};
answer = inputdlg(prompt,title,[1 40],definput);
outfile = [path filesep answer{:}];
disp(outfile);
save(outfile,'fretTraces','-v7.3');

%% End

