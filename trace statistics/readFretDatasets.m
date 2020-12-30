function cellData = readFretDatasets(tracesFile)
%--------------------------------------------------------------------------
%
% readFretDatasets.m:
%   A subfunction of lscriptFretTracesStat.mlx that imports FretTraces data
%   into the Matlab workspace.
% 
% Description:
%   The function combines data from different experiments into a single 
%   structure variable cellData.
%
% Syntax:  
%   cellData = readFretDatasets(tracesFile)
% 
% Inputs:
%   1. tracesFile - String variable that specifies which fretTraces file is
%      to be loaded. Only one of the following file names can be assigned: 
%      'fretTraces.mat' or 'fretTracesPst.mat'.
%   2. The file explorer prompts the user to select multiple pairs of 
%      #Ch1/#Ch2 folders. Each Ch2 folder must contain the file name 
%      assigned to the string variable 'tracesFile'. 
% 
% Outputs:
%   A structure variable cellData that contains the following fields:
%     .path          - path of each selected file
%     .acc           - a list of the selected acceptor project folders 
%                      using the naming convention #%%Ch1. The wildcard 
%                      characters %% refer to a two-digit integer 
%                      indicating the number of the experiment. 
%     .don           - a list of the selected donor project folders 
%                      using the naming convention #%%Ch2. The wildcard 
%                      characters %% refer to a two-digit integer 
%                      indicating the number of the experiment.
%     .nTraces       - List of integers giving the total number of traces
%                      per experiment. 
%     .nPstFlt       - List of integers representing the total number of 
%                      traces per experiment after post-filtering the data
%                      with scriptFretFilter.m.
%     .traceMetadata - List of structure variables containing the trace
%                      metadata of each experiment. 
%     .traceStatSel  - List of structure variables containing the filter 
%                      results and the selected trace metadata of each 
%                      experiment.
%     .traceStatRej  - List of structure variables containing the filter 
%                      results and the rejected trace metadata of each 
%                      experiment.
%     .density       - List of 6 element vectors. Each vector contains the
%                      density information of the acceptor (1: number of 
%                      spots, 2: cell area, 3: spot density) and the donor
%                      (4: number of spots, 5: cell area, 6: spot density). 
% 
% See also: 
%   lscriptFretTracesStat.mlx
%
% Authors: 
%   - P.G.
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% For test purposes
%   clear;
%   tracesFile = 'fretTracesPst.mat';

%% Load Data

% Get directories from user
dirInfo = getDirs('Select pairs of #...CH1 and #...Ch2 folders');

% Initialize CellData
nFolders = numel(dirInfo(1,:));
cellData=struct('path',[],...
                'acc',[],...
                'don',[],...
                'nTraces',[],...
                'nPstFlt',[],...
                'traceMetadata',[],...
                'traceStatSel',[],...
                'traceStatRej',[],...
                'density',[]);

% Read dir information and update cellData structure            
folders=1; row=1;
while folders <= nFolders
    % Set path  
    cellData(row).path = dirInfo{2,folders};
    % Set channel names
    if contains(dirInfo{1,folders},'Ch1')
    cellData(row).acc = dirInfo{1,folders};  
    elseif contains(dirInfo{1,folders},'Ch2')
    cellData(row).don = dirInfo{1,folders};
    end

    folders=folders+1;

    if folders<=nFolders && ...
       ~strncmp(dirInfo{1,folders},dirInfo{1,folders-1},5)
        row=row+1;
    end

end

if isempty(dirInfo); cellData=[]; return; end

% Display directory names in command window
disp([cellData(1).path]);
disp('   ');
disp([cellData.acc]);
disp([cellData.don]);

%% Collecting All FretTraces.mat in one file called CombFret

h = waitbar(0,'Importing Data...');

for row=1:size(cellData,2)

    % Get Density Data: Searches for the density file in each dir
    chName={'acc','don'};
    for ch=1:2
        fn=chName{ch};
        if ~isempty(cellData(row).(fn))
            path = [cellData(row).path filesep cellData(row).(fn)];
            file = [path filesep 'density.txt'];

            if exist(file,'file')
                fileID  = fopen(file,'r');
                density = textscan(fileID,'%f','Delimiter','\t');
                fclose(fileID);
                density=cell2mat(density);
                cellData(row).density=density';
            end
           
        end
    end

    % Load single Molecule Data
    switch tracesFile
        case 'fretTraces.mat'

            % Load Data
            path = [cellData(row).path filesep cellData(row).don];
            load([path filesep tracesFile]);
            % Update cellData structure
            cellData(row).traceMetadata = ... 
                fretTraces.Ch1.traceMetadata(:,:);
            cellData(row).nTraces = size(fretTraces.Ch1.traceMetadata,2);

        case 'fretTracesPst.mat'

            % Load Data
            path = [cellData(row).path filesep cellData(row).don];
            load([path filesep tracesFile]);
            
            % Update cellData structure
            traceStatSel = fretTraces.Fret.traceStatSel;
            traceStatRej = fretTraces.Fret.traceStatRej;
            cellData(row).nTraces = ...
                size(traceStatSel,2) + size(traceStatRej,2);

            pstFltIdx = [fretTraces.Fret.traceStatSel.filterVal];
            cellData(row).nPstFlt = sum(pstFltIdx);

            cellData(row).traceMetadata = fretTraces.Ch1.traceMetadata(:,:);
            cellData(row).traceStatSel= fretTraces.Fret.traceStatSel;
            cellData(row).traceStatRej=fretTraces.Fret.traceStatRej;

        case 'smTracesCh1.mat'

            % Load Data
            path = [cellData(row).path filesep cellData(row).acc];
            load([path filesep tracesFile]);
            
            % Density
            cellData(row).density=0;
            
            % Update cellData structure
            cellData(row).traceMetadata = ... 
                smTraces.Ch1.traceMetadata(:,:);
            cellData(row).nTraces = size(smTraces.Ch1.traceMetadata,2);
            
        case 'smTracesCh1_DSel.mat'

            % Load Data
            path = [cellData(row).path filesep cellData(row).acc];
            load([path filesep tracesFile]);

            % Density
            cellData(row).density=0;

            % Update cellData structure
            cellData(row).traceMetadata = ... 
                smTraces.Ch1.traceMetadata(:,:);
            cellData(row).nTraces = size(smTraces.Ch1.traceMetadata,2); 
            
            
        case 'smTracesCh2.mat'
            
            % Load Data
            path = [cellData(row).path filesep cellData(row).don];
            load([path filesep tracesFile])
            
            % Density
            cellData(row).density=0;
            
            % Update cellData structure
            cellData(row).traceMetadata = ... 
                smTraces.Ch2.traceMetadata(:,:);
            cellData(row).nTraces = size(smTraces.Ch2.traceMetadata,2);    

        case 'smTracesCh1Pst.mat'
            
            % Load Data
            path = [cellData(row).path filesep cellData(row).acc];
            load([path filesep tracesFile])
            
            % Density
            cellData(row).density=0;
            
            % Update cellData structure
            if isfield(smTraces.filter,'traceStatSel')
                traceStatSel = smTraces.filter.traceStatSel;
                traceStatRej = smTraces.filter.traceStatRej;
                cellData(row).nTraces = ...
                    size(traceStatSel,2) + size(traceStatRej,2);
            else
                traceStatRej = smTraces.filter.traceStatRej;
                cellData(row).nTraces = size(traceStatRej,2);
            end
                
            if isfield(smTraces.filter,'traceStatSel')
                pstFltIdx = [smTraces.filter.traceStatSel.filterVal];
                cellData(row).nPstFlt = sum(pstFltIdx);
            else
                cellData(row).nPstFlt = 0;
            end
            
            if isfield(smTraces.Ch1,'traceMetadata')
                cellData(row).traceMetadata = smTraces.Ch1.traceMetadata(:,:);
            else
                cellData(row).traceMetadata = [];
            end
            
            if isfield(smTraces.filter,'traceStatSel')
                cellData(row).traceStatSel  = smTraces.filter.traceStatSel;
            else
                cellData(row).traceStatSel  = [];
            end
            
            cellData(row).traceStatRej  = smTraces.filter.traceStatRej;

        case 'smTracesCh2Pst.mat'
            % Load Data
            path = [cellData(row).path filesep cellData(row).don];
            load([path filesep tracesFile])
            
            % Density
            cellData(row).density=0;
            
            % Update cellData structure
            if isfield(smTraces.filter,'traceStatSel')
                traceStatSel = smTraces.filter.traceStatSel;
                traceStatRej = smTraces.filter.traceStatRej;
                cellData(row).nTraces = ...
                    size(traceStatSel,2) + size(traceStatRej,2);
            else
                traceStatRej = smTraces.filter.traceStatRej;
                cellData(row).nTraces = size(traceStatRej,2);
            end
            
            if isfield(smTraces.filter,'traceStatSel')
                pstFltIdx = [smTraces.filter.traceStatSel.filterVal];
                cellData(row).nPstFlt = sum(pstFltIdx);
            else
                cellData(row).nPstFlt = 0;
            end
            
            
            if isfield(smTraces.Ch2,'traceMetadata')
                cellData(row).traceMetadata = smTraces.Ch2.traceMetadata(:,:);
            else
                cellData(row).traceMetadata = [];
            end
            
            
            if isfield(smTraces.filter,'traceStatSel')
                cellData(row).traceStatSel  = smTraces.filter.traceStatSel;
            else
                cellData(row).traceStatSel  = [];
            end
            
            cellData(row).traceStatRej=smTraces.filter.traceStatRej;
            
        otherwise
            disp('no file available');

    end

    waitbar(row / nFolders)
end
close(h)
    %% End of Function
   
    empty_elems = arrayfun(@(s) isempty(s.nTraces),cellData);
    cellData(empty_elems)=[];

%%
% end

