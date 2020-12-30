function traces2mat(dirArray, donorTracking)
%--------------------------------------------------------------------------
%
% traces2mat.m:
%   subfunction of scriptGetFretTraces.m 
% 
% Description:
%   Combines the content of 'singleTracesDon_phc.traces' and 
%   singleTracesAcc_phc.traces' into a single fretTraces structure and  
%   saves the result as a fretTraces.mat file.   
% 
% Syntax:  
%   traces2mat(dirArray, donorTracking)
% 
% Inputs:
%   1. dirArray      - 2D cell array containing the selected project  
%                      folders #Ch1 and #Ch2 in the first row and the 
%                      corresponding path of the project folders in the 2nd 
%                      row (see example below).  
%   2. donorTracking - Logical variable that can either be set to TRUE 
%                      (donor tracking was used) or FALSE (no donor 
%                      tracking).
%  
% 
% Outputs:
%   The function combines the fretTraces structure of the donor ( saved in 
%   '...\#Ch2\singleTracesDon_phc.traces' ) and the acceptor (saved in '...
%   \#Ch1\singleTracesAcc_phc.traces' ) to a single fretTraces structure 
%   and saves the result as a 'fretTraces.mat' file in the #Ch2 folder of
%   the respective project. 
%
%   NOTE: For a description of the new fretTraces structure variable, see
%   the function header of scriptGetFretTraces.m
% 
% Example: 
%   Cell array 'dirArray' of #Ch1 and #Ch2 project folders and their path 
%   as generated by the function getDirs.m
%   dirArray = {'#01Ch1',      '#01Ch2',      '#02Ch1',      '#02Ch2';
%               'path_#01Ch1', 'path_#01Ch2', 'path_#02Ch1', 'path_#02Ch2'}
% 
% See also: 
%   scriptGetFretTraces.m, getDirs.m
%
% Authors: 
%   - P.G. 2017
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


disp('convert traces files...'); 
nExperiments     = 1:2:length(dirArray(1,:));
nChannelsTracked = 1+donorTracking;

% Loop through the experiments
for i=nExperiments
    fretTraces=[];
    basepath=[dirArray{2,i}];
%     try
        % Loop through the Acc / Don  Channel of the i-th experiment
        for j=1:2
            % Select channel specific directory
            if j==1
                dataSetCh1 = dirArray{1,i};
                disp([dataSetCh1 '...']);
                currentDir = [basepath '\' dataSetCh1 '\'];
                filename   = 'singleTracesAcc_phc.traces';
                %idxLog=[];

            elseif j==2
                dataSetCh2   = dirArray{1,i+1};
                disp([dataSetCh2 '...']);
                currentDir   = [basepath '\' dataSetCh2 '\'];
                if nChannelsTracked==1
                    filename = 'singleTracesDonA_phc.traces';
                elseif nChannelsTracked==2
                    filename = 'singleTracesDon_phc.traces';
                end
            end

            % Load Data and save to a new structure
            data = loadTracesCell([currentDir filename]);
            ch   = ['Ch' num2str(j)];

            structNames = fieldnames(data);
            nTraces = size(data.x,1);
            %if isempty(idxLog); idxLog=true(1,nTraces);end
            fretTraces.(ch).time  = data.time';
            fretTraces.(ch).x     = data.x;
            fretTraces.(ch).y     = data.y;

            if isfield(data,'xCorr')
                fretTraces.(ch).xCorr = data.xCorr;
                fretTraces.(ch).yCorr = data.yCorr;
            end

            fretTraces.(ch).int   = data.int;
            if contains([structNames{:}],'snr')
                fretTraces.(ch).snr  = data.snr;
            end

            for l=1:nTraces
                fretTraces.(ch).traceMetadata(l).ids          = ...
                    data.traceMetadata(l).ids;
                fretTraces.(ch).traceMetadata(l).lenBackset   = ...
                    data.traceMetadata(l).lenBackset;
                fretTraces.(ch).traceMetadata(l).startOfTrace = ...
                    data.traceMetadata(l).startOfTrace;
                fretTraces.(ch).traceMetadata(l).endOfTrace   = ...
                    data.traceMetadata(l).endOfTrace;
                fretTraces.(ch).traceMetadata(l).lenBaseline  = ...
                    data.traceMetadata(l).lenBaseline;
                fretTraces.(ch).traceMetadata(l).traceLen     = ...
                    data.traceMetadata(l).traceLen;

                if isfield(data.traceMetadata,'meanInt')
                fretTraces.(ch).traceMetadata(l).meanInt     = ...
                    data.traceMetadata(l).meanInt;
                end

                if isfield(data.traceMetadata,'meanXIntAcc')
                fretTraces.(ch).traceMetadata(l).meanXIntAcc = ...
                    data.traceMetadata(l).meanXIntAcc;
                end

                if isfield(data.traceMetadata,'meanXIntDon')
                fretTraces.(ch).traceMetadata(l).meanXIntDon = ...
                    data.traceMetadata(l).meanXIntDon;
                end

                if isfield(data.traceMetadata,'meanXSnrAcc')
                fretTraces.(ch).traceMetadata(l).meanXSnrAcc     = ...
                    data.traceMetadata(l).meanXSnrAcc;
                end

                if isfield(data.traceMetadata,'meanXSnrDon')
                fretTraces.(ch).traceMetadata(l).meanXSnrDon     = ...
                    data.traceMetadata(l).meanXSnrDon;
                end

                if isfield(data.traceMetadata,'meanXCrt')
                fretTraces.(ch).traceMetadata(l).meanXCrt     = ...
                    data.traceMetadata(l).meanXCrt;
                end

            end

            % Save Fret traces in the folder of nameCh2
            if j==2
                save([currentDir 'fretTraces.mat'],'-mat',...
                    'fretTraces','-v7.3');
            end
        end
%     catch
%         disp(currentDir)
%         warning('missing File')
%         continue
end 

%% End 


 