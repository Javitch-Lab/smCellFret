%--------------------------------------------------------------------------
%
% scriptPostSyncTraces.m:
%   Synchronizes all FRET traces in a fretTraces structure to the start 
%   point of the movie. 
% 
% Description:
%   In live cell single-molecule movies, particle tracks and FRET traces
%   have different start and end times. To generate FRET distribution 
%   histograms, the individual fret traces are synchronized with 
%   respect to their starting point. 
% 
% Syntax:  
%   scriptPostSyncTraces
% 
% Inputs:
%   Running scriptPostSyncTraces.m will prompt the user to select a single
%   'fretTraces*.mat' input file (e.g 'fretTracesDif_SegFre.mat')
% 
% Outputs:
%   By default, the results are saved as a traces object under the same 
%   file name as the input file with the extension _sync.traces (e.g. 
%   fretTraces_SegFre_sync.traces). Saving the data as a traces object has 
%   the advantage that it is compatible with the data format used in the 
%   single molecule FRET software SPARTAN [1].
% 
% See also: 
%   alphaCorrectCell.m, deltaCorrectCell.m, gammaCorrectCell.m
%
% References:
%   [1] Juette, M. F., et al. (2016). "Single-molecule imaging of 
%   non-equilibrium molecular ensembles on the millisecond timescale." 
%   Nature Methods 13(4): 341-344.
%
% Authors: 
%   - P.G. 2018
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Load files

% Get DIRECTORIES
clear;
[fileName,pathName] = uigetfile('*.mat','Select a fretTraces*.mat file');

fretTracesFile = [pathName fileName];

if ~fileName==0 %% && exist(fretTracesFile,'file')
    load(fretTracesFile);
else
    warning('file not available');
    return
end

% Set flag for post-synchronizing 
isSynch=true;

% Rename variable
tmp=fretTraces;
clearvars fretTraces

%% Determine time resolution
dt_s  = tmp.Ch1.time(2)-tmp.Ch2.time(1);
dt_ms = dt_s*1000;

%% Post-Synchronize Traces
h = waitbar(0,'Postsynchronizing Traces....');

[nTraces,nFrames]=size(tmp.Ch1.x);
data=struct('time',(1:nFrames)'*dt_ms,...
            'donor',zeros(1,nFrames),...
            'acceptor',zeros(1,nFrames),...
            'total',zeros(1,nFrames),...
            'fret',zeros(1,nFrames),...
            'nTraces',nTraces,...
            'nFrames',nFrames,...
            'traceMetadata',struct('ids',[]));
        
% raw traces
for i=1:nTraces
    waitbar(i / nTraces)
    startOfTrace =  tmp.Ch2.traceMetadata(i).startOfTrace;
    endOfTrace   =  tmp.Ch2.traceMetadata(i).endOfTrace;
    lenBaseline  =  tmp.Ch2.traceMetadata(i).lenBaseline;
    
    data.channelNames={'donor' 'acceptor' 'fret'};
    
    % Post Synchronize to startOfTrace
    if isSynch==1
        pstSynTime = startOfTrace-1;
        lenBacksetAcc = 0;
        lenBacksetDon = 0;
    else
        pstSynTime = 0;
        lenBacksetAcc = tmp.Ch1.traceMetadata(i).lenBackset;
        lenBacksetDon = tmp.Ch2.traceMetadata(i).lenBackset;
    end
    
    data.donor(i,startOfTrace-pstSynTime:endOfTrace+lenBaseline-pstSynTime) = ...
                    tmp.Ch2.int(i,startOfTrace:endOfTrace+lenBaseline);
                
    data.acceptor(i,startOfTrace-pstSynTime:endOfTrace+lenBaseline-pstSynTime) = ...
                    tmp.Ch1.int(i,startOfTrace:endOfTrace+lenBaseline);
    
    data.traceMetadata(i).ids             = tmp.Ch1.traceMetadata(i).ids;
    data.traceMetadata(i).lenBacksetAcc   = lenBacksetAcc;
    data.traceMetadata(i).startOfTraceAcc = tmp.Ch1.traceMetadata(i).startOfTrace-pstSynTime;
    data.traceMetadata(i).endOfTraceAcc   = tmp.Ch1.traceMetadata(i).endOfTrace-pstSynTime;
    data.traceMetadata(i).lenBaselineAcc  = tmp.Ch1.traceMetadata(i).lenBaseline;
    data.traceMetadata(i).traceLenAcc     = tmp.Ch1.traceMetadata(i).traceLen;
    
    data.traceMetadata(i).lenBacksetDon   = lenBacksetDon;
    data.traceMetadata(i).startOfTraceDon = tmp.Ch2.traceMetadata(i).startOfTrace-pstSynTime;
    data.traceMetadata(i).endOfTraceDon   = tmp.Ch2.traceMetadata(i).endOfTrace-pstSynTime;
    data.traceMetadata(i).lenBaselineDon  = tmp.Ch2.traceMetadata(i).lenBaseline;
    data.traceMetadata(i).traceLenDon     = tmp.Ch2.traceMetadata(i).traceLen;
    
    if isfield(tmp.Ch1.traceMetadata,'meanXIntAcc')
        data.traceMetadata(i).meanXIntAcc     = tmp.Ch1.traceMetadata(i).meanXIntAcc;
        data.traceMetadata(i).meanXIntDon     = tmp.Ch1.traceMetadata(i).meanXIntDon;
        data.traceMetadata(i).meanXSnrAcc     = tmp.Ch1.traceMetadata(i).meanXSnrAcc;
        data.traceMetadata(i).meanXSnrDon     = tmp.Ch1.traceMetadata(i).meanXSnrDon;
        data.traceMetadata(i).meanXCrt        = tmp.Ch1.traceMetadata(i).meanXCrt;
    end
    
end
data.total = data.donor + data.acceptor;
data.fret = data.acceptor./data.total; 
close(h);


%% Replace NaNs in FRET for further processing with Spartan
data.fret(isnan(data.fret))=0;
 
%% Save traces file as traces object for Spartan 
[~,name] = fileparts(fileName);
newFileName  = inputdlg({'Enter File Name:'},...               %prompt
                         'Save File As: ',...                  %dlg titel
                          [1 50],...                           %width of line
                          {[name '_sync.traces']}); %default answer

if isempty(newFileName)
    return
else
    outfile      = [pathName newFileName{1}];
    disp(['save data in: ' outfile]);
    saveTraces(outfile, data);
end

%% End
disp('Finished...');





