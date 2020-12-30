%--------------------------------------------------------------------------
%
% scriptFretTraces2Seg:
%   Generate fret and intensity time traces with only one diffusion mode
%   per trace. 
% 
% Description:
%   The script uses the motion classification results of the DC-MSS 
%   algorithm and parses FRET traces with multiple diffusion modes per 
%   trace into a segmented trace with only one diffusion mode per trace. 
% 
% Syntax:  
%  scriptFretTraces2Seg:
% 
% Inputs:
%   The script asks the user to select one or more #Ch2 folders containing 
%   file names of type 'fretTracesDif*.mat'.
%   
%   Note: Change the file name extension of 'fretTracesDif.mat' to e.g.
%   'fretTracesDif_bestFret.mat' in the code section 'Initial Settings'. 
% 
% Outputs:
%   The script saves the results under the name of the input file using the 
%   following extensions: *_SegImb.mat, *_SegCnf.mat, *_SegFre.mat and 
%   *_SegSpr.mat. 
%   1. fretTracesDif_SegImb.mat - trace that contains only immobile 
%      segments                                       
%   2. fretTracesDif_SegCnf.mat - trace that contains only segments of
%      confined diffusion                                   
%   3. fretTracesDif_SegFre.mat - trace that contains only segments of free
%      diffusion                                   
%   4. fretTracesDif_SegSpr.mat - trace that contains only segments of
%      directed motion    
%   
% Other m-files required: 
%   Subfunctions: uses local call back functions
%   
% Authors: 
%   -P.G April 2019 
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------
%% Initial Settings
clear
infile = 'fretTracesDif_allFret.mat';

%% Get Directories (only Ch2 directories required)

dirArray = getDirs(['Select #Ch2 Directories with saved files of type: ' infile]);
nExperiments = length(dirArray(1,:));
if nExperiments==0,return,end

%% Loop through all experiments

% Create waitbar
h = waitbar(0,'Get Frettraces Diffusion Segments ...');

% Set current Directopry
for j=1:nExperiments 
    
    waitbar(j / (nExperiments))
    
    % Set current Directopry
    basepath = dirArray{2,j};
    dataSet  = dirArray{1,j};
    currentDir = [basepath filesep dataSet filesep];
    
    disp(dataSet);
    
    if contains(dataSet,'Ch1')
        continue
    
    % Load post filtered fretTraces & SegResultsFinal (in Dir Ch2) 
    elseif contains(dataSet,'Ch2') && ...
            exist([currentDir infile],'file')
        
        % Load fretTraces into workspace
        load([currentDir infile]);
        
        %------------------------------------------------------------------
        %--- Immobile Segments 
        %------------------------------------------------------------------
        outfile    = [currentDir strrep(infile,'.mat', '_SegImb.mat')];
        getFrettraceSegments(fretTraces,outfile,'idlImb');
        
        %------------------------------------------------------------------
        %--- Confined Segments 
        %------------------------------------------------------------------
        outfile    = [currentDir strrep(infile,'.mat', '_SegCnf.mat')];
        getFrettraceSegments(fretTraces,outfile,'idlCnf');
        
        %------------------------------------------------------------------
        %--- Free Segments 
        %------------------------------------------------------------------
        outfile    = [currentDir strrep(infile,'.mat', '_SegFre.mat')];
        getFrettraceSegments(fretTraces,outfile,'idlFre');
        
        %------------------------------------------------------------------
        %--- Super Segments 
        %------------------------------------------------------------------
        outfile    = [currentDir strrep(infile,'.mat', '_SegSpr.mat')];
        getFrettraceSegments(fretTraces,outfile,'idlSpr');
        
    end
    
end
% Close Waitbar
close(h)

%% =======================================================================%
%                                                                         %
%                        Local Callback Functions                         %
%                                                                         %
%=========================================================================%


%-------------------------------------------------------------------------%
%                                                                         %
% Get Segments from fretTracesDif*.mat                                    %
%                                                                         %
%   input:  1. fretTraces structure variable                              %
%           2. output file name: e.g 'fretTracesDif_.mat'                 %
%           3. name of idealized segments:                                %
%               -> immobile: 'idlImb'  or ...                             %
%               -> confined: 'idlCnf'  or ...                             %
%               -> free:     'idlFre'  or ...                             %
%               -> super:    'idlSpr'                                     %
%                                                                         %  
%   output: fretTraces file of the motion type specified in 3.            %
%           -> fretTracesDif_SegImb.mat or ...                            %
%           -> fretTracesDif_SegCnf.mat or ...                            %
%           -> fretTracesDif_SegFre.mat or ...                            %
%           -> fretTracesDif_SegSpr.mat                                   %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%


function getFrettraceSegments(fretTraces,outfile,seg) 

% Get number of Traces
[nTraces,nFrames] = size(fretTraces.Ch1.x);

% Field Names 
fretTracesFieldsCh1  = {'x','y','xCorr','yCorr','int','snr'};

% Data arrays
fretSeg.Ch1=struct ('x',[],...
                'y',[],...
                'xCorr',[],...
                'yCorr',[],...
                'int',[],...
                'snr',[],...
                'traceMetadata',[]);
            
fretSeg.Ch2=struct ('x',[],...
                'y',[],...
                'int',[],...
                'snr',[],...
                'traceMetadata',[]);
            
        
%% Load data

% load idealization matrix 
idl = fretTraces.Diff.Ch1.(seg);

for i=1:numel(fretTracesFieldsCh1)
    
    % Acceptor 
    fn=fretTracesFieldsCh1{i};
    if isfield(fretTraces.Ch1,fn)
        fretSeg.Ch1.(fn)=vertcat(fretSeg.Ch1.(fn),fretTraces.Ch1.(fn).*idl);
    end
    
    % Donor
    if isfield(fretTraces.Ch2,fn)
        fretSeg.Ch2.(fn)=vertcat(fretSeg.Ch2.(fn),fretTraces.Ch2.(fn).*idl);
    end
end

%% Metadata
chName={'Ch1','Ch2'};

for i=1:nTraces
    startOfTrace = find(fretTraces.Diff.Ch1.(seg)(i,:)==1,1,'first');
    endOfTrace   = find(fretTraces.Diff.Ch1.(seg)(i,:)==1,1,'last');
    
    for  c=1:2
        Ch=chName{c};
        % ids
        fretSeg.(Ch).traceMetadata(i).ids          = ...
            fretTraces.(Ch).traceMetadata(i).ids;
        
        % lenBackset
        if fretTraces.(Ch).traceMetadata(i).lenBackset <=5
            fretSeg.(Ch).traceMetadata(i).lenBackset   = ...
                fretTraces.(Ch).traceMetadata(i).lenBackset;
        else
            fretSeg.(Ch).traceMetadata(i).lenBackset   = 5;
        end
            
        % startOfTrace
        fretSeg.(Ch).traceMetadata(i).startOfTrace = startOfTrace;
        
        % endOfTrace
        fretSeg.(Ch).traceMetadata(i).endOfTrace   = endOfTrace;
        
        % lenBaseline
        if endOfTrace + 5 <=nFrames
            fretSeg.(Ch).traceMetadata(i).lenBaseline = 5;
        else
            fretSeg.(Ch).traceMetadata(i).lenBaseline = 0;
        end
        
        %traceLen
        fretSeg.(Ch).traceMetadata(i).traceLen = ...
            endOfTrace - startOfTrace + 1;
        
    end
end

%% Remove empty rows in structure
filter = arrayfun(@(x) isempty(x.startOfTrace),fretSeg.Ch1.traceMetadata);

for c=1:2
    Ch=chName{c};
    
    % Keep field that are not empty
    fretSeg.(Ch).x            = fretSeg.(Ch).x(~filter,:);
    fretSeg.(Ch).y            = fretSeg.(Ch).y(~filter,:);
    if isfield(fretSeg.(Ch),'xCorr')
        fretSeg.(Ch).xCorr    = fretSeg.(Ch).xCorr(~filter,:);
        fretSeg.(Ch).yCorr    = fretSeg.(Ch).yCorr(~filter,:);
    end
    fretSeg.(Ch).int          = fretSeg.(Ch).int(~filter,:);
    fretSeg.(Ch).snr          = fretSeg.(Ch).snr(~filter,:);
    
    % Delete empty fields 
    fretSeg.(Ch).traceMetadata(filter) = [];
    
    % Set time
    fretSeg.(Ch).time = fretTraces.(Ch).time;
    
    % Reorder fields 
    fretSeg.(Ch) = orderfields(fretSeg.(Ch),fretTraces.(Ch));
end

%% Save 
    fretTraces = fretSeg; %#ok<NASGU>
    save(outfile,'fretTraces');

end




