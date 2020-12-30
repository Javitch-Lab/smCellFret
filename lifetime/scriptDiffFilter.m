%--------------------------------------------------------------------------
%
% scriptDiffFilter.m:
%   The script applies a diffusion filter to fretTracesPst.mat files. 
% 
% Description:
%   The script uses the results of the motion classification to extract 
%   from fretTracesPst.mat those traces associated with entirely immobile,
%   entirely restricted, entirely free and entirely directed diffusion. 
% 
% Syntax:  
%   scriptDiffFilter
% 
% Inputs:
%   The script prompts the user to select pairs of #Ch1 and #CH2 
%   directories. #Ch1 directories must contain data files of type 
%   diffFilter.mat, while #Ch2 directories must contain data files of type 
%   fretTracesPst.mat (fretTracesPst.mat files are gnerated with 
%   scriptFretFilter.m using the FDT filter condition).
% 
% Outputs:
%   The results are saved in each of the selected #Ch2 folders under the
%   following file names: 
%   1. fretTracesImob.mat - FRET traces data associated with entirely 
%      immobile diffusion
%   2. fretTracesConf.mat - FRET traces data associated with entirely
%      confined diffusion
%   3. fretTracesFree.mat - FRET traces data associated with entirely free 
%      diffusion
%   4. fretTracesSupr.mat - FRET traces data associated with entirely 
%      directed diffusion
% 
% See also: 
%   scriptFretFilter.m, scriptDiffLT.m
%
% Authors: 
%   - P.G. Mar 2019
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------





%% Get Directories (Ch1 and Ch2 directories required)

clear
dirArray = getDirs('select multiple pairs of Ch1  A N D  Ch2 directories');
nExperiments = length(dirArray(1,:));
if nExperiments==0,return,end

%% Initialize Parameters
diffFilter = {'isImob','isConf','isFree','isSupr'};

%% Loop through all experiments

% Set current Directopry
for j=1:nExperiments 
    
    % Set current Directopry
    basepath = dirArray{2,j};
    dataSet  = dirArray{1,j};
    currentDir = [basepath filesep dataSet filesep];
    
    disp(dataSet);
    
    % Load Filter Index (in Dir Ch1) 
    if contains(dataSet,'Ch1')
        logicalFilter = load([currentDir 'diffFilter.mat']);

    % Load post filtered fretTraces (in Dir Ch2) 
    elseif contains(dataSet,'Ch2')
        
        data=load([currentDir 'fretTracesPst.mat']);
   
        %------------------------------------------------------------------
        %--- Apply Filter and save fretTraces
        %------------------------------------------------------------------
        disp('Run Diffusion Filters...');
        
        for s = 1:numel(diffFilter)
            state = diffFilter{s};
            idxFilter  = [logicalFilter.diffFilter.(state)].'; 
            fretTraces = applyFilterToFrettraces(data.fretTraces,idxFilter); 
            fretTraces.Fret.pstFltVal = data.fretTraces.Fret.pstFltVal;
            % Save data
            outfile=strrep('fretTracesPst.mat','Pst',state(3:6));
            save([currentDir outfile],'fretTraces');
        end
    end
    
    
end % End Experiment Loop
disp('Finished'); 



%% =======================================================================%
%                                                                         %
%                        Local Callback Functions                         %
%                                                                         %
%=========================================================================%


%-------------------------------------------------------------------------%
%                                                                         %
% Apply Filter To fretTracesPst.mat                                       %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

function fltTraces=applyFilterToFrettraces(fretTraces,idxFilter) 

    fltTraces=[];

    % Channel Names (Acc: Ch1, Don: Ch2)  
    chName = {'Ch1','Ch2'};

    % Metadata Field Names
    metaFN = fieldnames(fretTraces.Ch1.traceMetadata);
    %traceStatRejFN = fieldnames(fretTraces.Fret.traceStatRej);
    traceStatSelFN = fieldnames(fretTraces.Fret.traceStatSel);

    % Update Ch1(Acc) and Ch2(Don)
    nTraces=size(fretTraces.Ch1.x);   
    for c=1:2

        % Update data
        Ch=chName{c};
        fltTraces.(Ch).time         = fretTraces.(Ch).time;
        fltTraces.(Ch).x            = fretTraces.(Ch).x(idxFilter,:);
        fltTraces.(Ch).y            = fretTraces.(Ch).y(idxFilter,:);
        if isfield(fretTraces.(Ch),'xCorr')
            fltTraces.(Ch).xCorr    = fretTraces.(Ch).xCorr(idxFilter,:);
            fltTraces.(Ch).yCorr    = fretTraces.(Ch).yCorr(idxFilter,:);
        end
        fltTraces.(Ch).int          = fretTraces.(Ch).int(idxFilter,:);
        fltTraces.(Ch).snr          = fretTraces.(Ch).snr(idxFilter,:);

        % Update Metadata for each trace 
        m=0;
        for i=1:nTraces

            % Apply Filter
            if idxFilter(i)==0; continue;end
            m=m+1;

            %Update Metadata
            for k=1:numel(metaFN)
                fn=metaFN{k};
                if isfield(fretTraces.(Ch).traceMetadata,fn)
                    fltTraces.(Ch).traceMetadata(m).(fn) =...
                    fretTraces.(Ch).traceMetadata(i).(fn);
                end
            end

        end % End Update Metadata 
    end %End Update Ch1,Ch2

    % Update FRET fields
    fltTraces.Fret.idlTotal = fretTraces.Fret.idlTotal(idxFilter,:);
    fltTraces.Fret.total    = fretTraces.Fret.total(idxFilter,:);

    % Update Fret Trace Stat for each trace 
    n=0; 
    for i=1:nTraces
        % Apply Filter
        if idxFilter(i)==0; continue; end
        n=n+1;
        
        %Update traceStat (selected)
        for k=1:numel(traceStatSelFN)
            fn=traceStatSelFN{k};
            if isfield(fretTraces.Fret.traceStatSel,fn)
                fltTraces.Fret.traceStatSel(n).(fn) =...
                fretTraces.Fret.traceStatSel(i).(fn);
            end
        end
        
    end % End Update Metadata 
end