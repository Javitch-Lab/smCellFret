%--------------------------------------------------------------------------
%
% scriptTotalFilter.m:
%   Applies a total intensity filter to a fretTraces structure variable 
%   saved in a fretTraces*.mat file. 
% 
% Description:
%   Running scriptTotalFilter.m provides the possibility to use only FRET 
%   traces with total intensities in a predefined interval in the further
%   downstream data analysis. This step is particularly useful for cleaning 
%   up trace populations that have not yet been processed by other filters. 
% 
% Syntax:  
%   scriptTotalFilter.m
% 
% Inputs:
%   1. The script prompts the user to select one or more fretTraces*.mat 
%      files. 
%   2. In the code section 'Filter data', update the value for the mean
%      total intensity and its standard deviation. The default values for 
%      the mean and standard deviation here were calculated from time 
%      traces measured in fixed cells.
% 
% Outputs:
%   The filtered structure variable fretTraces is saved at the same 
%   location and under the same file name as the input file with the 
%   additional extension *_flt.mat. 
%
% Authors: 
%   - P.G. Jun 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% Get files from user
% -------------------------------------------------------------------------

% Get FILES
clear;
matfiles = getFiles('*.mat','Select multiple fretTraces*.mat files (Ch2 Directories)');

% Check if files array is empty
if any(cellfun('isempty',matfiles))== true; return;end

% Show selection
for z=1:numel(matfiles)
    [path,name,ext]=fileparts(matfiles{z});
    disp(path);
    disp([name ext]);
end

%% -------------------------------------------------------------------------
% Filter data
% -------------------------------------------------------------------------

% Total Intensity Filter Parameters 
fixedMean = 457.91;
fixedStdv = 218.144/2;
lwLim = fixedMean-2*fixedStdv;
upLim = fixedMean+3*fixedStdv;

% 
nCells = numel(matfiles);
data    =struct([]);
h = waitbar(0,'Load Data ...');

for j=1:nCells 
    clear tmp;
    waitbar(j / nCells)
    
    % load traces file
    tracesFile = matfiles{j};
    load(tracesFile);
    totalInt = fretTraces.Fret.total;
    [nTraces,nFrames]=size(totalInt);
    
    % save filename
    data(j).file = tracesFile;

    % logical Idl Array 
    % idl = fretTraces.Fret.idlTotal; 
    idl = false(nTraces,nFrames);
    donSOT= [fretTraces.Ch2.traceMetadata.startOfTrace].'; %unit is frames
    donEOT= [fretTraces.Ch2.traceMetadata.endOfTrace].'; %unit is frames
        for i=1:nTraces
            idl(i,donSOT(i):donEOT(i))=true;
        end
    
    % Calc Total intensity where trace is idealized 
    totalInt_idl2D = totalInt.*idl;
    totalInt_idl2D(totalInt_idl2D==false)=NaN;
    meanTotalPerTrace = nanmean(totalInt_idl2D,2);
    data(j).total = meanTotalPerTrace;
    
    % Get Filter Idx
    idx = meanTotalPerTrace > lwLim & meanTotalPerTrace < upLim;
    data(j).flt = idx;
    data(j).nMol = sum(idx);
    
    % Apply Filter to data
    fretTraces  = applyFilterToData(fretTraces,idx);
    data(j).fltTraces = fretTraces;
    
    % Save Traces to File
    [pathName, fileName, ext]= fileparts(tracesFile);
    newExt = strrep(ext,'.mat','_flt.mat');
    save([pathName filesep fileName newExt],'fretTraces')
    
    
end
close(h)

%% Disp Number of Molecules
nMol = sum([data.nMol].');
disp(['nMol = ' num2str(nMol)]);


%% =======================================================================%
%                                                                         %
%                        Local Callback Functions                         %
%                                                                         %
%=========================================================================%


%% filter matrix array

function fltTraces=applyFilterToData(fretData,idxFilter)

    fltTraces=[];
    
    % Channel Names (Acc: Ch1, Don: Ch2)  
    chName = {'Ch1','Ch2'};
    
    % Metadata Field Names
    metaFN = fieldnames(fretData.Ch1.traceMetadata);
    if ~isempty(fretData.Fret.traceStatSel)
        traceStatFN = fieldnames(fretData.Fret.traceStatSel);
    end
         
    % Update Ch1(Acc) and Ch2(Don)
    nTraces=size(fretData.Ch1.x);   
    for c=1:2
        
        % Update data
        Ch=chName{c};
        fltTraces.(Ch).time         = fretData.(Ch).time;
        fltTraces.(Ch).x            = fretData.(Ch).x(idxFilter,:);
        fltTraces.(Ch).y            = fretData.(Ch).y(idxFilter,:);
        if isfield(fretData.(Ch),'xCorr')
            fltTraces.(Ch).xCorr    = fretData.(Ch).xCorr(idxFilter,:);
            fltTraces.(Ch).yCorr    = fretData.(Ch).yCorr(idxFilter,:);
        end
        fltTraces.(Ch).int          = fretData.(Ch).int(idxFilter,:);
        fltTraces.(Ch).snr          = fretData.(Ch).snr(idxFilter,:);

        % Update Metadata for each trace 
        j=0;
        for i=1:nTraces
            
            % Apply Filter
            if idxFilter(i)==0; continue;end
            j=j+1;
            
            %Update Metadata
            for k=1:numel(metaFN)
                fn=metaFN{k};
                if isfield(fretData.(Ch).traceMetadata,fn)
                    fltTraces.(Ch).traceMetadata(j).(fn) =...
                    fretData.(Ch).traceMetadata(i).(fn);
                end
            end
            
        end % End Update Metadata 
    end %End Update Ch1,Ch2
    
    % Update FRET fields
    fltTraces.Fret.idlTotal = fretData.Fret.idlTotal(idxFilter,:);
    fltTraces.Fret.total    = fretData.Fret.total(idxFilter,:);
    
    % Update Fret Trace Stat for each trace 
    m=0;n=0; 
    for i=1:nTraces
        % Apply Filter
        if idxFilter(i)==0
            m=m+1;
%             %Update traceStat (rejected)
%             for k=1:numel(traceStatFN)
%                 fn=traceStatFN{k};
%                 if isfield(fretData.Fret.traceStatSel,fn)
%                     fltTraces.Fret.traceStatRej(m).(fn) =...
%                     fretData.Fret.traceStatSel(i).(fn);
%                 end
%             end
            
        else
            n=n+1;
            %Update traceStat (selected)
            if ~isempty(fretData.Fret.traceStatSel)
                for k=1:numel(traceStatFN)
                    fn=traceStatFN{k};
                    if isfield(fretData.Fret.traceStatSel,fn)
                        fltTraces.Fret.traceStatSel(n).(fn) =...
                        fretData.Fret.traceStatSel(i).(fn);
                    end
                end
            end
        end

    end % End Update Metadata 
    
    
    
end %End Function
%% End
