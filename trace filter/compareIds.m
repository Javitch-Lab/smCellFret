%--------------------------------------------------------------------------
%
% compareIds.m:
%   Script that filters a fretTraces structure saved as a .mat file by 
%   trace Ids. 
% 
% Description:
%   Filtering a 'fretTracesDif_*.mat' file by trace IDs. The script 
%   compares the trace IDs in 'fretTracesDif_*.mat' with a pool of best IDs 
%   (trace Ids that are e.g. manually selected). All traces in 
%   'fretTracesDif*.mat' that match the trace Ids in the pool are selected. 
% 
% Syntax:  
%   compareIds
% 
% Inputs:
%   The script is prompting the user to select:
%   1. fretTraces.mat file(s) - One or more fretTraces*.mat files with the 
%      track IDs to be selected by the filter. Use *.bestFret.mat and / or 
%      *.allFret.mat files to create a pool of IDs.
%   2. fretTraces.mat file - The file that is filtered by trace IDs. Only 
%      traces with Ids that are present in the ID pool will be selected.
% 
% Outputs:
%   The script prompts the user to enter a file name under which the 
%   filtered data should be saved.
%
% Other m-files required: 
%   Subfunctions: uses a local call back function to filter data
%
% See also: 
%   scriptGetfretTraces, scriptFretFilter.m, scriptFretTraces2Seg.m 
%
% Authors: 
%   - P.G. Aug 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Set initial Parameters
clear; 

%% ------------------------------------------------------------------------
% Get files from user
% -------------------------------------------------------------------------

%--- Select data set 1
mfiles = getFiles('*.mat','SELECT FILES TO GET A POOL OF IDs (Filter IDs)');

% Check if files array is empty
if any(cellfun('isempty',mfiles))== true; return;end

% Load selection to structure data1 
data = struct([]); idsPool=[];
for z=1:numel(mfiles)
    
    % Load data
    matFile = mfiles{z};
    [path,fn,ext]  = fileparts(matFile);
    data(z).path = path;
    data(z).file = fn;
    
    load(matFile);
    ids = {fretTraces.Ch1.traceMetadata.ids}.';
    
    % read Ids as real numbers
    nTraces = size(ids,1);
    idNumber=zeros(1,nTraces);
    for k=1:nTraces
        tmpId=textscan(ids{k}, '%*s %d','delimiter','_');
        idNumber(k)=cell2mat(tmpId);
    end
    data(z).ids = idNumber';
    data(z).nTraces = nTraces;
    
    % Concatenate Ids
    idsPool = [idsPool, idNumber]; 
    
end
% this fretTraces file will not be used any more
clear fretTraces ids
    
%% --- Select data set 2
[file,path] = uigetfile('*.mat','SELECT A FILE THAT WILL BE FILTERED BY IDs');

% Check if files array is empty
if file == 0; return;end

%% Load selection to structure data2 

% Load data
data(z+1).path = path;
data(z+1).file = fn;

load([path filesep file]);
ids = {fretTraces.Ch1.traceMetadata.ids}.';

% read Ids as real numbers
nTraces = size(ids,1);
fileIds=zeros(1,nTraces);
for k=1:nTraces
    tmpId=textscan(ids{k}, '%*s %d','delimiter','_');
    fileIds(k)=cell2mat(tmpId);
end
data(z+1).ids = fileIds';
data(z+1).nTraces = nTraces;


%% compare Ids
% Lia = ismember(A,B) returns an array containing logical 1 (true) where 
% the data in A is found in B. Elsewhere, the array contains logical 0 (false).
logicalFileIds=ismember(fileIds,idsPool);
nSelTraces=sum(logicalFileIds);

%% ------------------------------------------------------------------------
 %--- Apply Filter and save fretTraces / Tracking.mat
 %-------------------------------------------------------------------------
%  if nSelTraces > 0
%     newFile = strrep(file,'.mat','_allBst.mat');
%     fn=[path newFile];
%     disp(['nSelTraces = ' num2str(nSelTraces)]);
%     disp(['Saving ...' fn]);
%     
%     fretTraces = applyFilterToData(fretTraces,logicalFileIds);
%     
%     save(fn,'-mat','fretTraces','-v7.3');
%  else
%      disp('no traces were seletced')
%  end

if nSelTraces > 0
    name = strrep(file,'.mat','_entireImb.mat');      %txt default answer
    newFileName  = inputdlg({'Enter File Name:'},...  %prompt
                             'Save File As: ',...     %dlg titel
                              [1 50],...              %width of line
                              {name});                %default answer

    if isempty(newFileName)
        return
    else
        fn=[path newFileName{1}];
        disp(['nSelTraces = ' num2str(nSelTraces)]);
        disp(['Saving ...' fn]);

        fretTraces = applyFilterToData(fretTraces,logicalFileIds);

        save(fn,'-mat','fretTraces','-v7.3');
    end
else
     disp('no traces were seletced')
end

%% End


%-------------------------------------------------------------------------%
%                                                                         %
% Apply Filter To Fret Data                                               %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%
function fltTraces=applyFilterToData(fretData,idxFilter)

    fltTraces=[];
    
    % Channel Names (Acc: Ch1, Don: Ch2)  
    chName = {'Ch1','Ch2'};
    
    % Metadata Field Names
    metaFN = fieldnames(fretData.Ch1.traceMetadata);  
    
    if isfield(fretData,'Fret') && isfield(fretData.Fret,'traceStatSel')
        traceStatFN = fieldnames(fretData.Fret.traceStatSel);
    else 
        traceStatFN =[];
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
    if isfield(fretData,'Fret')
        fltTraces.Fret.idlTotal = fretData.Fret.idlTotal(idxFilter,:);
        fltTraces.Fret.total    = fretData.Fret.total(idxFilter,:);
    end
    
    % Update Fret Trace Stat for each trace 
    if ~isempty(traceStatFN) 
        m=0;n=0; 
        for i=1:nTraces
            % Apply Filter
            if idxFilter(i)==0
                m=m+1;
                %Update traceStat (rejected)
                for k=1:numel(traceStatFN)
                    fn=traceStatFN{k};
                    if isfield(fretData.Fret.traceStatSel,fn)
                        fltTraces.Fret.traceStatRej(m).(fn) =...
                        fretData.Fret.traceStatSel(i).(fn);
                    end
                end

            else
                n=n+1;
                %Update traceStat (selected)
                for k=1:numel(traceStatFN)
                    fn=traceStatFN{k};
                    if isfield(fretData.Fret.traceStatSel,fn)
                        fltTraces.Fret.traceStatSel(n).(fn) =...
                        fretData.Fret.traceStatSel(i).(fn);
                    end
                end
            end

        end % End Update Metadata 
    end
    
    
    
end %End Function