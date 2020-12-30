function rejList=selTraces( inputfile, selList, isID, outfileSel, outfileRej)
%--------------------------------------------------------------------------
%
% selTraces.m:
%   The function selects a subset of traces from a fretTraces structure 
%   specified by a selection list.
% 
% Syntax:  
%   rejList=selTraces( inputfile, selList, isID, outfileSel, outfileRej)
% 
% Inputs:
%   1. inputfile  - A fretTraces file, using the naming convention 
%                   'singleTraces*.traces'. 
%   2. selList    - List of integers, e.g. [1 2 5 10, ...], stating which 
%                   traces need to be selected.  
%   3. isID       - Logical variable indicating whether the integer numbers 
%                   in selList are trace indices (isID = 0) or real trace
%                   identifiers (isID = 1).
%   4. outfileSel - Path and name of the output file containing the 
%                   selected data.  
%   5. outfileRej - Path and name of the output file containing the 
%                   rejected data.  
% 
% Outputs:
%   1. rejList    - List of integers, indicating the traces that were 
%                   rejected.  
%   2. outfileSel - A fretTraces*.traces file containing the selected data 
%                   is saved in the location specified by the user. 
%   3. outfileRej - A fretTraces*.traces file containing the rejected data 
%                   is saved in the location specified by the user.
%
% See also: 
%   scriptGetFretTraces.m
%
% Authors: 
%   - P.G. May 2011
%   - P.G. Feb 2015
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Script testing
% clear;
% inputfile='singleTracesCh2_stc.traces';
% selList=1:10;
% isID=0;
% outfileSel='singleTracesTest_stcs.traces';
% outfileRej='singleTracesTest_stcr.traces';

%% Load data 
data    = loadTracesCell(inputfile);
nTraces = size(data.traceMetadata,2);
time    = data.time;
xCor    = data.x;
yCor    = data.y;
I       = data.int;
% Check fields in structure array data 
if isfield(data,'snr')
    snr     = data.snr;
end

if isfield(data,'locInt')
    locInt  = data.locInt;
end

if isfield(data,'locBgd')
    locBgd  = data.locBgd;
end

if isfield(data,'xCorr')
    xMapped = data.xCorr;
    yMapped = data.yCorr;
end

% Load ids
ids        = cell(nTraces,1);
traceLen   = cell(nTraces,1);
lenBackset = cell(nTraces,1);
startOfTrace = cell(nTraces,1);
endOfTrace   = cell(nTraces,1);
lenBaseline     = cell(nTraces,1);

for i=1:nTraces
    ids{i}        = data.traceMetadata(i).ids;
    if isfield(data.traceMetadata,'lenBackset')
        lenBackset{i}   = data.traceMetadata(i).lenBackset;
        startOfTrace{i} = data.traceMetadata(i).startOfTrace;
        endOfTrace{i}   = data.traceMetadata(i).endOfTrace;
        lenBaseline{i}  = data.traceMetadata(i).lenBaseline;
        traceLen{i}     = data.traceMetadata(i).traceLen;
    end
end

%%
nSel=size(selList,2);
switch isID
    case 0
        % Store selected data
        xCorSel       = xCor(selList',:);
        yCorSel       = yCor(selList',:);
        intSel        = I(selList',:);
        idsSel        = ids(selList);
        
        if isfield(data,'snr')
            snrSel          = snr(selList',:);
        end
        
        if isfield(data,'locBgd')
            locBgdSel       = locBgd(selList',:);
        end
        
        if isfield(data,'locInt')
            locIntSel       = locInt(selList',:);
        end
        
        if isfield(data,'xCorr')
            xMappedSel      = xMapped(selList',:);
            yMappedSel      = yMapped(selList',:);
        end
        
        if isfield(data.traceMetadata,'lenBackset')
            lenBacksetSel   = lenBackset(selList);
            startOfTraceSel = startOfTrace(selList);
            endOfTraceSel   = endOfTrace(selList);
            lenBaselineSel  = lenBaseline(selList);
            traceLenSel     = traceLen(selList);
        end
        
        size(selList);
      
        % Store rejected data
        rejList=1:nTraces;
        for i=1:nSel
            rejList(selList(i))=0;
        end
        rejList=rejList(rejList>0);
        size(rejList);
        xCorRej       = xCor(rejList',:);
        yCorRej       = yCor(rejList',:);
        intRej        = I(rejList',:);
        idsRej        = ids(rejList);
        
        if isfield(data,'snr')
            snrRej          = snr(rejList',:);
        end
        
        if isfield(data,'locInt')
            locIntRej       = locInt(rejList',:);
        end
        
        if isfield(data,'locBgd')
            locBgdRej       = locBgd(rejList',:);
        end
        
        if isfield(data,'xCorr')
            xMappedRej      = xMapped(rejList',:);
            yMappedRej      = yMapped(rejList',:);
        end
        
        if isfield(data.traceMetadata,'lenBackset')
            lenBacksetRej   = lenBackset(rejList);
            startOfTraceRej = startOfTrace(rejList);
            endOfTraceRej   = endOfTrace(rejList);
            lenBaselineRej  = lenBaseline(rejList); 
            traceLenRej     = traceLen(rejList);
        end
        
    case 1
        % Read IDs from the inputfile
        idNumber=zeros(1,nTraces);
        for i=1:nTraces
            tmpId=textscan(ids{i}, '%*s %d','delimiter','_');
            idNumber(i)=cell2mat(tmpId);
        end
        
        % Compare old and new IDs
        existIDs=zeros(1,nTraces);
        for i=1:nSel
            idIndex=idNumber==selList(i);
            existIDs=existIDs+idIndex;
        end
        idsSelIdx=find(existIDs==1);
        idsRejIdx=find(existIDs==0);
        
        % Store selected data
        xCorSel       = xCor(idsSelIdx',:);
        yCorSel       = yCor(idsSelIdx',:);
        intSel        = I(idsSelIdx',:);
        idsSel        = ids(idsSelIdx);
        
        if isfield(data,'snr')
            snrSel          = snr(idsSelIdx',:);
        end
        
        if isfield(data,'locInt')
            locIntSel       = locInt(idsSelIdx',:);
        end
        
        if isfield(data,'locBgd')
            locBgdSel       = locBgd(idsSelIdx',:);
        end
        
        if isfield(data,'xCorr')
            xMappedSel      = xMapped(idsSelIdx',:);
            yMappedSel      = yMapped(idsSelIdx',:);
        end
        
        if isfield(data.traceMetadata,'lenBackset')
            lenBacksetSel   = lenBackset(idsSelIdx);
            startOfTraceSel = startOfTrace(idsSelIdx);
            endOfTraceSel   = endOfTrace(idsSelIdx);
            lenBaselineSel  = lenBaseline(idsSelIdx);
            traceLenSel     = traceLen(idsSelIdx);
        end
       
        % Store rejected data
        xCorRej       = xCor(idsRejIdx',:);
        yCorRej       = yCor(idsRejIdx',:);
        intRej        = I(idsRejIdx',:);
        idsRej        = ids(idsRejIdx);
        
        if isfield(data,'snr')
            snrRej          = snr(idsRejIdx',:);
        end
        
        if isfield(data,'locInt')
            locIntRej       = locInt(idsRejIdx',:);
        end
        
        if isfield(data,'locBgd')
            locBgdRej       = locBgd(idsRejIdx',:);
        end
        
        if isfield(data,'xCorr')
            xMappedRej      = xMapped(idsRejIdx',:);
            yMappedRej      = yMapped(idsRejIdx',:);
        end
        
        if isfield(data.traceMetadata,'lenBackset')
            lenBacksetRej   = lenBackset(idsRejIdx);
            startOfTraceRej = startOfTrace(idsRejIdx);
            endOfTraceRej   = endOfTrace(idsRejIdx);
            lenBaselineRej  = lenBaseline(idsRejIdx);
            traceLenRej     = traceLen(idsRejIdx);
        end
        rejList=idsRej;
end
   
%% Save the modified traces files 
if nargin<5
    dataSel.x   = xCorSel;
    dataSel.y   = yCorSel;
    dataSel.int = intSel;
    if isfield(data,'snr')
        dataSel.channelNames = {'x','y','int','snr','locBgd'};
        dataSel.snr    = snrSel;
        dataSel.locBgd = locBgdSel;
    end
    
    if isfield(data,'locInt')
        dataSel.channelNames = {'x','y','int','snr','locInt','locBgd'};
        dataSel.locInt = locIntSel;
    end
    
    if isfield(data,'xCorr')&& ~isfield(data,'locInt')
        dataSel.channelNames = {'x','y','xCorr','yCorr','int','snr','locBgd'};
        dataSel.xCorr = xMappedSel;
        dataSel.yCorr = yMappedSel;
    end
    
    if isfield(data,'xCorr')&& isfield(data,'locInt')
        dataSel.channelNames = {'x','y','xCorr','yCorr','int','snr','locInt','locBgd'};
        dataSel.xCorr = xMappedSel;
        dataSel.yCorr = yMappedSel;
    end
    
    dataSel.time= time;   
    for i=1:numel(selList)
        dataSel.traceMetadata(i).ids = idsSel{i};
        if isfield(data.traceMetadata,'lenBackset')
            dataSel.traceMetadata(i).lenBackset   = lenBacksetSel{i};
            dataSel.traceMetadata(i).startOfTrace = startOfTraceSel{i};
            dataSel.traceMetadata(i).endOfTrace   = endOfTraceSel{i};
            dataSel.traceMetadata(i).lenBaseline  = lenBaselineSel{i};
            dataSel.traceMetadata(i).traceLen     = traceLenSel{i};
        end
    end
    saveTracesCell(outfileSel,'traces',dataSel);
    
elseif nargin==5
    dataSel.x   = xCorSel;
    dataSel.y   = yCorSel;
    dataSel.int = intSel;
    if isfield(data,'snr')
        dataSel.channelNames = {'x','y','int','snr','locBgd'};
        dataSel.snr    = snrSel;
        dataSel.locBgd = locBgdSel;
    end
    
    if isfield(data,'locInt')
        dataSel.channelNames = {'x','y','int','snr','locInt','locBgd'};
        dataSel.locInt = locIntSel;
    end
    
    if isfield(data,'xCorr') && ~isfield(data,'locInt')
        dataSel.channelNames = {'x','y','xCorr','yCorr','int','snr','locBgd'};
        dataSel.xCorr = xMappedSel;
        dataSel.yCorr = yMappedSel;
    end
    
    if isfield(data,'xCorr') && isfield(data,'locInt')
        dataSel.channelNames = {'x','y','xCorr','yCorr','int','snr','locInt','locBgd'};
        dataSel.xCorr = xMappedSel;
        dataSel.yCorr = yMappedSel;
    end
    
    dataSel.time= time;   
    for i=1:numel(selList)
        dataSel.traceMetadata(i).ids = idsSel{i};
        if isfield(data.traceMetadata,'lenBackset')
            dataSel.traceMetadata(i).lenBackset   = lenBacksetSel{i};
            dataSel.traceMetadata(i).startOfTrace = startOfTraceSel{i};
            dataSel.traceMetadata(i).endOfTrace   = endOfTraceSel{i};
            dataSel.traceMetadata(i).lenBaseline  = lenBaselineSel{i};
            dataSel.traceMetadata(i).traceLen     = traceLenSel{i};
        end
    end
    saveTracesCell(outfileSel,'traces',dataSel);
    
    dataRej.x   = xCorRej;
    dataRej.y   = yCorRej;
    dataRej.int = intRej;
    if isfield(data,'snr')
        dataRej.channelNames = {'x','y','int','snr','locBgd'};
        dataRej.snr    = snrRej;
        dataRej.locBgd = locBgdRej;
    end
    
    if isfield(data,'locInt')
        dataRej.channelNames = {'x','y','int','snr','locInt','locBgd'};
        dataRej.locInt = locIntRej;
    end
    
     if isfield(data,'xCorr') && ~isfield(data,'locInt')
        dataRej.channelNames = {'x','y','xCorr','yCorr','int','snr','locBgd'};
        dataRej.xCorr = xMappedRej;
        dataRej.yCorr = yMappedRej;
    end
    
    if isfield(data,'xCorr') && isfield(data,'locInt')
        dataRej.channelNames = {'x','y','xCorr','yCorr','int','snr','locInt','locBgd'};
        dataRej.xCorr = xMappedRej;
        dataRej.yCorr = yMappedRej;
    end
 
    dataRej.time= time;   
    for i=1:numel(rejList)
        dataRej.traceMetadata(i).ids = idsRej{i};
        if isfield(data.traceMetadata,'lenBackset')
            dataRej.traceMetadata(i).lenBackset   = lenBacksetRej{i};
            dataRej.traceMetadata(i).startOfTrace = startOfTraceRej{i};
            dataRej.traceMetadata(i).endOfTrace   = endOfTraceRej{i};
            dataRej.traceMetadata(i).lenBaseline  = lenBaselineRej{i};
            dataRej.traceMetadata(i).traceLen     = traceLenRej{i};
        end
    end
    saveTracesCell(outfileRej,'traces',dataRej);
    
end
%%



