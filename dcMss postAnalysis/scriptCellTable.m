%--------------------------------------------------------------------------
%
% scriptCellTable:
%   Script that compiles the tracks and diffusion segments statistics of 
%   several cells into a single table. 
% 
% Description:
%   The script returns for each analyzed cell a statistic of the analyzed
%   tracks and diffusion segments which are classified by their mode of 
%   diffusion. 
%
% Syntax:  
%   scriptCellTable
% 
% Inputs: The user is prompeted to select the following files:  
%   1. CellDataPst.mat  - is an output file generated by the script 
%      lscriptFretTracesStat.mlx
%   2. CombinedStats.mat - is an output file generated by 
%      scriptGetSegResults.m
%   Note: Both files must be located in the same directory. Select both 
%   files (multi-select option is active), then press the Open button. 
%
% Outputs:
%   1. A table with 13 columbs and n (n = number of cells) rows in the 
%      command window (see Example) 
%   2. A workspace variable traceStatTable which contains the same data. 
%
% Example: 
%       cellNo    nAcc    nDon    nTraces    nPstFlt    nDiffFlt    nDiffSeg  ...
%       ______    ____    ____    _______    _______    ________    ________  ...
%   
%          1        85      65       1650       321         321         458   ...  
%          2       264     161       3211       733         733        1020   ... 
%          3       212     165       3031       606         606         763   ...    
%          4       133     100       2421       478         478         625   ...   
%          5       147     115       2346       380         380         569   ...   
%          6       148     112       2025       278         278         411   ...      
%    Total:6       989     718      14684      2796        2796        3846   ...
%
%      nSegImbl    nSegConf    nSegFree    nSegSupr    occurence     totalLength
%      ________    ________    ________    ________    __________    ___________
%   
%         55         111          275         17       {4�1 cell}         36944 
%        145         199          617         59       {4�1 cell}    1.0785e+05 
%        183         151          392         37       {4�1 cell}         60460 
%         68         125          406         26       {4�1 cell}         55389 
%        138         120          278         33       {4�1 cell}         65026 
%         58          64          259         30       {4�1 cell}         46997 
%  Total: 647         770         2227        202       []            3.7267e+05
%
% 
% 
% See also: 
%   lscriptFretTracesStat.mlx,  scriptGetSegResults.m
%
% Author: 
%   - P.G. Jul 2018
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------
%% Select Files

clear;
try 
    [file,pathName] = uigetfile('*.mat',...
                            'Multi Select a cellData** & CombinedStats file',...
                            'MultiSelect','on');
catch
    % No file seleceted
    return
end

% Load Data 
try
    load([pathName filesep file{1}],'cellData');
    load([pathName filesep file{2}]);
catch 
    load([pathName filesep file{2}],'cellData');
    load([pathName filesep file{1}]);
end
    
%% Filter Statistics
statData = struct('cellNo',[],...
                  'nAcc',[],...
                  'nDon',[],...
                  'nTraces',[],...
                  'nPstFlt',[],...
                  'nDiffFlt',[],...
                  'nDiffSeg',[],...
                  'nSegImbl',[],...
                  'nSegConf',[],...
                  'nSegFree',[],...
                  'nSegSupr',[],...
                  'occurence',[],...
                  'totalLength',[]);
             
statDataSegFN = {'cellNo',...
                 'nAcc',...
                 'nDon',...
                 'nTraces',...
                 'nPstFlt',...
                 'nDiffFlt',...
                 'nDiffSeg',...
                 'nSegImbl',...
                 'nSegConf',...
                 'nSegFree',...
                 'nSegSupr',...
                 'occurence',...
                 'totalLength'};
             
nCells = size(cellData,2);
cellNo = zeros(nCells+1,1);
%%
for i=1:nCells
    
    % Read Cell Number
    txt = cellData(i).acc;
    cellNoStr = textscan(txt(2:3),'%2c');
    statData(i).cellNo = str2double(cellNoStr{1});
    

    % Read Acc and Donor Density
    if size(cellData(i).density,2)>1
        statData(i).nAcc = cellData(i).density(4);
        statData(i).nDon = cellData(i).density(1);
    else
        statData(i).nAcc = cellData(i).density(1);
        statData(i).nDon = cellData(i).density(1);
    end

    % Read nTraces
    statData(i).nTraces = cellData(i).nTraces;

    % Read nPstFlt
    if ~isempty(cellData(i).nPstFlt)
        statData(i).nPstFlt = cellData(i).nPstFlt;
    else
        statData(i).nPstFlt = 0;
    end

    % Get number of traces that are analaysed by the MSS-Diff Analysis
    path    = cellData(i).path;
    dataset = cellData(i).acc;
    
    % Load SegResultsFinal into workspace
    try
        filename   = [ path filesep dataset filesep 'SegResultsFinal.mat'];
        load(filename);
    catch
        newPath = strrep(path,'P:','F:\Disk E');
        filename   = [ newPath filesep dataset filesep 'SegResultsFinal.mat'];
        load(filename);
    end

    if ~isempty(segResultsFinal)
        statData(i).nDiffFlt = size(unique(segResultsFinal(:,23)),1);
        % Get Segment numbers for each state
        n=4*i-3;
        statData(i).nSegImbl = CombinedStats(n+0,2);
        statData(i).nSegConf = CombinedStats(n+1,2);
        statData(i).nSegFree = CombinedStats(n+2,2);
        statData(i).nSegSupr = CombinedStats(n+3,2);
        statData(i).nDiffSeg = sum(CombinedStats(n:n+3,2));

        idxClass0=segResultsFinal(:,3)==0;
        idxClass1=segResultsFinal(:,3)==1;
        idxClass2=segResultsFinal(:,3)==2;
        idxClass3=segResultsFinal(:,3)==3;

        Class0=segResultsFinal(idxClass0,:);
        Class1=segResultsFinal(idxClass1,:);
        Class2=segResultsFinal(idxClass2,:);
        Class3=segResultsFinal(idxClass3,:);

        % Occurence in % 
        statData(i).occurence   = {Class0(:,24); Class1(:,24);Class2(:,24);Class3(:,24)};
        statData(i).totalLength = sum(segResultsFinal(:,24));

    else
        for j=6:numel(statDataSegFN)
            fn=statDataSegFN{j};
            statData(i).(fn) = 0;
        end
        
    end

end

% Calculate Total Numbers
statData(nCells+1).cellNo = nCells;
for j=2:numel(statDataSegFN)
    if j==12, continue, end
    fn=statDataSegFN{j};
    statData(nCells+1).(fn) = sum([statData(1:nCells).(fn)].');
end

% Make Table
traceStatTable=struct2table(statData);
head(traceStatTable,nCells+1)


%% Calc Occurence (all together)
totalTotalLength=statData(nCells+1).totalLength;
segLenS0 = 0;
segLenS1 = 0;
segLenS2 = 0;
segLenS3 = 0;

segLenTmpS0 = [];segLenTmpS1 = [];segLenTmpS2 = [];segLenTmpS3 = [];
for j=1:nCells
    if iscell(statData(j).occurence)
        segLenS0    = vertcat(segLenTmpS0,statData(j).occurence{1});
        segLenTmpS0 = segLenS0;

        segLenS1    = vertcat(segLenTmpS1,statData(j).occurence{2});
        segLenTmpS1 = segLenS1;

        segLenS2    = vertcat(segLenTmpS2,statData(j).occurence{3});
        segLenTmpS2 = segLenS2;

        segLenS3    = vertcat(segLenTmpS3,statData(j).occurence{4});
        segLenTmpS3 = segLenS3;
    end
end

fract_Imob = nansum(segLenS0)/totalTotalLength;
fract_Conf = nansum(segLenS1)/totalTotalLength;
fract_Free = nansum(segLenS2)/totalTotalLength;
fract_Supr = nansum(segLenS3)/totalTotalLength;

chk=sum([fract_Imob, fract_Conf, fract_Free, fract_Supr]);

%% Clear workspace
clearvars -except traceStatTable chk

%%
% cellNo = 7;
% diffState = 3;
% totalDiffLen = statData(cellNo).totalLength;
% lenDiffTrack = sum(statData(cellNo).occurence{diffState});
% fraction = lenDiffTrack/totalDiffLen 







