%--------------------------------------------------------------------------
%
% getTraceInt:
%   Get intensities and lifetimes of fret traces objects 
% 
% Description:
%   The Script compiles intensity and lifetime data of fret traces objects 
%   (Spartan) into histograms. 
%
%   NOTE: Time traces, can be segmented after reading out only one specific
%   diffusion mode (here the mode of free diffusion). This means that the 
%   value of their amplitude (intensity or FRET ) is set to zero at the 
%   points where the diffusion mode is either confined, immobile or
%   directional.  
%     Calculation of the mean total intensity: For segmented time traces
%     with one diffusion mode, the total intensity is only calculated if
%     the signal amplitude is greater than zero.
%     The lifetime distributions are based on the total length of the time 
%     trace. For segmented time traces with only one diffusion mode, the
%     lifetime distributions are calculated based on the whole trace length, 
%     i.e. segments where the signal amplitude is zero are included in the 
%     lifetime calculation. TODO: parse out segments with zero signal 
%     amplitude. 
%     NOTE: Therefore, the lifetime analysis described above differs from 
%     the lifetime analysis where only traces that are entirely in one 
%     diffusion mode are analyzed. 
%
% Syntax:  
%   getTraceInt
% 
% Inputs:
%   The script prompts the user to select one or multiple Spatan traces 
%   objects e.g. fretTracesDif_all_SegFre01_sync.traces
% 
% Outputs:
%   1. The script generates the following workspace variables: 
%       a. histDataTotal - total intensity histogram array of N x data sets.
%       The last three columns contain a) the total intensity distribution, 
%       averaged over all data sets, b) the standard deviations and c) the
%       cumulative distribution.   
%       b. histDataAcc - acceptor intensity histogram array of N x data sets.
%       The last three columns contain a) the acceptor intensity distribution, 
%       averaged over all data sets, b) the standard deviations and c) the
%       cumulative distribution.  
%       c. histDataAccLt - Lifetime histogram array of N x datasets. The last
%       two columns contain a) the cumulative distribution over all datasets 
%       and b) the normalized cummulative distribution.    
%       d. histDataAccLtFlt - As in point 3, but with additional total 
%       intensity filter.  
%   2. The following figures will be generated: 
%       a. Array of N total intensity histograms
%       b. Array of N acceptor intensity histograms
%       c. Array of N acceptor lifetime histograms
%       d. Array of N acceptor total intensity filtered lifetime histograms
% 
% See also: 
%   loadTraces.m (Spartan), scriptFretTraces2Seg.m, getTraceSnr.m 
%
% Authors: 
%   - P.G. 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Set initial Parameters
clear;

%% -------------------------------------------------------------------------
% Get files from user
% -------------------------------------------------------------------------

% Get FILES
tfiles = getFiles('*.traces','Select *.traces files');

% Check if files array is empty
if any(cellfun('isempty',tfiles))== true; return;end

% Show selection
for z=1:numel(tfiles)
    [path,name,ext]=fileparts(tfiles{z});
    disp(path);
    disp([name ext]);
end

%% Get Intensities

nCells = numel(tfiles);
data    =struct([]);

% Initialize Intensity Histogram
edges = -200:100:3500;
nBins=numel(edges)-1;
histDataTotal=zeros(1,nCells+1);
histDataAcc=zeros(1,nCells+1);


% Initialize Intensity Histogram
edgesLt = 0:4:250;
nBinsLt=numel(edgesLt)-1;
histDataAccLt=zeros(1,nCells+1);

% Total Intensity Filter Parameters 
fixedMean = 457.91;
fixedStdv = 218.144/2;
lwLim = fixedMean-2*fixedStdv
upLim = fixedMean+3*fixedStdv


%% Analyze Data 

% Initialize figures and waitbar
figure(1); %total
figure(2); %accInt
figure(3); %accLT 
figure(4); %accLT filtered 
f = waitbar(0,'Load Data ...');

% Loop over selected experiment folders
for j=1:nCells
    clear tmp;
    waitbar(j / nCells,f)
    
    % load traces file
    tracesFile = tfiles{j};
    tmpData=loadTraces(tracesFile);
    totalInt = tmpData.total;
    accInt   = tmpData.acceptor;
    [nTraces,nFrames]=size(totalInt);
    
    % save filename
    data(j).file = tracesFile;
    
    %---
    %--- GET TOTAL INTENSITY
    %---
    
    % logical Idl Array 
    idlAcc = false(nTraces,nFrames);
    accSOT= [tmpData.traceMetadata.startOfTraceAcc].'; %unit is frames
    accEOT= [tmpData.traceMetadata.endOfTraceAcc].'; %unit is frames
    
    idlDon = false(nTraces,nFrames);
    donSOT= [tmpData.traceMetadata.startOfTraceDon].'; %unit is frames
    donEOT= [tmpData.traceMetadata.endOfTraceDon].'; %unit is frames
    
    for i=1:nTraces
        idlAcc(i,accSOT(i):accEOT(i))=true;
        idlDon(i,donSOT(i):donEOT(i))=true;
    end
    
    % Calc Total intensity per Trace over the region where trace is idealized / tracked 
    totalInt_idl2D = totalInt.*idlDon;
    % Don not include zeros in the calculation of the mean total intensity
    totalInt_idl2D(totalInt_idl2D==false)=NaN;
    meanTotalPerTrace = nanmean(totalInt_idl2D,2);
    data(j).total = meanTotalPerTrace;
    
    % Plot Total Intensity Histogramms
    figure(1)
    subplot(2,round(nCells/2),j)
    XTotal=data(j).total;
    histTotal(j) = histogram(XTotal,edges); 
    title('Total Intensity')
    % Offset the bin width by half the bin width to plot data at bin
    % centers
    histDataTotal(1:nBins,1) = (histTotal(j).BinEdges(1:nBins)+histTotal(j).BinWidth/2)';
    histDataTotal(1:nBins,j+1) = histTotal(j).Values';
    
    %---
    %--- GET ACCEPTOR INTENSITY
    %---
    
    % Calc Acc intensity per Trace over the region where trace is idealized  / tracked 
    accInt_idl2D = accInt.*idlAcc;
    accInt_idl2D(accInt_idl2D==false)=NaN;
    meanAccPerTrace = nanmean(accInt_idl2D,2);
    data(j).acc = meanAccPerTrace;
    
    % Number of Molecules
    data(j).nTraces = size(accInt,1);
    
    % Plot Acc Intensity Histogramms
    figure(2)
    subplot(2,round(nCells/2),j)
    XAcc=data(j).acc;
    histAcc(j) = histogram(XAcc,edges);
    histDataAcc(1:nBins,1) = (histAcc(j).BinEdges(1:nBins)+histAcc(j).BinWidth/2)';
    histDataAcc(1:nBins,j+1) = histAcc(j).Values';
    title('Acceptor Intensity')
    %---
    %--- GET ACCEPTOR LIFETIMES OF TOTAL INTENSITY FILTEERED DATA 
    %--- (S E G M E N T Lifetime)
    %---
    
    % Get lifetime (no total filter)
    figure(3)
    dt=tmpData.sampling/1000;
    data(j).accLt = dt*[tmpData.traceMetadata.traceLenAcc].';
    subplot(2,round(nCells/2),j)
    accLt=data(j).accLt;
    histAccLt(j) = histogram(accLt,edgesLt);
    histDataAccLt(1:nBinsLt,1) = (histAccLt(j).BinEdges(1:nBinsLt)+histAccLt(j).BinWidth/2)';
    histDataAccLt(1:nBinsLt,j+1) = histAccLt(j).Values';
    title('Acc LT (no total Intensity filter')
    
    % Get Filter Idx
    idx = meanTotalPerTrace > lwLim & meanTotalPerTrace < upLim;
    data(j).fltIdx = idx;
    data(j).nTracesFlt = sum(idx);

    % Apply Filter to data
    tracesData  = applyFilterToData(tmpData,idx);
    data(j).fltTraces = tracesData;
    
    % Get lifetime after applying total intensity filtered
    figure(4)
    data(j).accLtFlt = dt*[data(j).fltTraces.traceMetadata.traceLenAcc].';
    subplot(2,round(nCells/2),j)
    histAccLtFlt(j) = histogram(data(j).accLtFlt,edgesLt);
    histDataAccLtFlt(1:nBinsLt,1) = (histAccLtFlt(j).BinEdges(1:nBinsLt)+histAccLtFlt(j).BinWidth/2)';
    histDataAccLtFlt(1:nBinsLt,j+1) = histAccLtFlt(j).Values';
    title('Acc LT (with total Intensity filter')
    
end
close(f)

%% Disp Number of Molecules
nAllTraces = sum([data.nTraces].');
nFltTraces = sum([data.nTracesFlt].');
disp(['nAllTraces = ' num2str(nAllTraces)]);
disp(['nFltTraces = ' num2str(nFltTraces)]);

%%  Mean Hist
% total mean
histDataTotal(1:nBins,j+2) = nanmean(histDataTotal(:,2:nCells+1),2);
% total std
histDataTotal(1:nBins,j+3) = nanstd(histDataTotal(:,2:nCells+1),0,2);
% total cummulative
histDataTotal(1:nBins,j+4) = sum(histDataTotal(:,2:nCells+1),2);

% Acc mean
histDataAcc(1:nBins,j+2) = nanmean(histDataAcc(:,2:nCells+1),2);
% Acc std
histDataAcc(1:nBins,j+3) = nanstd(histDataAcc(:,2:nCells+1),0,2);
% Acc cummulative
histDataAcc(1:nBins,j+4) = sum(histDataAcc(:,2:nCells+1),2);

%--- Acc LT
% sum
histDataAccLt(1:nBinsLt,j+2) = sum(histDataAccLt(:,2:nCells+1),2);
% normalize [0 1]
histDataAccLt(1:nBinsLt,j+3) = histDataAccLt(1:nBinsLt,j+2)./histDataAccLt(1,j+2);
% for cum fitting
ltValuesFlt=[];
for i=1:size(data,2)
    ltValuesFlt = vertcat(ltValuesFlt,data(i).accLt);
end
ltValuesFlt=ltValuesFlt/dt;

% Acc LT Flt
% Sum 
histDataAccLtFlt(1:nBinsLt,j+2) = sum(histDataAccLtFlt(:,2:nCells+1),2);
% normalize [0 1]
histDataAccLtFlt(1:nBinsLt,j+3) = histDataAccLtFlt(1:nBinsLt,j+2)./histDataAccLtFlt(1,j+2);
% for cum fitting
ltValuesFlt=[];
for i=1:size(data,2)
    ltValuesFlt = vertcat(ltValuesFlt,data(i).accLtFlt);
end
ltValuesFlt=ltValuesFlt/dt;

%%

%% =======================================================================%
%                                                                         %
%                        Local Callback Functions                         %
%                                                                         %
%=========================================================================%


%% filter matrix array

function fltTraces=applyFilterToData(fretData,idxFilter)

    
    % Update data
    [nTraces,nFrames]=size(fretData.donor); 
    
    % Create traces object
    fltTraces = TracesFret(nTraces,nFrames);
    
    % Metadata Field Names
    metaFN = fieldnames(fretData.traceMetadata);
             
    % Update data
    fltTraces.donor            = fretData.donor(idxFilter,:);
    fltTraces.acceptor         = fretData.acceptor(idxFilter,:);
    fltTraces.fret             = fretData.fret(idxFilter,:);
    fltTraces.time             = fretData.time;
    
    fltTraces.traceMetadata    = [];

    % Update traceMetadata for each trace 
    j=0;
    for i=1:nTraces

        % Apply Filter
        if idxFilter(i)==0; continue;end
        j=j+1;

        %Update Metadata
        for k=1:numel(metaFN)
            fn=metaFN{k};
            if isfield(fretData.traceMetadata,fn)
                fltTraces.traceMetadata(j).(fn) =...
                fretData.traceMetadata(i).(fn);
            end
        end

    end % End Update Metadata
end
%% End
