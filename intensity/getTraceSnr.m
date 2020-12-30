%--------------------------------------------------------------------------
%
% getTraceSnr:
%   Retrieve signal-to-noise (SNR) values from fretTraces*Seg*.mat files.
% 
% Description:
%   Compile SNR data of fretTraces*Seg*.mat files into histograms. Use
%   scriptFretTraces2Seg.m to generate time *Seg*.mat traces with one
%   diffusion mode only. 
%   NOTE: Time traces, can be segmented after reading out only one specific
%   diffusion mode (here the mode of free diffusion). This means that the 
%   value of their amplitude (intensity or SNR ) is set to zero at the 
%   points where the diffusion mode is either confined, immobile or
%   directional.  
%     Calculation of the mean SNR: For segmented time traces with one
%     diffusion mode, the SNR is only calculated if the signal amplitude is 
%     greater than zero.
%   
% Syntax:  
%   getTraceSnr
% 
% Inputs:
%   The script prompts the user to select one or multiple 
%   fretTraces*Seg*.mat files e.g. fretTracesDif_all_SegFre01.mat
% 
% Outputs:
%   1. histSnrAcc - Acceptor SNR histogram array of N x data sets. The last 
%      three columns contain a) the snr distribution, averaged over all 
%      data sets, b) the SNR standard deviations and c) the SNR cumulative 
%      distribution.   
%   2. histSnrDon - Donor SNR histogram array of N x data sets. The last 
%      three columns contain a) the snr distribution, averaged over all 
%      data sets, b) the SNR standard deviations and c) the SNR cumulative 
%      distribution.
%   3. The following figures will be generated: 
%       a. Array of N acceptor SNR histograms
%       b. Array of N donor SNR histograms
% 
% Other m-files required: 
%   Subfunctions: none
%   MAT-files required: none
% 
% See also: 
%   scriptFretTraces2Seg.m, getTraceInt.m
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
tfiles = getFiles('*Seg*.mat','Select *.mat files');

% Check if files array is empty
if any(cellfun('isempty',tfiles))== true; return;end

% Show selection
for z=1:numel(tfiles)
    [path,name,ext]=fileparts(tfiles{z});
    disp(path);
    disp([name ext]);
end

%% Get Snr

nCells = numel(tfiles);
data    =struct([]);

% Initialize Intensity Histogram
edges = 1:0.05:4;
nBins=numel(edges)-1;
histSnrDon=zeros(1,nCells+1);
histSnrAcc=zeros(1,nCells+1);


%% Analyze Data 

% Initialize figures and waitbar
figure(1); %snr acc
figure(2); %snr don
f = waitbar(0,'Load Data ...');

% Loop over selected experiment folders
for j=1:nCells
    clear tmp;
    waitbar(j / nCells,f)
    
    % load traces file
    tracesFile = tfiles{j};
    load(tracesFile);
    snrAcc = fretTraces.Ch1.snr;
    snrDon = fretTraces.Ch2.snr;
    [nTraces,nFrames]=size(snrAcc);
    
    % save filename
    data(j).file = tracesFile;
    
    %---
    %--- GET SNR
    %---
    
    % logical Idl Array 
    idlAcc = false(nTraces,nFrames);
    accSOT= [fretTraces.Ch1.traceMetadata.startOfTrace].'; %unit is frames
    accEOT= [fretTraces.Ch1.traceMetadata.endOfTrace].'; %unit is frames
    
    idlDon = false(nTraces,nFrames);
    donSOT= [fretTraces.Ch2.traceMetadata.startOfTrace].'; %unit is frames
    donEOT= [fretTraces.Ch2.traceMetadata.endOfTrace].'; %unit is frames
    
    for i=1:nTraces
        idlAcc(i,accSOT(i):accEOT(i))=true;
        idlDon(i,donSOT(i):donEOT(i))=true;
    end
    
    % Calc SNR per Trace over the region where trace is idealized / tracked 
    snrAcc_idl2D = snrAcc.*idlAcc;
    snrAcc_idl2D(snrAcc_idl2D==false)=NaN;
    meanSnrAccPerTrace = nanmean(snrAcc_idl2D,2);
    data(j).snrAcc = meanSnrAccPerTrace;
    
    snrDon_idl2D = snrDon.*idlDon;
    snrDon_idl2D(snrDon_idl2D==false)=NaN;
    meanSnrDonPerTrace = nanmean(snrDon_idl2D,2);
    data(j).snrDon = meanSnrDonPerTrace;
    
    % Plot SNR Histogramms
    % Acc
    figure(1)
    subplot(2,round(nCells/2),j)
    XsnrAcc=data(j).snrAcc;
    histAcc(j) = histogram(XsnrAcc,edges); 
    % Offset the bin width by half the bin width to plot data at bin
    % centers
    histSnrAcc(1:nBins,1) = (histAcc(j).BinEdges(1:nBins)+histAcc(j).BinWidth/2)';
    histSnrAcc(1:nBins,j+1) = histAcc(j).Values';
    title('Acceptor SNR Distribution');
    
    % Don
    figure(2)
    subplot(2,round(nCells/2),j)
    XsnrDon=data(j).snrDon;
    histDon(j) = histogram(XsnrDon,edges); 
    % Offset the bin width by half the bin width to plot data at bin
    % centers
    histSnrDon(1:nBins,1) = (histDon(j).BinEdges(1:nBins)+histDon(j).BinWidth/2)';
    histSnrDon(1:nBins,j+1) = histDon(j).Values';
    title('Donor SNR Distribution');
    
    % Number of Molecules
    data(j).nTraces = size(snrAcc,1);
    
end
close(f)

%% Disp Number of Molecules
nAllTraces = sum([data.nTraces].');
disp(['nAllTraces = ' num2str(nAllTraces)]);

%%  Stat Acc
% Snr mean
histSnrAcc(1:nBins,j+2) = nanmean(histSnrAcc(:,2:nCells+1),2);
% Snr std
histSnrAcc(1:nBins,j+3) = nanstd(histSnrAcc(:,2:nCells+1),0,2);
% total cummulative
histSnrAcc(1:nBins,j+4) = sum(histSnrAcc(:,2:nCells+1),2);

% Snr mean
histSnrDon(1:nBins,j+2) = nanmean(histSnrDon(:,2:nCells+1),2);
% Snr std
histSnrDon(1:nBins,j+3) = nanstd(histSnrDon(:,2:nCells+1),0,2);
% total cummulative
histSnrDon(1:nBins,j+4) = sum(histSnrDon(:,2:nCells+1),2);
