function [trackIndex,data] = trackMate2traces( filename, dt, nFrames, roiData )
%--------------------------------------------------------------------------
%
% trackMate2traces:
%   Subfunction called by imageRegistrationTrackMate.m that imports 
%   a grid calibration file. 
% 
% Description:
%   The function imports tracking data measured with TrackMate [1] (file 
%   Spots_in_tracks_statistics.txt) and converts it into a trace file 
%   format.       
% 
% Syntax:  
%   [trackIndex,data] = trackMate2traces( filename, dt, nFrames, roiData )
% 
% Inputs:
%   1. filename - file name of the calibration file 
%      'Spots in tracks statistics.txt'
%   2. dt - time resolution in s 
%   3. nFrames - Number of frames of the calibration movie 
%      (nFrames = rows*cols) 
%   4. roiData - Polygon area described by x-coordinates (1st column) and y 
%      coordinates (2nd column). Only tracks within the ROI are analyzed. 
% 
% Outputs:
%   1. trackIndex - Structure with track IDs (this output is obsolete and 
%      no longer in use)
%   2. data - traces file structure
% 
% Other m-files required: 
%   Subfunctions: txtReadTrackMate.m
% 
% See also: 
%  imageRegistrationTrackMate.m
%
% References:
%   [1] Tinevez, J. Y., et al. (2017). "TrackMate: An open and extensible 
%   platform for single-particle tracking." Methods 115: 80-90.
%
% Authors: 
%   - P.G. Feb 2013
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% Import data from spreadsheet
disp('importing trackMate data...');

path=fileparts(filename);
tracesNum = txtReadTrackMate(filename);
% find and delete rows with zero entries
%idx=find(tracesNum(:,10)==0);
%for i=idx
    %delete empty row
    %tracesNum(i,:)    = [];
%end

%% Determine nTraces in data set
tmpInd=0;nTraces=0;
maxInd=size(tracesNum,1);
while tmpInd<maxInd
    tmpInd=find(tracesNum(:,4)==tracesNum(tmpInd+1,4),1,'last'); 
    nTraces=nTraces+1;
end
disp(['singleTraces :' num2str(nTraces)]);

% Determine traces index: start and end
segStart=1;traceSeg=zeros(nTraces,2);
for i=1:nTraces
    segEnd=find(tracesNum(:,4)==tracesNum(segStart,4),1,'last'); 
    traceSeg(i,1)=segStart;
    traceSeg(i,2)=segEnd;
    segStart=segEnd+1;
end

%% Extract traces
% Translate coordinates
xTransl=1;
yTransl=1;
% Allocate memory  
x=zeros(nTraces,nFrames);
y=zeros(nTraces,nFrames);
int=zeros(nTraces,nFrames);
T=dt*(1:nFrames)';

for i=1:nTraces
    [frames,idx]  = sort(tracesNum(traceSeg(i,1):traceSeg(i,2),10));
    frames=frames+1;
    xUnsorted     = tracesNum(traceSeg(i,1):traceSeg(i,2),6);
    x(i,frames)   = xUnsorted(idx)+xTransl;
    yUnsorted     = tracesNum(traceSeg(i,1):traceSeg(i,2),7);
    y(i,frames)   = yUnsorted(idx)+yTransl;
    intUnsorted   = tracesNum(traceSeg(i,1):traceSeg(i,2),16);
    int(i,frames) = intUnsorted(idx);
end

%% Only analyze the subset of traces located in the ROI
if ~isempty(roiData)
    tracesStart=zeros(nTraces,2);
    for k=1:nTraces 
        ind=find(x(k,:)>0,1,'first');
        tracesStart(k,1)=x(k,ind);
        tracesStart(k,2)=y(k,ind);
    end
    xCor=tracesStart(:,1);
    yCor=tracesStart(:,2);
    xRoi=roiData(:,1);
    yRoi=roiData(:,2);
    in=inpolygon(xCor, yCor, xRoi, yRoi);
    %set(gcf, 'CurrentAxes',handles.axContourSel);
    %cla;
    %plot(xRoi,yRoi,xCor(in),yCor(in),'r+',xCor(~in),yCor(~in),'bo')
    %xlim(gca,[0 256]);
    %ylim(gca,[0 256]);
    x=x(in,:);
    y=y(in,:);
    int=int(in,:);
    nTraces=size(x,1);
    disp(['singleTraces (ROI) :' num2str(nTraces)]);
end

%% Interpolate traces 
for k=1:nTraces
   xi=interpGaps(x(k,:));
   yi=interpGaps(y(k,:));
   x(k,:)=xi;
   y(k,:)=yi;
end

%% Save Traces in a structure
data.time = T;
data.x    = x;
data.y    = y;
data.int  = int;

% create ids
name='singleTracesCh1';
idList=(1:nTraces)';
for j=1:nTraces
    data.traceMetadata(j).ids = sprintf('%s_%d', name, idList(j));
end

%save data
outfile=[path '\' name '.traces'];
saveTracesCell(outfile,'traces',data);

%% Generate variable: trackIndex 
trackIndex.single.ids=idList;
trackIndex.single.nChannels=1;

%%
