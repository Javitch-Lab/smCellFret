function segments = getSegments(traceSet,trace_i)
%--------------------------------------------------------------------------
%
% getSegments.m:
%   Subfunction of scriptFretFilter.m     
% 
% Description:
%   The function finds all tracks in a track set that cross the 
%   neighborhood region of a given track i at any given time. The 
%   neighborhood region R is a square region defined by the outer 
%   boundaries of track i. 
%
% Syntax:  
%   segments = getSegments(trackSet,track_i)
% 
% Inputs:
%   1. trackSet - is a fretTraces structure, as is saved in fretTraces*.mat 
%      files with the following fields:
%      .time
%      .x
%      .y
%      .xCorr
%      .yCorr
%      .int
%      .snr
%      .traceMetadata
%   2. track_i - Integer number specifying a track_i in the trackset. The 
%      function evaluates which tracks in the trackset cross the square 
%      range R of track_i at any point in time (in the example below
%      track_i=10).
% 
% Outputs:
%   1. segments - Cell array with entries:
%   No. 	SOT(s)  EOT(s)     Centroid (pixels)	     Trace ID
%   10      0.040	4.200   [ 133.954	165.398	]	'singleTracesCh1_10'
%   9       0.040	4.560	[ 148.876	170.771	]	'singleTracesCh1_9'
%   95      4.280	7.120	[ 143.266	167.642	]	'singleTracesCh1_95'
%   99      4.440	22.320	[ 133.262	156.051	]	'singleTracesCh1_99'
%   153     7.120	9.640	[ 143.822	162.194	]	'singleTracesCh1_153'
%   856	    4.560	64.680	[ 133.890	168.041	]	'singleTracesCh1_856'
%   1207	104.720	104.760	[ 132.724	170.769	]	'singleTracesCh1_1207'
%   1320	116.800	116.880	[ 130.500	163.797	]	'singleTracesCh1_1320'
%   1524	139.560	139.600	[ 134.509	161.655	]	'singleTracesCh1_1524'
%
%   No :      track number 
%   SOT:      Start of track 
%   EOT:      End of track
%   Centroid: xy centroid coordinates of the track in pixels
%   Trace ID: Unique identifier of the track
%   
%   The list starts with track_i (No=10). All following tracks in the list 
%   are tracks that cross the neighborhood region R of track_i.          
% 
% See also: 
%   repaetingEventsFilter.m, scriptFretFilter.m
%
% Authors: 
%   - P.G. Sep 2018
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


% Initialize Variables
dt           = traceSet.time(2)-traceSet.time(1);
nTraces      = size(traceSet.x);
segments     = [];  
startOfTrace = traceSet.traceMetadata(trace_i).startOfTrace;
endOfTrace   = traceSet.traceMetadata(trace_i).endOfTrace;

% Calculate square region R defined by the boundary of the trace with
% identifier number trace_i. Use corrected coordinates if available.
if isfield(traceSet,'xCorr')
    r1 = [traceSet.xCorr(trace_i,startOfTrace:endOfTrace);... % r1(x1,x2,x3,...)
          traceSet.yCorr(trace_i,startOfTrace:endOfTrace)];   %   (y1,y2,y3,...)
else
    r1 = [traceSet.x(trace_i,startOfTrace:endOfTrace);...     
          traceSet.y(trace_i,startOfTrace:endOfTrace)];       
end

% The square region R
delta = 2; % unit in pixels
r1_Min  = [min(r1(1,:)-delta);min(r1(2,:))-delta]; % lower left corner of R 
                                                   % r1_Min(xMin) 
                                                   %       (yMin) 
                                                   
                                                   
r1_Max  =[max(r1(1,:)+delta);max(r1(2,:))+delta];  % upper right corner of R
                                                   % r1_Max(xMax)
                                                   %       (yMax) 

% Calculate Centroid of trace r1
r1Centroid = [(nanmean(r1(1,:))),(nanmean(r1(2,:)))];

% Generate ID        
traceID = {trace_i,...                            % Trace Number
           startOfTrace*dt,...                    % SOT
           endOfTrace*dt,...                      % EOT
           r1Centroid,...                         % Centroid Coord
           traceSet.traceMetadata(trace_i).ids};    % ID Name

% Loop over all traces and check if segments exist in square region R
for i=1:nTraces

    if i==trace_i, continue, end        
            
    startOfTrace = traceSet.traceMetadata(i).startOfTrace;
    endOfTrace   = traceSet.traceMetadata(i).endOfTrace;
    
    if isfield(traceSet,'xCorr')
        r2 = [traceSet.xCorr(i,startOfTrace:endOfTrace);... % r2(x1,x2,x3,...)
             traceSet.yCorr(i,startOfTrace:endOfTrace)];    %   (y1,y2,y3,...)
    else
        r2 = [traceSet.x(i,startOfTrace:endOfTrace);...     
             traceSet.y(i,startOfTrace:endOfTrace)];        
    end
    
    % Check if any of the r2-coord are in region R
    inRegion   = (r2>r1_Min & r2<r1_Max); 
    % r2(x1,x2,x3,...) > r1_Min(xMin)  AND  r2(x1,x2,x3,...) < r1_Max(xMax)
    %   (y1,y2,y3,...) >       (yMin)         (y1,y2,y3,...) <       (yMax) 
    
    % inRegion (  0	0 0 0 1 1 1 0 0 0 ...)
    %          (  0 0 0 0 1 1 1 0 0 0 ...)

    isInRegion = any(sum(inRegion,1)==2);

    if isInRegion     
        % Calculate Centroid of Trace
        r2Centroid = [(nanmean(r2(1,:))),(nanmean(r2(2,:)))];
        % Generate Ids         
        tmpSeg = {i, ...                            % Trace Number
                  startOfTrace*dt,...               % SOT
                  endOfTrace*dt,...                 % EOT
                  r2Centroid,...                    % Centroid Coord
                  traceSet.traceMetadata(i).ids};     % ID Name
        segments = [segments; tmpSeg]; %#ok<AGROW>
    else
        continue
    end
end

% Generate Output 
segments = [traceID; segments];
end
