function [ xGapsClose , seg ] = interpGaps( xGaps )
%--------------------------------------------------------------------------
%
% interpGaps.m:
%   The function interpolates zero entries (gaps) in 1dim intensity time 
%   trace data. Leading zero entries and trailing zero entries in the time
%   trace are not considered. 
%
% Description:
%   The function:  
%   a) finds the start and end of a trace with gaps (zeros): 
%      e.g. xGaps = [ 0 0 0 1 4 3 6 7 0  8 5 4 0 0 0 10 12 0 0 0 0 0]
%      start at idx 4
%      end   at idx 17
%   b) locates segments (zero entries) within the trace:
%      zeros at seg [9 9] [13 15]
%   c) interpolates these segments
%
% Syntax:  
%   [ xGapsClose , seg ] = interpGaps( xGaps )
% 
% Inputs:
%   xGaps - Row vector with intensity values from a time trace.
% 
% Outputs:
%   1. xGapsClose - the interpolated intensity time trace, e.g.:
%      xGapsClose:   [ 0  0  0  1  4  3  6  7  7.5  8  5  4  5.5  7  8.5  10  12  0  0  0  0  0]
%      idx             1  2  3  4  5  6  7  8  9   10 11 12  13  14  15   16  17  18 19 20 21 22 
%   2. seg - Location of the segments: e.g. 
%      seg = 1     3
%            9     9
%            13    15
%            18    21
% 
% Other m-files required: 
%   Subfunctions: uses matlab function 'interp1.m'
%
% Authors: 
%   - P.G. Jul 2013
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% localize gaps in data sequence 
b=xGaps==0;
i=0;j=0;s=0;
nFrames=size(xGaps,2);
while i<nFrames
    i=i+1;
    while b(i+j)>0 && i+j<=nFrames
        j=j+1;
        if i+j>nFrames;break;end
    end
    if j>0
        s=s+1;
        seg(s,1)=i;
        seg(s,2)=i+j-1;
        i=i+j;
        j=0;
    end
end

%% Interpolate segments

% Do not interpolate zero entries at the beginning and the end of the trace
nSeg=size(seg,1);
if seg(1,1) == 1 && nSeg>1; seg=seg(2:end,:);nSeg=nSeg-1; end
if seg(nSeg,2) == nFrames && nSeg>1; seg=seg(1:nSeg-1,:);nSeg=nSeg-1; end

% Interpolate remaining segments
if nSeg>1
    for k=1:nSeg
        t  = [seg(k,1)-1 seg(k,2)+1];
        ti = seg(k,1) : seg(k,2);
        xGapInterval  = [xGaps(seg(k,1)-1) xGaps(seg(k,2)+1)];
        xi = interp1(t,xGapInterval,ti);
        xGaps(seg(k,1) : seg(k,2))= xi;
    end
    xGapsClose=xGaps;
else
    xGapsClose=xGaps;
end

