function baseline = makeBaseline(filename, lengthBackset,...
                                 lengthBaseline, stdBgd) 
%--------------------------------------------------------------------------
%
% makeBaseline.m:
%   Generates baseline coordinates for each track before the real track 
%   starts and after the real track ends.
% 
% Description:
%   The function adds a defined number of (x,y) coordinates (backset) to 
%   the real track coordinates which are temporally before the start point 
%   of the track. It also adds a defined number of (x,y) coordinates 
%   (baseline) to the real track coordinates that are time-wise behind the 
%   end point of the track. The intensity values measured at the backset 
%   and baseline positions reflect the background intensities of the 
%   particular track.
% 
% Syntax:  
%   baseline = makeBaseline(filename, lengthBackset,...
%                                 lengthBaseline, stdBgd) 
% 
% Inputs:
%   1. filename       - files of type 'singleTracesAcc.traces'  
%   2. lengthBackset  - Integer number specifying the length of the backset 
%                       in frames. 
%   3. lengthBaseline - Integer number specifying the length of the
%                       baseline in frames. 
%   4. stdBgd         - No longer in use
% 
% Outputs:
%   The function returns the variable baseline and saves two files in the
%   folder #Ch2:
%   1. baseline - workspace variable with 4 columns and N (N = number of 
%      traces) rows:
%       a) lenBackset (1st column: length of the backset).
%       b) startOfTrace (2nd column: start of the current track)  
%       c) endOfTrace (3rd column: end of current track)
%       d) lenBaseline (4th column: length of the baseline)     
%   2. baseline.txt - ascii-text file with 4 columns and N (N = number of 
%      traces) rows:
%       a) lenBackset (1st column: length of the backset).
%       b) startOfTrace (2nd column: start of the current track)  
%       c) endOfTrace (3rd column: end of current track)
%       d) lenBaseline (4th column: length of the baseline)
%   3. 'singleTracesAcc_bgd.traces' - traces file containing the updated 
%       baseline information. 
% 
% See also: 
%   scriptGetFretTraces.m
%
% Authors: 
%    -P.G. Aug 2011
%    -P.G. Aug 2014
%    -P.G. Sep 2014
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Load the traces files
data  = loadTracesCell(filename);
xCor  = data.x;
yCor  = data.y;

%Create an output filename
[path,name,~]=fileparts(filename);
outfile = [path filesep name '.traces'];
outfile = strrep(outfile,'.traces','_bln.traces');

%% Calculate length of Baseline and Backset for each trace and save in 
%  variable baseline
[nTraces,nFrames]=size(xCor);
baseline=zeros(nTraces,4);
trackLogical=abs((xCor==0)-1); % Tracks have logical 1 when tracked else 0.

for i=1:nTraces
    % Calculate startOfTrace and lenBackset
    tmpStartOfTrace = find(trackLogical(i,:)==1,1,'first');
    if tmpStartOfTrace<=lengthBackset
        tmpLenBackset=tmpStartOfTrace-1;
    else
        tmpLenBackset=lengthBackset;
    end
    baseline(i,1)= tmpLenBackset;    % save lenBackset of track i 
    baseline(i,2)= tmpStartOfTrace ; % save startOfTrace of track i
    
    % Calculate endOfTrace and lenBaseline
    tmpEndOfTrace=find(trackLogical(i,:)==1,1,'last');
    maxLen = nFrames-tmpEndOfTrace;
    if maxLen<lengthBaseline
       tmpLenBaseline = maxLen;
    else
       tmpLenBaseline=lengthBaseline;
    end 
    baseline(i,3)= tmpEndOfTrace;  % save endOfTrace of track i
    baseline(i,4)= tmpLenBaseline; % save lenBaseline of track i
    
end

%% Calculate trace length and append baseline
for n=1:nTraces
    lenBackset   = baseline(n,1);
    startOfTrace = baseline(n,2); 
    endOfTrace   = baseline(n,3);
    lenBaseline  = baseline(n,4);
    
    % append backset to track 
    if startOfTrace > 1
        randVectorX = xCor(n,startOfTrace); %+ stdBgd.*randn(lenBackset,1);
        xCor(n,startOfTrace-lenBackset:startOfTrace-1)= randVectorX';
        randVectorY = yCor(n,startOfTrace); %+ stdBgd.*randn(lenBackset,1);
        yCor(n,startOfTrace-lenBackset:startOfTrace-1)= randVectorY';
    end
    
    % append baseline to track
    if endOfTrace < nFrames
        randVectorX = xCor(n,endOfTrace); %+ stdBgd.*randn(lenBaseline,1);
        xCor(n,endOfTrace+1:endOfTrace+lenBaseline)= randVectorX';
        randVectorY = yCor(n,endOfTrace); %+ stdBgd.*randn(lenBaseline,1);
        yCor(n,endOfTrace+1:endOfTrace+lenBaseline)= randVectorY';
    end

end

%% Save the modified traces to file 

data.x      = xCor;
data.y      = yCor;
saveTracesCell(outfile,'traces',data);
dlmwrite([path '\baseline.txt'],baseline);

