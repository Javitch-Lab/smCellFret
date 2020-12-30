function [ meanExtLocBgd, stdvExtLocBgd,result] = ...
         extendedLocBgd( S,traceNo,nFrames,minTraceLen )
%--------------------------------------------------------------------------
%
% extendedLocBgd.m:
%   Calculates the mean local spot background intensity for short time 
%   traces.
% 
% Description:
%   This is a subfunction of statFretTraces.m
% 
% Syntax:  
%   [ meanExtLocBgd, stdvExtLocBgd, result] = extendedLocBgd(S,traceNo,...
%                                            nFrames,minTraceLen )
% 
% Inputs:
%   1. S           - Structure variable that must contain the fields of a 
%                    trace structure as used in a *_phc.traces file format. 
%   2. traceNo     - Integer number identifying the specific trace being 
%                    stored in the structure S.  
%   3. nFrames     - Integer number that specifies the number of movie
%                    frames that are saved in the structure S.   
%   4. minTraceLen - If the length of a trace (in frames) is below this 
%                    limit, an extended background range is used to 
%                    calculate the value of the mean local background.
%
% Outputs:
%   1. meanExtLocBgd - Value of the mean extended local background 
%                      intensity. 
%   2. stdvExtLocBgd - Value of the standard deviation of the extended mean 
%                      local background intensity. 
%   3. result        - Not used (for test purposes only)
% 
% See also: 
%   statFretTraces.m, scriptGetFretTraces.m
%
% Authors: 
%   - P.G. May 2017
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Initialize Parameters
i            = traceNo;
lenBackset   = S.traceMetadata(i).lenBackset;
lenBaseLine  = S.traceMetadata(i).lenBaseline;
startOfTrace = S.traceMetadata(i).startOfTrace;
endOfTrace   = S.traceMetadata(i).endOfTrace;
traceLen     = S.traceMetadata(i).traceLen;   
leftLim      = lenBackset;
rightLim     = lenBaseLine;

%% Calculate mean Bgd / stdv over the extended region 
try
    % The the length of the actual trace has to be smaller than the minTraceLen (pulse length, e.g 10 frames) 
    assert(traceLen < minTraceLen)
    
    % Find Background Limits
    if leftLim < round(0.5*minTraceLen)
        % Pulse starts near the left border  
        winLeft=0;
        while startOfTrace - winLeft >1 && winLeft <(minTraceLen-traceLen)
            winLeft=winLeft+1;
        end

        winRight=0;
        while endOfTrace + winRight < nFrames && winRight < (minTraceLen-traceLen-winLeft)
            winRight=winRight+1;
        end

    elseif rightLim < round(0.5*minTraceLen)
        % Pulse starts near the right border 
        winRight=0;
        while endOfTrace + winRight < nFrames && winRight < (minTraceLen-traceLen)
            winRight=winRight+1;
        end

        winLeft=0;
        while startOfTrace - winLeft >1 && winLeft <(minTraceLen-traceLen-winRight)
            winLeft=winLeft+1;
        end

    else
        % Pulse is distant from the left and right border
        winLeft = round(0.5*(minTraceLen-traceLen));
        winRight = minTraceLen - traceLen - winLeft;

    end

    % Calculate local background over the extended region
    meanExtLocBgd = nanmean(S.locBgd(i,startOfTrace-winLeft:endOfTrace+winRight));
    stdvExtLocBgd = nanstd(S.locBgd(i,startOfTrace-winLeft:endOfTrace+winRight));
    
    % Output: For Test Purposes
    % result = [leftLim startOfTrace endOfTrace rightLim traceLen winLeft winRight];
    % disp(result);
    
catch
    warning('trace length is >= min trace length')
    meanExtLocBgd    = nanmean(S.locBgd(i,startOfTrace:endOfTrace));
    stdvExtLocBgd    = nanstd(S.locBgd(i,startOfTrace:endOfTrace));
end

       
%% END
 
