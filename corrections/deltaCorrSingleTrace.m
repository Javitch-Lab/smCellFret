function [cIntAcc, cIntDon, dcorr, cFret] = deltaCorrSingleTrace(...
    delta, time, intAcc, intDon, accRange, totalAccRange, idl, meanTotal)
%--------------------------------------------------------------------------
%
% deltaCorrSingleTrace:
%   This is a sub-function called by the script exportFretTracesDif2Origin 
%
% Description:
%   The function applies the direct excitation (delta) correction to the 
%   acceptor data of a particular trace and recalculates the FRET efficiency.  
%
% Syntax:  
%   [cIntAcc, cIntDon, dcorr, cFret] = ...
%    deltaCorrSingleTrace(delta, time, intAcc, intDon, accRange,...
%    totalAccRange, idl, meanTotal)
% 
% Inputs:
%   1. delta - value of the direct excitation factor
%   2. time -  time range of the acceptor from track start to end 
%      (including backset and baseline data points, unit in sec). 
%   3. intAcc - acceptor intensities from the beginning to the end of the
%      acceptor track (including backset and baseline data points)
%   4. intDon - donor intensities from the beginning to the end of the  
%      donor track (including backset and baseline data points)
%   5. accRange - binary vector with the elements 0 (no tracking) and 1 
%      (tracking), which defines the track range of the acceptor (in frames) 
%   6. totalAccRange - range of the acceptor from track start to end 
%      (including backset and baseline data points, unit in frames).
%   7. idl -  binary vector with the elements 0 (no FRET) and 1 (FRET) from 
%      the beginning to the end of the donor track (including backset and
%      baseline data points)
%   8. meanTotal - nanmean of the total intensity over the track range of
%      the donor   
% 
% Outputs:
%   1. cIntAcc - corrected acceptor intensity 
%   2. cIntDon - corrected donor intensity
%   3. dcorr - the calculated acceptor intensity correction
%   4. cFret - corrected FRET value
% 
% See also: 
%   exportFretTracesDif2Origin.m
%
% Author: 
%   - P.G. Sep 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

% range
dimAcc = size(totalAccRange);
range = dimAcc(1):dimAcc(2);

% corrected donor
cIntDon = intDon;

% corrected acceptor
notZero = intAcc(accRange)>0; 
dcorr   = delta*meanTotal;
cIntAcc = intAcc;
cIntAcc(accRange) = (intAcc(accRange) - dcorr ).*notZero;

% corrected fret
cFret  = (cIntAcc./(intDon(range) + cIntAcc));
      
cFret = cFret.*idl(range);
      
% plot
figure;
subplot(2,1,1)
plot(time(range), cIntDon(range),'g');hold on; grid minor;
title('Direct Excitation Correction');

subplot(2,1,1)
plot(time(range), cIntAcc(range),'r')

subplot(2,1,2)
plot(time(range), cFret(range)); grid minor;

end

