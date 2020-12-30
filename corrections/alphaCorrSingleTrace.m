function [cIntAcc, cIntDon, acorr, cFret] = ...
    alphaCorrSingleTrace(alpha,time, intAcc, intDon, totalAccRange, idl)
%--------------------------------------------------------------------------
%
% alphaCorrSingleTrace:
%   This is a sub-function called by the script exportFretTracesDif2Origin 
%
% Description:
%   The function applies the crosstalk (alpha) correction to the acceptor  
%   data of a particular trace and recalculates the FRET efficiency.  
%
% Syntax:  
%   [cIntAcc, cIntDon, acorr, cFret] = ...
%    alphaCorrSingleTrace(alpha,time, intAcc, intDon, totalAccRange, idl)
% 
% Inputs:
%   1. alpha - value of the crosstalk factor
%   2. time -  time range of the acceptor from track start to end 
%      (including backset and baseline data points, unit in sec). 
%   3. intAcc - acceptor intensities from the beginning to the end of the
%      acceptor track (including backset and baseline data points)
%   4. intDon - donor intensities from the beginning to the end of the  
%      donor track (including backset and baseline data points)
%   5. totalAccRange - range of the acceptor from track start to end 
%      (including backset and baseline data points, unit in frames)
%   6. idl -  binary vector with the elements 0 (no FRET) and 1 (FRET) from 
%      the beginning to the end of the donor track (including backset and
%      baseline data points) 
% 
% Outputs:
%   1. cIntAcc - corrected acceptor intensity 
%   2. cIntDon - corrected donor intensity
%   3. acorr - the calculated acceptor intensity correction
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
acorr = alpha*intDon(range);
cIntAcc = intAcc - acorr;

% corrected fret
cFret  = (intAcc - acorr )./(intDon(range) + (intAcc - acorr));
cFret = cFret.*idl(range);
      
% plot 
figure
subplot(2,1,1)
plot(time(range), cIntDon(range),'g');hold on; grid minor
title('Crosstalk Correction');

subplot(2,1,1)
plot(time(range), cIntAcc(range),'r')

subplot(2,1,2)
plot(time(range), cFret(range)); grid minor

end

