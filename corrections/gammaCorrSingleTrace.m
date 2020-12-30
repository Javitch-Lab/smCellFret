function [cIntAcc, cIntDon, gcorr, cFret] = ...
    gammaCorrSingleTrace(gamma,time, intAcc, intDon, totalAccRange,idl)
%--------------------------------------------------------------------------
%
% gammaCorrSingleTrace:
%   This is a sub-function called by the script exportFretTracesDif2Origin 
%
% Description:
%   The function applies the brightness (gamma) correction to the 
%   donor data of a particular trace and recalculates the FRET efficiency.  
%
% Syntax:  
%   [cIntAcc, cIntDon, gcorr, cFret] = ...
%    gammaCorrSingleTrace(gamma, time, intAcc, intDon, totalAccRange, idl)
% 
% Inputs:
%   1. gamma - value of the brightness correction factor
%   2. time -  time range of the acceptor from track start to end 
%      (including backset and baseline data points, unit in sec). 
%   3. intAcc - acceptor intensities from the beginning to the end of the
%      acceptor track (including backset and baseline data points)
%   4. intDon - donor intensities from the beginning to the end of the  
%      donor track (including backset and baseline data points)
%   5. totalAccRange - range of the acceptor from track start to end 
%      (including backset and baseline data points, unit in frames).
%   6. idl -  binary vector with the elements 0 (no FRET) and 1 (FRET) from 
%      the beginning to the end of the donor track (including backset and
%      baseline data points)  
% 
% Outputs:
%   1. cIntAcc - corrected acceptor intensity 
%   2. cIntDon - corrected donor intensity
%   3. gcorr - the calculated donor intensity correction
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

% total acc range
dimAcc = size(totalAccRange);
range = dimAcc(1):dimAcc(2);

% corrected donor
gcorr = gamma.*intDon(range);
cIntDon = gcorr;

% corrected acceptor
cIntAcc = intAcc;

% corrected fret
cFret  = intAcc./(gcorr + intAcc);
cFret = cFret.*idl(range);
      
% plot 
figure
subplot(2,1,1)
plot(time(range), cIntDon(range),'g');hold on; grid minor
title('Gamma Correction');

subplot(2,1,1)
plot(time(range), cIntAcc(range),'r')

subplot(2,1,2)
plot(time(range), cFret(range)); grid minor

end

