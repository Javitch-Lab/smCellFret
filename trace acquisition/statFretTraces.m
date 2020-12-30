function  statFretTraces( accFile, donFile, flagExtBgd )
%--------------------------------------------------------------------------
%
% statFretTraces:
%   Updates the fretTraces traceMetadata field with additional statistical
%   quantities. 
%
% Description:
%   Subfunction of scriptGetFretTraces.m
% 
% Syntax:  
%   statFretTraces( accFile, donFile, flagExtBgd )
% 
% Inputs:
%   1. accFile    - Acceptor traces file with intensities in units of 
%                   photons/frame (files with extension 'Acc_phc.traces')
%   2. donFile    - Donor traces file with intensities in units of 
%                   photons/frame (files with extension 'Don_phc.traces')
%   3. flagExtBgd - binary variable with  
%                   FALSE = switch off extended background calculation 
%                   TRUE  = switch on extended background calculation
% 
% Outputs:
%   The following statistical measures are calculated for each trace:    
%     meanInt      - Mean value of local acceptor intensity (locIntAcc) 
%                    minus mean value of local acceptor background  
%                    (locBgdAcc). See also '.int' and '.snr' in 
%                    scriptGetFretTraces.m for the definition of locInt and 
%                    locBgd.
%     meanXIntAcc  - locIntAcc minus the mean value of the extended local 
%                    acceptor background (locXBgdAcc). The extended local  
%                    background is used for traces that have a length of  
%                    less than 10 frames (see function statFret.Traces.m)
%     meanXIntDon  - mean of locIntDon minus mean of locXBgdDon. 
%     meanXSnrAcc  - mean of locIntAcc diveded by mean of locXBgdAcc
%     meanXSnrDon  - mean of locIntDon diveded by mean of locXBgdDon
%     meanXCrt     - meanXIntAcc diveded by meanXIntDon
%      
%   These measures are added to the donor and acceptor traceMetadata field 
%   of the fretTraces structure. The updated fretTraces of the donor and 
%   acceptor structures are saved under the same filenames as the two input
%   files.  
% 
% Other m-files required: 
%   Subfunctions: extendedLocBgd.m
% 
% See also: 
%   scriptGetFretTraces.m
%
% Authors: 
%   - P.G. Sep 2018
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% Load files for test purposes
% clear;
% currentDir = 'F:\Disk E\#BackupMicData\2017 May-Jun\170511_WA-MH_SNAPf-mGlu2-Cy3AC-Cy5AC_w&wo Glu_Live FRET\333 nM Cy3AC_666 nM Cy5AC_PassX_100 ngmL tet_50 PCA PCD\Track_Opt28-TW5_52418\LT20FreeNoCrt\';
% dataSetCh1='#01Ch1\';
% accName='singleTracesAcc_phc.traces';
% dataSetCh2='#01Ch2\';
% donName='singleTracesDon_phc.traces';
% accFile=[currentDir dataSetCh1 accName];
% donFile=[currentDir dataSetCh2 donName];
% flagExtBgd=1;

%% Load files

rawTraces.Acc=loadTracesCell(accFile);
rawTraces.Don=loadTracesCell(donFile);
try
    assert(size(rawTraces.Acc.x,1)==size(rawTraces.Don.x,1))
catch
    warning('Acceptor and Donor file have different sizes');
end

%% Initialize Variables

[nTraces,nFrames] = size(rawTraces.Acc.x);
pulseLength = 10;

%% Calculate Trace Statistics
           
h=waitbar(0,'calc Trace Statistics ...');

for i=1:nTraces

    %--- Load Acceptor Data
    startOfTraceAcc = rawTraces.Acc.traceMetadata(i).startOfTrace;
    endOfTraceAcc   = rawTraces.Acc.traceMetadata(i).endOfTrace;
    traceLenAcc     = rawTraces.Acc.traceMetadata(i).traceLen;
    meanLocIntAcc   = nanmean(rawTraces.Acc.locInt(i,...
                              startOfTraceAcc:endOfTraceAcc));
    meanLocBgdAcc   = nanmean(rawTraces.Acc.locBgd(i,...
                              startOfTraceAcc:endOfTraceAcc));
    rawTraces.Acc.traceMetadata(i).meanInt = meanLocIntAcc - meanLocBgdAcc;
    
    % Calculate mean background and mean standard deviation
    if traceLenAcc < pulseLength && flagExtBgd==1
        % Calculate the meanLocBgd and stdvLocBgd over an extended region
        % for short traces
        meanLocBgdAcc = extendedLocBgd(rawTraces.Acc,i,nFrames,pulseLength);
    end
  
    %--- Load Donor Data 
    startOfTraceDon = rawTraces.Don.traceMetadata(i).startOfTrace;
    endOfTraceDon   = rawTraces.Don.traceMetadata(i).endOfTrace;
    traceLenDon     = rawTraces.Don.traceMetadata(i).traceLen;
    meanLocIntDon   = nanmean(rawTraces.Don.locInt(i,...
                              startOfTraceDon:endOfTraceDon)); 
    meanLocBgdDon   = nanmean(rawTraces.Don.locBgd(i,...
                              startOfTraceDon:endOfTraceDon));
%   This is the mean int. calc bf joining donor tracks (don and acc
%   tracks have the same length)
%   rawTraces.Don.traceMetadata(i).meanInt = meanLocIntDon - meanLocBgdDon;
    
    % Calculate mean background and mean standard deviation
    if traceLenDon < pulseLength && flagExtBgd==1
        % Calculate the meanLocBgd and stdvLocBgd over an extended region
        % for short traces
        meanLocBgdDon = extendedLocBgd(rawTraces.Don,i,nFrames,pulseLength);
    end
    
    %--- Save additional data in FretTraces Metadata field
    meanXIntAcc = meanLocIntAcc - meanLocBgdAcc;
    rawTraces.Acc.traceMetadata(i).meanXIntAcc = meanXIntAcc;
    
    %This is the extended int. calc bf joining donor tracks (don and acc
    %tracks have the same length)
    meanXIntDon = meanLocIntDon - meanLocBgdDon;
    rawTraces.Acc.traceMetadata(i).meanXIntDon = meanXIntDon; 
    
    meanXSnrAcc = (meanLocIntAcc/meanLocBgdAcc);
    rawTraces.Acc.traceMetadata(i).meanXSnrAcc = meanXSnrAcc;
    
    %This is the extended Snr. calc bf joining donor tracks (don and acc
    %tracks have the same length)
    meanXSnrDon = (meanLocIntDon/meanLocBgdDon);
    rawTraces.Acc.traceMetadata(i).meanXSnrDon = meanXSnrDon; 
    
    meanXCrt = (meanLocIntAcc - meanLocBgdAcc)/...
               (meanLocIntDon - meanLocBgdDon);
    rawTraces.Acc.traceMetadata(i).meanXCrt = meanXCrt;   
    
    %--- Update Waitbar
    waitbar(i / nTraces)
    
end
close(h);

%% Save to File
disp('Save traces ...')
saveTracesCell(accFile,'traces',rawTraces.Acc);
saveTracesCell(donFile,'traces',rawTraces.Don);


%% End





