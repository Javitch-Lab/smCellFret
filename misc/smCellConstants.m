function constants = smCellConstants(currentDir)
%--------------------------------------------------------------------------
%
% snCellConstants.m:
%   The function stores parameters needed by functions and scripts used
%   within the smCellFret pipeline.
% 
% Syntax:  
%   constants = smCellConstants(currentDir)
% 
% Inputs:
%   currentDir - This variable is no longer used. 
% 
% Outputs:
%   constants - A structure variable that contains the parameters of the 
%               smCellFret pipeline. 
% 
% See also: 
%   scriptGetFretTraces.m
%
% Authors: 
%   - P.G. 2012 - 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%
%---
%--- Location of Data 
%---
%constants.currentDir=currentDir;

%---
%--- General parametrs
%---

constants.dt        = 0.040;         % Time resolution.
 
%---
%--- Camera Settings
constants.QE = 0.97;                 % Quantum Efficiency
constants.photonConvFactorDelta = 1; % 100/1.3; % 20MHz, Gain 2, 
constants.pixelSize    = 160;        % 160 nm/pixel 

%---
%--- Parameters for particle detection (scriptDetectGeneral)
%--- 
constants.psfSigmaCh1  = 1.05;       % 100xNa1.49: 1.00994 (1x1)
constants.psfSigmaCh2  = 1.05;                                  
constants.intgrWinCh1  = 1;          % (0 for beads | 3 for mol.) 
constants.intgrWinCh2  = 2;          % nFrames before and after a frame for  
                                     % time avgeraging 
                                                                 
constants.testAlphaCh1   = ...
struct('alphaR',0.2,'alphaA',0.2,...
      'alphaD',0.2,'alphaF',0);
constants.testAlphaCh2   = ...
struct('alphaR',0.2,'alphaA',0.2,...
      'alphaD',0.2,'alphaF',0);
                                
% .alphaR: For the residuals test, comparing N+1-kernel fit to
%  N-kernal fit. Optional. Default: 0.05.
% .alphaA: For amplitude test. Optional. Default: 0.05.
% .alphaD: For distance test. Optional. Default: 0.05.
% .alphaF: Final residuals test, comparing residuals from final
%  fit to estimated background noise.
%  Optional. Default: 0.05.
% constants.alphaLocMaxCh1 =  0.0350;  % Alpha-value for initial detection of
% constants.alphaLocMaxCh2 =  0.0250; 

constants.alphaLocMaxCh1 =  0.0250;  % Alpha-value for initial detection of
constants.alphaLocMaxCh2 =  0.0250;  % local maxima 0.025 (0.000001)

constants.nKernelsCh1 = 1;           % 0: Fit only 1 kernel per local max.
constants.nKernelsCh2 = 1;           % 0: Fit only 1 kernel per local max.
                                     % 1: Fit > 1 kernel per local max. 
%---
%--- Parameters for Frame-To-Frame Linking (scriptTrackGeneral)
%--- 
constants.linMotionCh1  = 0;       % 
constants.linMotionCh2  = 0;
                                     
constants.minSrRadCh1 = 0.25;      % 
constants.minSrRadCh2 = 0.25;

constants.maxSrRadCh1 = 3;         % 
constants.maxSrRadCh2 = 2;          
                                     
constants.srRadFactCh1  = 3;       % 
constants.srRadFactCh2  = 2;        

constants.locDenFr2FrCh1  = 1;     %  
constants.locDenFr2FrCh2  = 0;                                       

%---
%--- GapClosing (scriptTrackGeneral)
%---
constants.gcTimeWindowCh1      = 10;   
constants.gcTimeWindowCh2      = 20; 

constants.gcMergeSplitCh1      = 0;
constants.gcMergeSplitCh2      = 0;

constants.gcMinTrkLenCh1     = 3;
constants.gcMinTrkLenCh2     = 2;

constants.gcLinMotionCh1    = 0;
constants.gcLinMotionCh2    = 0;

constants.gcMinSrRadCh1 = 2;
constants.gcMinSrRadCh2 = 2;

constants.gcMaxSrRadCh1 = 3;
constants.gcMaxSrRadCh2 = 2;

constants.gcBrownStdMultCh1    = 4;
constants.gcBrownStdMultCh2    = 3;

constants.gcBrownScalingCh1    = 0 ;
constants.gcBrownScalingCh2    = 0 ;

constants.gcLenForClassifyCh1  = 5;
constants.gcLenForClassifyCh2  = 5;

constants.gcUseLocalDensityCh1 = 1;
constants.gcUseLocalDensityCh2 = 0;

constants.gcLinStdMultCh1      = 1;
constants.gcLinStdMultCh2      = 1;

constants.gcLinScalingCh1      = 0.25;
constants.gcLinScalingCh2      = 0.25;

constants.gcMaxAngleVVCh1      = 45;
constants.gcMaxAngleVVCh2      = 45;

constants.gcGapPenaltyCh1      = 1.5;
constants.gcGapPenaltyCh2      = 1.5;

                                  
%---
%--- Parameters for 'utracks2traces'
%---
constants.amplFactor=1;              % raises the amplitude of the 
                                     % intensity time traces
                                     
%---
%--- Parameters for 'makeBaseline.m'
%---
constants.lenBackset = 25;           % 25 Num. of points tracked randm before
                                     % the actual particle is located for
                                     % the first time
constants.lenBaseline = 1500;        % 1500 max length of baseline. Append 
                                     % n-times random values of the xyCord. 
                                     % around the last tracking position to 
                                     % the data. 
constants.stdBgd      = 0.15;        % Standard deviation of the normal  
                                     % distribution which is used to  
                                     % randomly select values around the
                                     % last detected xyCord.
%---
%--- Parameters for 'nPhotons.m'
%---                              
constants.widthPSF   = 7;            % The psf area in pixels can be for
                                     % a) Binning 1x1, 100x Obj(1.49):7x7 
                                     % b) Binning 1x1, 100x Obj(1.49):5x5
                                     % c) Binning 2x2, 100x Obj(1.49):3x3
                                     % to capture >80% of the ligth choose
                                     % 7x7 (5x5 intPsf, sourounding area
                                     % used for bgd calculation (see
                                     % function sumPixelArea.m, nPhotons.m)

%---
%--- Parameters for 'filterTraces.m' 
%---
constants.ltLimit  = 1;              % Unit in frames. Acceptor traces with 
                                     % lifetimes >= ltLimit
                                     % will be accepted in cellFRETGettraces. 
                                     % If ltFlt=0, the variable will be set 
                                     % to the default value 1.
                                     
constants.nSigma = 0.5;              % number of local background standard 
                                     % deviations.
                                     % (meanInt-meanBgd)/stdvBgd >= nSigma

constants.maxSnrLimit = 0.5;         % Acceptor traces with a
                                     % maxSnr >= maxSnrLimit
                                     % will be accepted in cellFRETGettraces
                                     
constants.minMeanXSnrAcc = 1.4;      % Acceptor traces with 
constants.maxMeanXSnrAcc = 2.2;      % meanXSnrAcc > minMeanXSnrAcc and
                                     % meanXSnrAcc < maxMeanXSnrAcc 
                                     % will be accepted in cellFRETGettraces
                                     
constants.minCrtValue = 0.3;         % Acceptor traces with
constants.maxCrtValue = 10;          % meanLocAccInt/meanLocDonInt > minCrtValue
                                     % and
                                     % meanLocAccInt/meanLocDonInt < maxCrtValue
                                     % will be accepted in cellFRETGettraces
                                     
                                     
constants.accLwLim = 118;            % Acceptor traces with
                                     % meanLocInt - meanLocBgd >= accLwLim
                                     % will be accepted in cellFRETGettraces
                                     
constants.accUpLim = 50000;          % Acceptor traces with 
                                     % meanLocInt - meanLocBgd <= accUpLim
                                     % will be accepted in cellFRETGettraces
                                     
                                     % Multiple Events Filter Settings:
                                     % A series of segments (sgm) is filtered out if: 
constants.trackDist  = 1.5;          % The centrDist(sgm_i|sgm 1...N) < trackDist (in pixels) AND
constants.nEvents    = 25;           % The sgm_i belongs to series of sgm 1...N  > nEvents AND
constants.offTimeLim = 0.5;          % MeanOffTime of the sgms 1...N is  < offTimeLim
                                     
%---
%--- Parameters for 'tracehist_norm.m'
%---
constants.conLevel   = 0:0.0002:0.001; % normalized Data
constants.bin        = 20;
constants.hstBackset = 0;%17;          % Baseline before post-synch.
constants.fretBin    = 0.1;            % bin for frethist_norm
constants.synch1     = 200;            % Intensity value interval between 
                                       % which events are 
constants.synch2     = 1000;           % post-synchr:synch2 > synch1

constants.xLwLim = constants.dt;       % unit in s
constants.xUpLim = 20;                 % unit in s
constants.yLwLim = -150;
constants.yUpLim = 1000;

%---                          
%--- Parameters for 'cont2Hst.m'
%---
constants.nFramesIntgr=50;

%---
%--- Parameters for 'TDP.m'
%---
constants.tdpIntAxis  = -125:125:5500; 
constants.tdp_max     = 0.0025*1;
constants.conLevelTdp = 0:1:5; 

%---
%--- Parameters for trace statistics
%---
constants.gamma = 0.7;   
constants.min_int = 30;             % above which consider data as signal
constants.fretEventTreshold = 0.14; % FRET val. used for detect FRET events
constants.rle_min = 5;              % count alive when above min_fret for 
                                    % # frames
constants.dimerThresh = 2;          % Dimer intensity is n times abouve the
                                    % bgnd intensity
constants.minStepSize = 1.0;        % Multiplication factor to define min 
                                    % height between two pb-steps

constants.NBK=100;                  % Window size used for calc. bgd stat.

constants.NSTD=2;                   % PB detect threshold (-> CalcPBSteps)
constants.TAU=3;                    % median filter window size (PB detect)

constants.gettracesThresholdStd = 8;%see gettraces.m

constants.overlap_nstd=5;           % multiple PB detection threshold
constants.blink_nstd=4;             % set FRET=0 below threshold (donor is 
                                    % blinking)
constants.nStdJumps=4;              % set upper/lower limit for particle 
                                    % acceleration (moving particles) 
constants.limJumps=10;              % set upper/lower limit for particle 
                                    % acceleration (resting particles)
constants.intMonomer=1000;          % Monomer Intensity

%---
%--- Parameters for 1D-Canny Edge Detector
%---
constants.scales     = [1, 3];      % Gauss Kernel
constants.threshold  = [.2,.5];     % Threshold 

%---
%--- Parameters for smCellColocalization
%---
constants.nnLimit  = 125;           % below this limit two neighboring 
                                    % molecules are counted as colocalized
                                    % Unit nm; 
                                    % Rayleigh Criterion R=lambda/2*NA
%---
%--- Parameters for calculating the mapping function from bead-grid 
%---                                    
constants.gridSize  = 1.5;          % unit um                      
constants.rows      = 26;           % number of rows in the grid
constants.cols      = 24;           % number of cols in the grid
constants.ampThreshRed = 0.001;     % Detect red beads only above this thresh.
constants.ampThreshGre = 0.001;     % Detect green beads only above this thresh.
constants.xOrigin  = 1;             % roiDataMapping(1)
constants.yOrigin  = 128;           % roiDataMapping(1)
constants.imWidth  = 256;           % roiDataMapping(1)
constants.imHeight = 150;           % roiDataMapping(1)
constants.xTranslation = 0;         % correction of a xCoordinate shift 
constants.yTranslation = 0;         % correction of a yCoordinate shift

end

