%--------------------------------------------------------------------------
%
% chkCorrections:
%  Error checking script for the alphaCorrCell, deltaCorrCell and 
%  gammaCorrCell functions in smCellFret.  
% 
% Description:
%   Script that calculates corrected FRET values from raw donor and acceptor  
%   data using equation 20 ( online methods [1] ) and compares the results
%   with the FRET values calculated with the alpha, delta and gamma 
%   correction functions used in smCellFret.
%   [1] Hellenkamp, B., et al. (2018). "Precision and accuracy of 
%   single-molecule FRET measurements-a multi-laboratory benchmark study." 
%   Nature Methods 15(9): 669-676.
% 
% Syntax:  
%   chkCorrections
% 
% Inputs:
%   Update the file names of the type 'fretTracesFif*.traces' in code  
%   section: 'Load Results from the smCellFRET pipeline into workspace' 
%   1. fretTracesDif*.traces file - uncorrected traces file
%   2. fretTracesDif*_a.traces file - alpha corrected traces file 
%   3. fretTracesDif*_a_d.traces file - alpha and delta corrected traces
%      file
%   4. fretTracesDif*_a_d_g.traces file - alpha, delta and gamma corrected
%      traces file
%
%   In code section: 'Specify range and correction factors'
%   1. set the range in which the corrections should be mathematically
%      compared
%   2. specify the correction factors
% 
% Outputs:
%   The script outputs the following workspace variables:
%   1. diff_1 - Difference map, in which the Fret-values calculated from 
%      the raw data are compared with the Fret-values in smCellFRET.
%   2. diff_2 - Difference map, in which the alpha corrected Fret-values 
%      calculated from the raw data are compared with the alpha corrected 
%      Fret-values in smCellFRET.
%   3. diff_3 - Difference map, in which the alpha & delta corrected 
%      Fret-values calculated from the raw data are compared with the alpha 
%      & delta corrected Fret-values in smCellFRET.
%   4. diff_4 - Difference map, in which the alpha & delta & gamma corrected 
%      Fret-values calculated from the raw data are compared with the alpha 
%      & delta & gamma corrected Fret-values in smCellFRET.
%     
% 
% See also: 
%   alphaCorrectCell, deltaCorrectCell, gammaCorrectCell
%
% Author: 
%   - P.G. Sep 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% Load results from the smCellFRET pipeline into workspace

clear;
% sync traces + total filter
flt  =loadTraces('fretTracesDif_allBst_Fre02_sync_flt.traces');

% sync traces + total filter + crosstalk corr
flt_a=loadTraces('fretTracesDif_allBst_Fre02_sync_flt_a.traces');

% sync traces + total filter + crosstalk + direct extn corr
flt_d=loadTraces('fretTracesDif_allBst_Fre02_sync_flt_a_d.traces');

% sync traces + total filter + crosstalk + direct extn + gamma corr
flt_g=loadTraces('fretTracesDif_allBst_Fre02_sync_flt_a_d_g.traces');


%% Define range and correction factors

nTraces  = 1:10; % range of traces
nFrames = 1:5;   % range of frames

alpha = 0.09417; % Crosstalk Factor
delta = 0.056;   % Direct Excitation Factor
gamma = 0.83;    % Donor Scaling Factor

%% raw Fret

acc = flt.acceptor(:,:);
don = flt.donor(:,:);
nMol = size(acc,1);

fret    = acc(nTraces,nFrames) ./(don(nTraces,nFrames) + acc(nTraces,nFrames));

% Check: calculate difference map
diff_1 = round(flt.fret(nTraces,nFrames)-fret,4);

%% crosstalk correction

acorr = alpha*don(nTraces,nFrames);
acc_a = acc(nTraces,nFrames)-acorr;

fret_a  = (acc(nTraces,nFrames) - acorr )./...
          (don(nTraces,nFrames) + (acc(nTraces,nFrames) - acorr));
      
% Check: calculate difference map
diff_2 = round(flt_a.fret(nTraces,nFrames)-fret_a,4);

%% direct excitation correction 

dExtn = [flt_d.traceMetadata.directExtn];
total = dExtn(2:3:3*max(nTraces))';
total = repmat(total,1,max(nFrames));
dcorr = delta.*total;
acc_d = acc(nTraces,nFrames) - acorr - dcorr;

fret_d  = (acc(nTraces,nFrames) - acorr - dcorr )./...
          (don(nTraces,nFrames) + (acc(nTraces,nFrames) - acorr - dcorr));
% fret_d  = acc_d./(don(nTraces,nFrames) + acc_d);

% Check: Calculate difference map
% in the pipeline fret is set to 0 if if the value is out of range. In the
% control calculation this is not the case. 
diff_3 = round(flt_d.fret(nTraces,nFrames)-fret_d,4);      

%% gamma correction

gcorr = gamma.*don(nTraces,nFrames);

fret_g  = (acc(nTraces,nFrames) - acorr - dcorr )./...
          (gamma.*don(nTraces,nFrames) + (acc(nTraces,nFrames) - acorr - dcorr));
      
%fret_g  = acc_d./(gcorr + acc_d);
      
% Check: calculate difference map
diff_4 = round(flt_g.fret(nTraces,nFrames)-fret_g,4);
        
%%
     