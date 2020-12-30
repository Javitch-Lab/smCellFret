%--------------------------------------------------------------------------
%
% scriptGetSegTrackingAnalysis:
%   Script that performs a DC_MSS transient motion analysis. 
% 
% Description:
%   This script is a batch processing script that applies the function 
%   'basicTransientDiffusionAnalysisv1' to multiple user selected tracking 
%   data. The function 'basicTransientDiffusionAnalysisv1.m' is part  
%   of the Divide-and-Conquer Moment Scaling Spectrum transient diffusion
%   analysis (DC-MSS) [1]
%
%   [1] Vega, A. R., et al. (2018). "Multistep Track Segmentation and 
%   Motion Classification for Transient Mobility Analysis." Biophysical 
%   Journal 114(5): 1018-1025.
% 
% Syntax:  
%   scriptGetSegTrackingAnalysis
% 
% Inputs:
%   The file explorer prompts the user to select a list of directories. Each 
%   selected directory need to contain the file name 'Tracking*.mat'.
%   a) for post filtered data 'TrackingPst.mat' (default)
%      Note: select only #Ch1 directories that contain the data file
%      'TrackingPst.mat'
%   b) for raw data 'Tracking.mat' Note: Change the default file name in
%      code section 'Get Results from the DC-MSS Analysis' 
% 
% Outputs:
%   1. Saves the output of the DC-MSS analysis in the results.mat file 
%      located in the current batch file folder.(For a description of 
%      results.mat file see header of basicTransientDiffusionAnalysisv1.m)
%   2. Figure with color-coded track segments. The colors indicate 
%      different diffusion modes (brown: immobile, blue: confined diffusion, 
%      cyan: normal diffusion, magenta: super diffusion, black: 
%      unclassified).
% 
% 
% Other m-files required: 
%   Subfunctions: getdirs.m, basicTransientDiffusionAnalysisv1.m
% 
% See also: 
%   scriptGetSegResults.m,  scriptDiff2Frettraces.m
%
% Author: 
%   - S. Mathiasen 2018
%   - P.G. Jul 2018: Modified script for handling TrackingPst files where
%     the variable tracksFinal is empty: assign empty to the structure 
%     variable 'results'.
%
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Get Results from the DC-MSS Analysis
clear;
dirs=getDirs('Highlight and Open ACCEPTOR (#Ch1) Folders (HIT CANCEL TO END SELECTION)'); 
folders=numel(dirs(1,:));

h = waitbar(0,'Processing Data...');

for n=1:folders
   
    % load Tracking.mat -> tracksFinal and get results
    fn=[dirs{2,n} filesep dirs{1,n} filesep 'TrackingPst.mat']; 
    load(fn,'tracksFinal'); 
    if ~isempty(tracksFinal(1).tracksCoordAmpCG)
        results = basicTransientDiffusionAnalysisv1(tracksFinal);
    else
        results = [];
    end
    
    % Save results
    filename=[dirs{2,n} filesep dirs{1,n} filesep 'results.mat' ];
    save(filename,'results');
    
    % Save figure
    figname=['figCell' dirs{1,n}];
    savefig(gcf,[dirs{2,n} filesep dirs{1,n} filesep figname]);
    
    % Show Waitbar
    waitbar(n / folders)

end

close(h) 

%% End



