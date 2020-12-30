%--------------------------------------------------------------------------
%
% scriptContourPlot:
%   Generate a contour plot and contour profile from single molecule FRET 
%   traces data. USE 'scriptContourPlot_flt.m' for batch processing.   
% 
% Description:
%   Running the script will open the file explorer and prompt the user to 
%   select a *.traces files (located in #Ch2 directories). The script 
%   uses the SPARTAN [1] function 'makecplot.m' and generates two ASCII data
%   files which will be saved under the same file name as the input file 
%   with the additional extension *._2dHst.txt and *._1dHst.txt. ASCII 
%   text files can be imported into a data analysis software to generate a 
%   2d-FRET contour plot and 1d-FRET histogram. 
% 
% Syntax:  
%   scriptContourPlot
% 
% Inputs:
%   1. *_sync*.traces files - Fret-traces files that have been 
%      post-synchronized using 'scriptPostSyncTraces.m'. The script prompts 
%      the user to select a traces-file.
%   2. 'integrWindwow' parameter (see code section 'Load Options') 
%      - An integer value n that defines the time interval [0, n] (unit of  
%      n is frames) of the 2d-FRET histogram (pseudo-color contour plot).   
%      The 1d-FRET histogram is calculated by integration over the interval  
%      [0, n] of the 2d-FRET histogram.  
% 
% Outputs:
%   1. '*1dHst.txt' file - containg the 1d FRET contour profile data
%   2. '*2dHst.txt' file - containg the 2d FRET contour plot data
% 
% Other m-files required: 
%   Subfunctions: SPARTAN functions makecplot.m and cascadeConstants.m [1] 
% 
% See also: 
%   scriptContourPlot_flt.m, scriptPostSyncTraces.m
%
% References:
%   [1] Juette, M. F., et al. (2016). "Single-molecule imaging of 
%   non-equilibrium molecular ensembles on the millisecond timescale." 
%   Nature Methods 13(4): 341-344.
%
% Permissions:
%
% Authors: 
%   - P.G. 2018
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Generate Contour Plot and Contour Profile

% Get filename
clear;

[fileName,pathName,idxCancel] = uigetfile('*.traces',...
                                            'Select fretTraces*.traces file');
tracesFile = [pathName fileName];

if idxCancel>0 && exist(tracesFile,'file')   
    tracesFile  = [pathName fileName];
    data  = loadTraces(tracesFile);
    nMol = size(data.donor,1);
    disp(['nTraces = ', num2str(nMol)]);
else
    warning('no file available')
    return
end


%% Load options
constants = cascadeConstants;% cascade constants from spartan v5.3.1 
options   = constants.defaultMakeplotsOptions;

% normalize contour plot
contourLength = 125; %unit is frames

% normalize to the number of traces
options.cplot_normalize_to_max = false;
cHst2D = makecplot(tracesFile,options); % makecplot from spartan v5.3.1 
                                        % cHst2D is normalized by number of traces
                                        % the sum of every colum in cHst2D=1 

%% Calculate Contour Profile and nTraces 
% In this normalization the sum of cHst1D gives 1
popHist_contourLength = nansum(cHst2D(2:end,2:contourLength+1),2);
cHst1D  = [ cHst2D(2:end,1) 100*(popHist_contourLength/contourLength)];

% % % normalize to the number of frames
% % % options.cplot_normalize_to_max = 4000;
% % % cHst2D = makecplot(filename,options); %makecplot from spartan v5.3.1 
% % % popHist_nTraces = sum(cHst2D(2:end,2:end),2);
% % % cHst1D  = [ cHst2D(2:end,1) popHist_nTraces];

% Cntrl Step: Integrate the Histogram to get the number of traces or ...
% frames (contour length)
sum1DHst = sum(cHst1D(1:end,2));
disp(['sum1DHst :' num2str(sum1DHst)]);

%% Save Histogram data to txt file

% Input Filename for 2d-Histogram
name2dHst  = inputdlg(...
             {'Enter File Name 2dHst:'},...              %prompt
             'Save File As: ',...                        %dlg titel
             [1 50],...                                  %width of line
             {strrep(fileName,'.traces','_2dHst.txt')}); %default answer

if isempty(name2dHst)
    return
else
    disp([pathName name2dHst{1}]);
    dlmwrite([pathName name2dHst{1}],cHst2D);
end

% Input Filename for 1d-Histogram
name1dHst  = inputdlg(...
             {'Enter File Name 1dHst:'},...              %prompt
             'Save File As: ',...                        %dlg titel
             [1 50],...                                  %width of line
             {strrep(fileName,'.traces','_1dHst.txt')}); %default answer

if isempty(name1dHst)
    return
else 
    disp([pathName name1dHst{1}]);
    dlmwrite([pathName name1dHst{1}],cHst1D);
end

%% End
disp('Finished...')

%%


