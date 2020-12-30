
%--------------------------------------------------------------------------
%
% smCellFretDemo.m:
%   The script demonstrates the workflow of the smCellFret pipeline. 
% 
% Syntax:  
%   smCellFretDemo.m
% 
% Inputs:
%   Input files or input data are described separately in each code
%   section.
% 
% Outputs:
%   Output files or output data are described separately in each code 
%   section.
%
% Authors: 
%   - Jozsef Meszaros 2020
%   - P.G. Oct 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
%                           USER INPUT REQUIRED
% 
%             Update the name of the installation directory
%                    including the smCellFret version. 
clear; 
installationDir = 'D:\#smCellFRET_v1.3';
%% ------------------------------------------------------------------------
% 
%                            START DEMO
%
%                    Use 'Run and Advance' Mode 
%
% -------------------------------------------------------------------------

%% Add the smCellFret to Matlab path
restoredefaultpath;
addpath( genpath( installationDir ) );

%% Delete Files from previous Demo runs.

pathName = [installationDir filesep 'examples']; cd(pathName);

if exist([pathName filesep 'mGlu2 data\#06Ch1'], 'dir') == 7
    rmdir([pathName filesep 'mGlu2 data\#06Ch1'],'s');
end

if exist([pathName filesep 'mGlu2 data\#06Ch2'], 'dir') == 7
    rmdir([pathName filesep 'mGlu2 data\#06Ch2'],'s');
end

if exist([pathName filesep 'microscope data\#06Ch1'], 'dir') == 7
    rmdir([pathName filesep 'microscope data\#06Ch1'],'s');
end

if exist([pathName filesep 'microscope data\#06Ch2'], 'dir') == 7
    rmdir([pathName filesep 'microscope data\#06Ch2'],'s');
end

if exist([pathName filesep 'mGlu2 data\40msG3mW80Pd100G_06.stk'],'file') == 2
    delete([pathName filesep 'mGlu2 data\40msG3mW80Pd100G_06.stk'])
end

if exist([pathName filesep 'mGlu2 data\CombinedStats.mat'],'file') == 2
    delete([pathName filesep 'mGlu2 data\CombinedStats.mat'])
end

%% Read Section 5.1 of the smCellFret Documentation
% Get familiar with the Project File Folder Structure
clc;
%% Read Section 5.2.1 of the smCellFret Documentation
% Acquire control point data for registering donor and acceptor emission 
% channels

% Example data is located in:
disp('----- 1: Get Calibration File');
save([installationDir filesep 'installationDir.mat'],'installationDir');
cd([installationDir filesep 'examples\calibration\Grid #01']);

% The file ‘grid xy.stk’ is a stack of TIFF images which can be played as a
% movie using the image processing software ImageJ. The movie shows how a 
% pair of control points translate stepwise across the acceptor (left 
% control point) and donor emission channel (right control point) thus 
% forming a grid of 26 rows and 24 columns in each emission channel. 
% The location of each control point was measured using the ImageJ plugin
% TrackMate and saved in the same folder under the name
% ‘Spots in tracks statistics.txt’. 

%% Read Section 5.2.2 of the smCellFret Documentation 
% Acquire smFRET movies and save movie data as .tif files

% Example data is located in:
cd([installationDir filesep 'examples\microscope data']);

% The example movie '40msG3mW80Pd100G_06.stk' is a stack of TIFF images 
% which can be played as a movie using the image processing software ImageJ.
% In order to track single donor and acceptor fluorescence spots using 
% u-track the movie file needs to be first converted into a sequence of 
% individual tiff images. Second, the two spectrally separated regions in 
% each tiff image must be split and saved as a series of numbered images in
% the acceptor (#%%Ch1/ImageData) and donor (#%%Ch2/ImageData) image
% directories. Convert and split the example smFRET movie  
% '40msG3mW80Pd100G_06.stk' into TIFF images by running the script
% ‘scriptStk2Tiff.m’. 
disp('----- 2: scriptStk2tiff.m');
run('scriptStk2tiff.m'); 

% Results will be save in project folder mGlu2 data
load('installationDir.mat');
source = [installationDir filesep 'examples\microscope data'];
destination = [installationDir filesep 'examples\mGlu2 data'];
copyfile(source, destination);

% Goto project folder
cd([installationDir filesep 'examples\mGlu2 data'])

%% Read Section 5.2.3 of the smCellFret Documentation
% Process TIFF stacks using u-track to extract single molecule tracks
 
% Run the tracking software u-track and use the Tiff images stored in the 
% experiment folders ...\#06Ch1\ImageData and ... \#06Ch2\ImageData as  
% input data. 

% In this demo, tracking results are copied to the project folders 
disp('----- 3: Analyze Data with u-track');
source = [installationDir filesep 'examples\u-track results'];
destination = [installationDir filesep 'examples\mGlu2 data'];
copyfile([source '\#06Ch1'], [destination '\#06Ch1']);
copyfile([source '\#06Ch2'], [destination '\#06Ch2']);

%% Read Section 5.2.4 of the smCellFret Documentation
% Run scriptGetFretTraces to extract fluorescence time traces from single
% molecule tracks

% The script extracts donor and acceptor fluorescence time traces from the
% u-track output files detectionMaxProject.mat and Tracking.mat. Running 
% the script prompts the user to select the calibration file:
%
% ‘Spots in tracks statistics.txt’, located in: 
%  ....\examples\calibration\Grid#01
%
% This file contains the coordinates of the control points required to 
% infer the transformation function. If the user clicks ‘Cancel’ the  
% program will continue without using a transformation function (no 
% mapping) and prompts the user to select pairs of #Ch1 – #Ch2 experiment 
% folders. After all folders are selected (in the demo only one folder pair
% is available), click 'Cancel' to close the File Explorer and continue
% running the Matlab script. Final results are saved in the structure array 
% fretTraces.mat located in the #Ch2 experiment folder. 
disp('----- 4: scriptGetFretTraces.m');
run('scriptGetFretTraces.m')

%% Read Section 5.2.5 of the smCellFret Documentation
% Filter data sets to discard unwanted donor and acceptor traces

% In this demo only a lifetime filter is applied. All acceptor traces  
% with a track lifetime >=20 frames will be selected (frame rate 25Hz)
% Running this script will open the file explorer and prompts the user to
% select pairs of #Ch1 – #Ch2 experiment folders. After all desired 
% pairs of folders are selected, click 'Cancel' to close the File Explorer 
% and continue running the Matlab script. 
%
% The script outputs two new workspace files, TrackingPst.mat 
% (#Ch1 experiment folder) with the filtered tracking data and 
% fretTracesPst.mat (#Ch2 experiment folder) with the filtered FRET data.
disp('----- 5: scriptFretFilter.m');
run('scriptFretFilter.m')

%% Read Section 5.2.6 of the smCellFret Documentation 
% Analyze and classify motion within acceptor tracks (requires DC-MSS)

% Motion analysis of acceptor tracks is performed to detect diffusion 
% segments in individual tracks and classify them as free, confined,
% directed and immobile diffusion. The analysis utilizes functions
% published in the transient mobility analysis framework, termed 
% ‘‘divide-and-conquer moment scaling spectrum’’ (DC-MSS) and can be
% performed in smCellFRET by running the following scripts in sequence: 

% --- Step 1 
% The script will save a data file 'results.mat ' (Classification 
% results and MSS data for each track segment) and a figure file 
% 'figCell*.fig' showing the color-coded diffusion segments in the 
% respective Acceptor File Directory 
disp('----- 6.Step1: scriptGetSegTrackingAnalysis.m');
run('scriptGetSegTrackingAnalysis.m');

%% --- Step 2
% The script converts the structure array 'results.mat' into a matrix 
% 'SegResultsFinal.mat' and performs a statistical analysis on the four 
% segment populations: free, confined, directed and immobile diffusion. 
% The results of the statistical analysis are saved in the data file 
% 'SegStats.mat' 
disp('----- 6.Step2: scriptGetSegResults.m');
run('scriptGetSegResults.m');   

%% --- Step 3
% In the last step of the motion analysis track segments which are 
% classified as either free, confined, directed or immobile diffusion will
% be linked to the smFRET data. Running the script opens the file explorer
% and prompts the user to select pairs of Acceptor and Donor File
% Directories (#Ch1 and #Ch2 directories). 
% Updated smFRET data containing motion dynamics which are linked to FRET
% Traces are saved as data files of type fretTracesDif.mat 
disp('----- 6.Step3: scriptDiff2Frettraces.m');
run('scriptDiff2Frettraces.m');
load('installationDir.mat');
delete([installationDir filesep 'examples\mGlu2 data\#06Ch2\fretTraces.mat']);
delete([installationDir filesep 'examples\mGlu2 data\#06Ch2\fretTracesPst.mat'])

%% Read Section 5.2.7 of the smCellFret Documentation
% Visualize individual molecule tracks and intensity time traces using 
% cellFretViewtraces

% The program cellFretViewtraces displays the motion of individual 
% donor-acceptor-labeled molecules, their tracks with their corresponding 
% intensity and FRET time traces and allows the user to manually correct
% and sort traces into the three different categories.

% --- Step1
% Sort races manually. In this demo pre-sorted traces are copied to the
% project demo folder  
disp('----- 7.Step1: Manually Sort Traces');
cd([installationDir filesep 'examples\mGlu2 data\#06Ch2']);
source = [installationDir filesep 'examples\sorted data'];
destination = [installationDir filesep 'examples\mGlu2 data\#06Ch2'];
copyfile(source, destination);

%% --- Step2 
% The gui 'cellFretViewtraces' is opened. Load the file fretTracesDif.mat
% (located in the folder ...\examples\mGlu2 data\#06Ch2) and view a 
% selection of traces from the demo data set. 


disp('----- 7.Step2: cellFretViewtraces.m');
run('cellFretViewtraces.m')

%% End of Demo
