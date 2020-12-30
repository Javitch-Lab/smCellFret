%--------------------------------------------------------------------------
%
%  scriptGetFretTraces:
%   Batch processing script that extracts donor and acceptor fluorescence 
%   time traces from the u-track output files detectionMaxProject.mat and 
%   Tracking.mat.  
% 
% Description:
%   To extract the donor and acceptor time traces from tracks, the script
%   performs the following steps in order: 
%    1. Get Ch1 and Ch2 directory pairs from user.
%    2. Calculate mapping function. 
%    3. Import donor and acceptor tracks from u-track and assign  
%       identifiers to each track. 
%    4. Create intensity baselines before and after each track. Each track 
%       consists of three regional segments: Baseline positions before the 
%       start of tracking, the positions of the tracked particle from the 
%       start to the termination of the tracker and the baseline positions 
%       after photobleaching or track terminateion.  
%    5. Determine the local donor and acceptor spot intensity at each time 
%       point along the track and baselines. Each spot intensity is 
%       calculated as the sum of the pixel intensities in a 5x5 rectangular 
%       array centered around the particle position at time t.
%    6. Apply the inverse geometric transformation of tform to the u,v 
%       (imCh1, Acceptor) coordinates and calculate the transformed 
%       coordinates imCh1_trans (x,y). The transfromed coordinates can now  
%       be colocalized with the coordinates in imCh2(Donor)
%    7. Subtract the background signal from the spot intensities to bring 
%       the baseline of the donor and acceptor fluorescence close to zero.
%    8. Calculate statistical measures for each trace and update the field 
%       traceMetadata.         
%    9. Continuation of the corresponding donor trajectory. After the loss
%       of the acceptor signal for each particle, the mapped donor 
%       trajectory was concatenated with the nearest matching donor 
%       trajectory within the search radius at the position of the acceptor
%       loss.  
%   10. Compile data into a fretTraces structure. 
% 
% Syntax:  
%   scriptGetFretTraces.m
% 
% Inputs:
%   The script prompts the user to select the calibration file 
%   (‘Spots in tracks statistics.txt’, see fct imageRegistrationTrackMate)
%   which contains the coordinates of the control points required to infer 
%   the transformation function. If the user clicks Cancel the program will
%   continue without using a transformation function (no mapping) and 
%   prompts the user to select pairs of #Ch1 – #Ch2 experiment folders. 
%   After all desired pairs of folders are selected, click Cancel to close 
%   the File Explorer and continue running the MATLAB script. Not providing 
%   donor-acceptor pairs of experiment folders will lead to errors. 
% 
% Outputs:
%   The script outputs various files (files with extension *.traces) most
%   of which are used to track down possible errors in the pipeline. The 
%   final results are saved in a structure array called fretTraces which is 
%   saved as a .mat file in the #Ch2 experiment folder. 
% 
%   The structure varibale fretTraces has the following two main fields:
%   Ch1 - acceptor data in emission channel 1 
%   Ch2 - donor data in emission channel 2
%   Each emission channel Ch has the following subfields:
%   .time - time in units of seconds
%   .x     - x coordinate of the localized spot 
%   .y     - y coordinate of the localized spot 
%   .xCorr - x coordinate after applying the transformation function
%   .yCorr - y coordinate after applying the transformation function
%   .int   - Integrated local spot intensities (locInt, the unit is photons  
%            / frame). The integration was performed over an area of 5x5
%            pixels covering the spot's point spread function (psf). 
%   .snr   - The local signal-to-background ratio is calculated by dividing 
%            the integrated spot intensity by the local background signal.
%            The local background signal (locBgd) is calculated by 
%            integrating the pixel intensities over the adjacent area of 
%            the psf (see fct. nPhotons.m). 
%   .traceMetadata has the following subfields:
%     .ids          - Trace identifier 
%     .lenBackset   - Length of the baseline before tracking starts (in 
%                     frames). 
%     .startOfTrace - Starting point of tracking.
%     .endOfTrace   - End point of tracking. 
%     .lenBaseline  - Length of the baseline in frames. The baseline begins  
%                     one frame after the point where the tracking ended.
%     .traceLen     - Track length equals endOfTrack - startOfTrack + 1   
%     .meanInt      - Mean value of local acceptor intensity (locIntAcc) 
%                     minus mean value of local acceptor background  
%                     (locBgdAcc). See above '.imt' and '.snr' for the 
%                     definition of locInt and locBgd.
%     .meanXIntAcc  - locIntAcc minus the mean value of the extended local 
%                     acceptor background (locXBgdAcc). The extended local  
%                     background is used for traces that have a length of  
%                     less than 10 frames (see function statFret.Traces.m)
%     .meanXIntDon  - mean of locIntDon minus mean of locXBgdDon. 
%     .meanXSnrAcc  - mean of locIntAcc diveded by mean of locXBgdAcc
%     .meanXSnrDon  - mean of locIntDon diveded by mean of locXBgdDon
%     .meanXCrt     - meanXIntAcc diveded by meanXIntDon
%
% Other m-files required: 
%   Subfunctions: imageRegistrationTrackMate.m, utracks2traces.m, 
%   makeBaseline.m, nPhotons.m, loadTracesCell.m, saveTracesCell.m,
%   correctBgnd.m, statFretTraces.m, findDonorTracks.m, traces2mat.m
%
% Authors: 
%   - P.G. Sep 2019
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% ========================================================================
%=== Initialize 
%==========================================================================

clear;
constants = smCellConstants;

% If donorTracking is true each mapped donor track will be concatenated 
% with the nearest matching donor track within the search radius at the 
% acceptor photobleaching position to give a joint track of the 
% observed donor signal.
donorTracking = true; 

% Use extended local background
flagExtLocBgd = true;

% Load calibration file
[fileName,pathName,idxCancel] =...
    uigetfile('*.txt',...
              'Select calibration file (Spots in tracks statistics.txt -> Go to folder: examples\calibration (HIT CANCEL TO END SELECTION)');
          
calFile = [pathName fileName];

if idxCancel>0 && exist(calFile,'file')   
    calFile  = [pathName fileName];
    % Use mapping
    calcImTransform = true;
    disp(calFile);    
    
else
    warning('no calibration file available')
    % No mapping
    calcImTransform = false;
end
                                       

%% ========================================================================
%=== Get data directories from user
%==========================================================================
dirArray = getDirs('Highlight and Open Pairs of #Ch1 and #Ch2 Folders (HIT CANCEL TO END SELECTION)');
if isempty(dirArray); return; end 

%% ========================================================================
%=== Calculate mapping function
%==========================================================================
% Only use the centerpart of the image emission splitter
if calcImTransform
    
    roiDataMapping(1) = constants.xOrigin; % xOrigin
    roiDataMapping(2) = constants.yOrigin; % yOrigin
    roiDataMapping(3) = constants.imWidth; % imWidth
    roiDataMapping(4) = constants.imHeight;% imHeight
    rows = constants.rows; % nrows of the calibration grid
    cols = constants.cols; % ncols of the calibration grid

    calcTRE   = true;
    [tform_lwm, ~, TRE] = imageRegistrationTrackMate(calFile,calcTRE,...
                            rows,cols,roiDataMapping);

    % Display target registration error
    if ~isempty(TRE);disp(['TRE :' num2str(TRE) ' nm']);end
    
end


%% ========================================================================
%=== Get Acceptor – Donor Intensity Time Traces
%==========================================================================

% Initialize
j=0;ch1Dirs=cell(1,1);

% Get the number of Ch1 directories
for i=1:numel(dirArray(1,:))
    if ~isempty(cell2mat(strfind(dirArray(1,i),'Ch1')))
        j=j+1;
        ch1Dirs{1,j} = dirArray{1,i}; % experiment
        ch1Dirs{2,j} = dirArray{2,i}; % path
    end
end

nExperiments=length(ch1Dirs(1,:));

% Loop over all Fret experiments
for j=1:nExperiments
    
    %----------------------------------------------------------------------
    %--- Import Tracks and add ids to each trace
    %----------------------------------------------------------------------
     disp('import tracking data...');
    % Set experiment folder
    dataSet    = ch1Dirs{1,j};
    basepath   = ch1Dirs{2,j};
    currentDir = [basepath '\' dataSet '\'];
    
    % Get number of movie frames
    if j==1
        tmpVar=load([currentDir 'detectionMaxProject.mat']);
        allFrames = length(tmpVar.movieInfo);
    end
    
    % if nChannelsTracked = 2: Concatenate mapped donor tracks with the 
    % nearest matching donor track after acceptor photobleaching
    nChannelsTracked = 1 + donorTracking;
    
    for k=1:nChannelsTracked
        
        chName={'Ch1','Ch2'};
        dataSetTmp = strrep(dataSet,'Ch1',chName{k});
        currentDir = [basepath '\' dataSetTmp '\'];

        % Import traces tracked with utrack
        disp('converting utrack traces...');
        utracks2traces([currentDir 'Tracking.mat'], ...
            constants.dt, allFrames);

        % Rename file in Ch1
        % (uTrack2traces saves data to a file 'singleTracesCh1.traces') 
        if k==1 
            oldFile=[currentDir 'singleTracesCh1.traces'];
            newFile=[currentDir 'singleTracesAcc.traces'];
            movefile(oldFile, newFile);

        % Rename file in Ch2
        elseif k==2
            oldFile=[currentDir 'singleTracesCh1.traces'];
            newFile=[currentDir 'singleTracesDonB.traces'];
            movefile(oldFile, newFile);
        end
    end
    
    %----------------------------------------------------------------------
    %--- Create Baseline for each trace
    %----------------------------------------------------------------------
    disp('makeBaseline...');
    baseline  = cell(nChannelsTracked,1);
    for k=1:nChannelsTracked
        chName       = {'Ch1','Ch2'};
        dataSetTmp   = strrep(dataSet,'Ch1',chName{k});
        currentDir   = [basepath '\' dataSetTmp '\'];
        filename     = 'singleTracesAcc.traces';
        if k==1      %Ch1 baseline of tracked acceptor traces (Acc)
            baseline{k} = makeBaseline([currentDir filename],...
                                        constants.lenBackset,...
                                        constants.lenBaseline,...
                                        constants.stdBgd);
    
        elseif k==2 %Ch2 baseline of tracked donor traces (DonB)
            % make baseline for donor traces tracked in Ch2
            filename = strrep(filename,'Acc','DonB');
            baseline{k} = makeBaseline([currentDir filename],...
                                        constants.lenBackset,...
                                        constants.lenBaseline,...
                                        constants.stdBgd);
        end
    end
    
    %----------------------------------------------------------------------
    %--- Calculate number of Photons in Ch1 (Acceptor) and Ch2 (Donor)
    %----------------------------------------------------------------------
    for k=1:nChannelsTracked
        chName       = {'Ch1','Ch2'};
        dataSetTmp   = strrep(dataSet,'Ch1',chName{k});
        currentDir   = [basepath '\' dataSetTmp '\'];
        filename     = 'singleTracesAcc_bln.traces';
        if k==1     %Ch1, tracked acceptor traces (Acc)
            disp('nPhotons Acc...');
            [~,~,selListCh1] = nPhotons([currentDir filename],...
                                        constants.QE);
        elseif k==2 %Ch2, tracked donor traces (DonB)
            disp('nPhotons DonB...');
            filename = strrep(filename,'Acc','DonB');
            [~,~,selListCh2b] = nPhotons([currentDir filename],...
                                        constants.QE);
        end
    end
    
    %----------------------------------------------------------------------
    %--- Apply the inverse geometric transformation of tform to the u,v 
    %--- (imCh1, Acceptor) coordinates and calculate the transformed 
    %--- coordinates imCh1_trans (x,y). The transfromed coordinates can now  
    %--- be colocalized with the coordinates in imCh2(Donor)
    %----------------------------------------------------------------------
    currentDir = [basepath '\' dataSet '\'];
    data1   = loadTracesCell([currentDir 'singleTracesAcc_bln.traces']);
 
    x=data1.x;
    y=data1.y;
    [nTraces,nFrames]=size(data1.x);
    
    if calcImTransform
        disp('Calc Transformed Coordinates...');
        xCorr=zeros(nTraces,nFrames);
        yCorr=zeros(nTraces,nFrames);
       
        parfor i=1:nFrames
            imCh1=[x(:,i), y(:,i)];

            % Calculate inverse transformation (switch off warnings)
            imCh1_trans=transformPointsInverse(tform_lwm,imCh1);
            warnStruct = warning('query', 'last');
            
            if ~isempty(warnStruct)
                msgid_tforminv = warnStruct.identifier;
                warning('off', msgid_tforminv);
            end

            % Check for incorrect registered datapts
            indX=find(abs(imCh1(:,1)-imCh1_trans(:,1))>10); 
            if ~isempty(indX)
                for k=1:numel(indX)
                    imCh1_trans(indX(k)) = imCh1(indX(k));
                end
            end

            % Set negativ entries to zero
            tmpX=imCh1_trans(1:nTraces,1);
            tmpY=imCh1_trans(1:nTraces,2);
            idx=tmpX<0;
            tmpX(idx)=0;
            tmpY(idx)=0;

            % Save corrected values
            xCorr(:,i) = tmpX;
            yCorr(:,i) = tmpY;  
        end
    end
    
    %----------------------------------------------------------------------
    %--- Save transformation results in data structure
    %----------------------------------------------------------------------
    if calcImTransform == 1
        % Save the corrected coordinates in Don and Acc file
        data2 = loadTracesCell([currentDir 'singleTracesAcc_pht.traces']);
        data2.channelNames = {'x','y','xCorr','yCorr','int','snr','locBgd'};
        data2.xCorr = xCorr; 
        data2.yCorr = yCorr;
        outfile = 'singleTracesAcc_pht.traces';
        saveTracesCell([currentDir outfile],'traces',data2);
        
        data1.x = xCorr; 
        data1.y = yCorr;
           
    end
    
    dataSetTmp = strrep(dataSet,'Ch1','Ch2');
    currentDir = [basepath '\' dataSetTmp '\'];
    outfile = 'singleTracesDonA_bln.traces'; %(DonA)
    saveTracesCell([currentDir outfile],'traces',data1);
    
    %----------------------------------------------------------------------
    %--- Calculate number of Photons in Ch2 (DonA)
    %----------------------------------------------------------------------
    disp('nPhotons DonA...');
    filename='singleTracesDonA_bln.traces';
    [~,~,selListCh2] = nPhotons([currentDir filename],...
                                constants.QE);
    
    %----------------------------------------------------------------------
    %--- Synchronize the selection lists (DonA)
    %----------------------------------------------------------------------
    disp('synchronize...')
    % Chk if tracks in selListCh2 are found in selListCh1 
    idx=ismember(selListCh2,selListCh1);
    selList=selListCh2(idx);
    
    %----------------------------------------------------------------------
    %--- Correct trace background in Ch2: mapped donor traces (DonA)
    %----------------------------------------------------------------------
    disp('correct trace bgd...');
    filename='singleTracesDonA_pht.traces';
    selTraces([currentDir filename],selList,1,[currentDir filename]);
    correctBgnd([currentDir filename],baseline{1}(selList',:));
    
    %----------------------------------------------------------------------
    %--- Correct trace background in Ch1 and Ch2 (Ch2, if donor is tracked)
    %----------------------------------------------------------------------
    for k=1:nChannelsTracked
        chName       = {'Ch1','Ch2'};
        dataSetTmp   = strrep(dataSet,'Ch1',chName{k});
        currentDir   = [basepath '\' dataSetTmp '\'];
        filename     = 'singleTracesAcc_pht.traces';
        if k==1     %Ch1: tracked acceptor traces (Acc)
            selTraces([currentDir filename],selList,1,[currentDir filename]);
            correctBgnd([currentDir filename],baseline{k}(selList',:));
        elseif k==2 %Ch2: tracked donor traces (DonB)
            filename = strrep(filename,'Acc','DonB');
            selTraces([currentDir filename],selListCh2b,1,[currentDir filename]);
            correctBgnd([currentDir filename],baseline{k}(selListCh2b',:));
        end
    end
    
    %----------------------------------------------------------------------
    %--- Search for corresponding donor traces
    %----------------------------------------------------------------------
    try
        % Load Donor / Acceptor data
        currentDir=cell(2,1);
        for k=1:2
            chName       = {'Ch1','Ch2'};
            dataSetTmp   = strrep(dataSet,'Ch1',chName{k});
            currentDir(k)= {[basepath '\' dataSetTmp '\']};

            if k==2 && nChannelsTracked==1
                % ---------------------------------------------------------
                % --- Filter traces (no donor tracking)
                % --------------------------------------------------------- 
                % Use the result (filter Idx) of the acceptor filtering
                disp('Calc trace Statistics...');
                donFile1 = 'singleTracesDonA_phc.traces';
                accFile  = 'singleTracesAcc_phc.traces';
                statFretTraces([currentDir{1} accFile],...
                               [currentDir{2} donFile1],...
                               flagExtLocBgd);

            elseif k==2 && nChannelsTracked==2
                
                % ---------------------------------------------------------
                % --- Filter traces (with donor tracking)
                % ---------------------------------------------------------
                disp('Calc trace Statistics...');
                donFile1 = 'singleTracesDonA_phc.traces';
                accFile  = 'singleTracesAcc_phc.traces';
                statFretTraces([currentDir{1} accFile],...
                               [currentDir{2} donFile1],...
                               flagExtLocBgd );

                %----------------------------------------------------------
                %--- Search for corresponding donor traces
                %----------------------------------------------------------
                disp('find Donor Traces...');
                donFile2 = 'singleTracesDonB_phc.traces';

                % Get data for calculating Fret 
                findDonorTracks([currentDir{2} donFile1],...
                                [currentDir{2} donFile2]);                                
            end
        end
    catch
        disp(currentDir);
        warning('Check Fret Filter Criteria:');
        continue
    end
end

%% ========================================================================
%=== Convert traces files
%==========================================================================

traces2mat(dirArray, donorTracking);

%% End
disp('END');    

%%






