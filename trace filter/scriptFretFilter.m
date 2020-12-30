%--------------------------------------------------------------------------
%
% scriptFretFilter.m:
%   The script applies a selected set of filters to a fretTraces structure 
%   and removes unwanted donor and acceptor traces.
% 
% Description:
%   Optional file header info
% 
% Syntax:  
%   scriptFretFilter
% 
% Inputs:
%   1. The script prompts the user to select one or more pairs of #Ch1/#Ch2
%      project folders (select e.g. #05Ch1 & #05Ch2; #06Ch1 & #06Ch2; ...)
%
%   2. crtValue: parameter that corrects for donor leakage into the 
%      acceptor emission channel (change this value in in the code section 
%      "Parameters" below). The default value is 0 (no correction is 
%      applied).
%
%   3.Select filters and set filter limits by updating the variable 
%     'fltPrm' in the code section 'Local Callback Functions'. The 
%     structure variable 'fltPrm' has different property fields that 
%     activate the filter and set the filter parameters:
%       .[property]FilterOn = enables (1) or disables (o) the filter. 
%       .min[property] = minVal (sets the lower limit of the filter value) 
%       .max[property] = maxVal (sets the upper limit of the filter value)
%
%   Note: A trace is selected if its property p is in the range:
%   minVal < p < maxVal. All activated filters are linked with AND 
%   logic, i.e. a trace must pass all filter criteria to be accepted. 
%      
%   The following filters are available (specified numbers are default 
%   values):
%   1. Signal2BgdRatio Filter 
%      fltPrm.snrFilterOn = true;
%      fltPrm.minSnrVal   = 1.4;
%      fltPrm.maxSnrVal   = 2.2;
%      Time averaged ratio of the acceptor intensity of the diffraction
%      limited single molecule spot and its bordering background intensity. 
%      Time-averaging is performed over the range of the acceptor track
%      time.
%   2. Crt-Filter
%      fltPrm.crtFilterOn = true;
%      fltPrm.minCrtVal   = 0.3;
%      fltPrm.maxCrtVal   = 10;
%      Average Crosstalk Value: Ratio of the time averaged acceptor 
%      intensity of the diffraction limited spot and its corresponding time 
%      averaged donor intensity. The lower limit of the crosstalk filter 
%      (crt_min) can be calculated from the ratio of the acceptor intensity 
%      and the donor intensity measured after acceptor photobleaching. This 
%      value estimates the spectral bleed-through of the donor signal into 
%      the acceptor emission channel. Applying this filter removes  
%      primarily acceptor traces with low intensity.
%   3. Acc-Intensity
%      fltPrm.intAccFilterOn = false;
%      fltPrm.minIntAcc   = 140;
%      fltPrm.maxIntAcc   = inf;
%      Time averaged acceptor intensity: Time-averaging is performed over 
%      the range of the acceptor track time.
%   4. Don-Intensity
%      fltPrm.intDonFilterOn = false;
%      fltPrm.minIntDon   = 160;
%      fltPrm.maxIntDon   = 800;
%      Time averaged donor intensity: Time-averaging is performed over the 
%      range of the donor track time.
%   5. Total Intensity Filter
%      fltPrm.totalFilterOn = true;
%      fltPrm.minTotalVal = 450;            
%      fltPrm.maxTotalVal = inf;
%      Time averaged total intensity: The total intensity is the sum of the 
%      donor and acceptor intensity. Time-averaging is performed over the 
%      range of the donor track time.
%   6. Average FRET Filter
%      fltPrm.avgFretFilterOn = false;
%      fltPrm.minAvgFret  = 0.1;
%      fltPrm.maxAvgFret  = 1.0;
%      Time averaged Fret efficiency. Time-averaging is performed 
%      over the range of the acceptor track time.
%   7. Lt-Filter (filter value in frames)
%      fltPrm.ltFilterOn = false;
%      fltPrm.minLtVal   = 20; 
%      All tracks with a track lifetime >= 20 will be accepted
%   8. ME-Filter (8)
%      fltPrm.meFilterOn = true;
%      fltPrm.trackDist   = 1.5;
%      fltPrm.nEvents     = 25;
%      fltPrm.offTimeLim  = 0.75;
%      Multiple Event Filter: If a trajectory reappears 25 times 
%      (nEvent = 25) at the same location (i.e. centroid distance d < 1.5
%      pixels) AND the mean track lifetime is less than t_limit (unit of 
%      t_limit is s) then the analysis returns: logical 0 (trajectory is 
%      rejected). Otherwise the filter returns logical 1 ( trajectory is 
%      selected).
% 
% Outputs:
%   1. Updated structure variable fretTraces which is saved in the file 
%      fretTracesPst.mat. 
%   2. Updated u-track structure variable tracksFinal which is saved in the
%      file TrackingPst.mat.
% 
% Other m-files required: 
%   Subfunctions: getSegments.m, repeatingEventsFilter.m
% 
% See also: 
%   scriptGetFretTraces.m
%
% Authors: 
%   - P.G. Jun 2018
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% Get Directories

clear
dirArray = getDirs('Highlight and OPEN Pairs of #Ch1 and #Ch2 Folders (HIT CANCEL TO END SELECTION)');
nExperiments = length(dirArray(1,:));

if nExperiments==0,return,end

%% Parameters
crtValue = 0.00; %0.08

%% Apply FRET Filter

% Concatenate data arrays
for j=1:nExperiments

    % Set current Directopry
    basepath = dirArray{2,j};
    dataSet  = dirArray{1,j};
    currentDir = [basepath filesep dataSet filesep];
    
    disp(dataSet);

    % Load Ch1 Data (Acceptor) 
    if contains(dataSet,'Ch1')
        trackingAcc  = load([currentDir 'Tracking.mat']);
        accFile      = [currentDir 'singleTracesAcc_phc.traces'];
        fretData.Ch1 = loadTracesCell(accFile);

    % Load Ch2 Data (Donor)    
    elseif contains(dataSet,'Ch2')     
        donFile = [currentDir 'singleTracesDon_phc.traces'];
        fretData.Ch2 = loadTracesCell(donFile);
        
        %------------------------------------------------------------------
        %--- Crosstalk correction
        %------------------------------------------------------------------
        if ~crtValue==0
            disp('Crosstalk Correction ...')
            fretData = crosstalkCorr(crtValue, fretData); 
        end
        
        %------------------------------------------------------------------
        %--- Idealize Traces using total Intensity
        %------------------------------------------------------------------
        disp('Idealizing Total Intensity ...');
        [idlTotal,total] = idlTotalIntensity(fretData);
        fretData.Fret.idlTotal = idlTotal;
        fretData.Fret.total    = total;
        
        %------------------------------------------------------------------
        %--- Filter Traces
        %------------------------------------------------------------------
        disp('Filter Traces ...');
        [traceStat,pstFltVal] = pstFilter(fretData);
        fretData.Fret.traceStat = traceStat;
        % Display Filter Results
        nTraces   = size(fretData.Ch1.x,1);
        idxFilter = [traceStat.filterVal];
        mefFilter = [traceStat.fltNotME];
        nMefRej   = numel(find(mefFilter==0));
        nSelTraces = sum(idxFilter,2);
        
        disp(['nSel: ' num2str(nSelTraces) ' of (' num2str(nTraces) ')' ]);  
        disp(['nTraces rej by MEF : ' num2str(nMefRej)]);  
        
        %------------------------------------------------------------------
        %--- Apply Filter and save fretTraces / Tracking.mat
        %------------------------------------------------------------------
        
        disp('Saving ...\fretTracesPst.mat');
        fretTraces     = applyFilterToData(fretData,idxFilter);
        fretTraces.Fret.pstFltVal = pstFltVal;
        save([currentDir filesep 'fretTracesPst.mat'],...
            '-mat','fretTraces','-v7.3');
        
        disp('Saving ...\TrackingPst.mat');
        TrackingPst = applyFilterToTrackingMatrix(trackingAcc,idxFilter);
        if isfield(fretTraces.Fret,'traceStatSel')
            ids = {fretTraces.Fret.traceStatSel.ids}.';
            [TrackingPst.tracksFinal.ids] = ids{:};
        else
            TrackingPst.tracksFinal.ids=[];
            TrackingPst.tracksFinal.tracksFeatIndxCG=[];
            TrackingPst.tracksFinal.tracksCoordAmpCG=[];
            TrackingPst.tracksFinal.seqOfEvents=[];
        end
        
        tracksFinal = orderfields(TrackingPst.tracksFinal,...
            {'ids','tracksFeatIndxCG','tracksCoordAmpCG','seqOfEvents'});
        TrackingPst.tracksFinal = tracksFinal; 
        
        currentDir = [dirArray{2,j-1} filesep dirArray{1,j-1} filesep];
        save([currentDir filesep 'TrackingPst.mat'],...
            '-struct','TrackingPst','-v7.3');    
       
    end 
end



%% =======================================================================%
%                                                                         %
%                        Local Callback Functions                         %
%                                                                         %
%=========================================================================%

%-------------------------------------------------------------------------%
%                                                                         %
% Filter Fret Traces                                                      %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

function [traceStat, fltPrm]=pstFilter(fretData)
    [nTraces,~] = size(fretData.Ch1.x);
    traceStat=struct('ids',        [],...
                     'filterVal',  [],...
                     'meanXSnrAcc',[],...
                     'flt1',       [],...
                     'meanXCrt',   [],...
                     'flt2',       [],...
                     'meanXIntAcc',[],...
                     'flt3',       [],...
                     'meanXIntDon',[],...
                     'flt4',       [],...
                     'meanTotal',  [],...
                     'flt5',       [],...
                     'meanFret',   [],...
                     'flt6',       [],...
                     'lt',         [],...
                     'flt7',       [],...
                     'MESegm',     [],...
                     'MEDist',     [],...
                     'MESortTimes',[],...
                     'MEOffTimes', [],...
                     'fltNotME',   [],...
                     'fretInRange',[],...
                     'totalDwell',[],...
                     'flt9',[]);
                 
    % Filter Limits             
    %--- Signal2BgdRatio (1)
    fltPrm.snrFilterOn = false;
    fltPrm.minSnrVal   = 1.4;
    fltPrm.maxSnrVal   = 2.2;
    
    %--- Crt-Filter (2)
    fltPrm.crtFilterOn = false;
    fltPrm.minCrtVal   = 0.3;
    fltPrm.maxCrtVal   = 10;
    
    %--- Acc-Intensity (3)
    fltPrm.intAccFilterOn = false;
    fltPrm.minIntAcc   = 140;
    fltPrm.maxIntAcc   = inf;
    
    %--- Don-Intensity (4)1
    fltPrm.intDonFilterOn = false;
    fltPrm.minIntDon   = 160;
    fltPrm.maxIntDon   = 800;
    
    %--- Total Intensity Filter (5)
    fltPrm.totalFilterOn = false;
    fltPrm.minTotalVal = 450;            
    fltPrm.maxTotalVal = inf;
    
    %--- Average FRET Filter (6)
    fltPrm.avgFretFilterOn = false;
    fltPrm.minAvgFret  = 0.1;
    fltPrm.maxAvgFret  = 1.0;
    
    %--- Lt-Filter (filter value in frames)(7)
    fltPrm.ltFilterOn = true;
    fltPrm.minLtVal    = 20; 
    
    %--- ME-Filter (8)
    fltPrm.meFilterOn = false;
    fltPrm.trackDist   = 1.5;
    fltPrm.nEvents     = 25;
    fltPrm.offTimeLim  = 0.75;
    
    %--- FRET Range Filter (9)
    fltPrm.fretRangeFilterOn = false;
    fltPrm.absFret   = 0.8; % Fret Value
    fltPrm.fretRange = 0.1; % 0.8 +/- 0.1
    fltPrm.totalDwt  = 15;  % Frames
    
    h=waitbar(0,'Post-Filter Traces...');
    for i=1:nTraces
        %--- Update Waitbar
        waitbar(i / nTraces)
        
        %------------------------------------------------------------------
        %--- Get IDs
        %------------------------------------------------------------------
        traceStat(i).ids = fretData.Ch1.traceMetadata(i).ids;
        
        %------------------------------------------------------------------
        %--- Snr-Filter: locAcc/locBgd Ratio (1)
        %------------------------------------------------------------------
        meanXSnrAcc = fretData.Ch1.traceMetadata(i).meanXSnrAcc;
        if fltPrm.snrFilterOn==true
              idxSnr = meanXSnrAcc > fltPrm.minSnrVal &...
                       meanXSnrAcc < fltPrm.maxSnrVal;
        else; idxSnr = true; 
        end
        traceStat(i).meanXSnrAcc = meanXSnrAcc;
        traceStat(i).flt1        = idxSnr;
        
        %------------------------------------------------------------------
        %--- D/A CrossTalk Filter (2)
        %------------------------------------------------------------------
        meanXCrt = fretData.Ch1.traceMetadata(i).meanXCrt;
        if fltPrm.crtFilterOn==true
              idxCrt = meanXCrt > fltPrm.minCrtVal &...
                       meanXCrt < fltPrm.maxCrtVal;
        else; idxCrt = true;
        end 
        traceStat(i).meanXCrt = meanXCrt;
        traceStat(i).flt2     = idxCrt;
        
        %------------------------------------------------------------------
        %--- Acceptor Intensity Filter (3)
        %------------------------------------------------------------------
        meanXIntAcc = fretData.Ch1.traceMetadata(i).meanXIntAcc;
        if fltPrm.intAccFilterOn==true
              idxIntAcc = meanXIntAcc > fltPrm.minIntAcc &...
                          meanXIntAcc < fltPrm.maxIntAcc;
        else; idxIntAcc = true;
        end 
        traceStat(i).meanXIntAcc = meanXIntAcc;
        traceStat(i).flt3        = idxIntAcc;
      
        %------------------------------------------------------------------
        %--- Donor Intensity Filter (4)
        %------------------------------------------------------------------
        meanXIntDon = fretData.Ch1.traceMetadata(i).meanXIntDon;
        if fltPrm.intDonFilterOn==true
              idxIntDon = meanXIntDon > fltPrm.minIntDon &...
                          meanXIntDon < fltPrm.maxIntDon;
        else; idxIntDon = true;
        end 
        traceStat(i).meanXIntDon = meanXIntDon;
        traceStat(i).flt4        = idxIntDon;
        
        %------------------------------------------------------------------
        %--- Total Intensity Filter (5)
        %------------------------------------------------------------------
        % select only the part that is idealized to state 1 over the range 
        % of the tracked donor 
        idl   = logical(fretData.Fret.idlTotal(i,:));
        total = fretData.Fret.total(i,:);
        meanTotal = nanmean(total(idl));
        if fltPrm.totalFilterOn==true
              idxTotal = meanTotal > fltPrm.minTotalVal &...
                         meanTotal < fltPrm.maxTotalVal;
        else; idxTotal = true;
        end
        traceStat(i).meanTotal = meanTotal;
        traceStat(i).flt5      = idxTotal;
   
        %------------------------------------------------------------------
        %--- Average FRET (6)
        %------------------------------------------------------------------
        acceptor = fretData.Ch1.int(i,:);
        fret     = acceptor./total;
        startOfTraceAcc = fretData.Ch1.traceMetadata(i).startOfTrace;
        endOfTraceAcc   = fretData.Ch1.traceMetadata(i).endOfTrace;
        avgFret = nanmean(fret(startOfTraceAcc:endOfTraceAcc));
        if fltPrm.avgFretFilterOn==true
              idxAvgFret = avgFret > fltPrm.minAvgFret &...
                           avgFret < fltPrm.maxAvgFret;
        else; idxAvgFret = true;
        end
        traceStat(i).meanFret = avgFret;
        traceStat(i).flt6     = idxAvgFret;
       
        %------------------------------------------------------------------
        %--- Lifetime Filter (7)
        %------------------------------------------------------------------
        lt  = fretData.Ch1.traceMetadata(i).traceLen;
        if fltPrm.ltFilterOn==true
              idxFretLt = lt >= fltPrm.minLtVal;
        else; idxFretLt = true;
        end
        traceStat(i).lt   = lt;
        traceStat(i).flt7 = idxFretLt;
        
        %------------------------------------------------------------------
        %--- fretRange Filter (9)
        %------------------------------------------------------------------
        acceptor = fretData.Ch1.int(i,:);
        fret     = acceptor./total;
        startOfTraceAcc = fretData.Ch1.traceMetadata(i).startOfTrace;
        endOfTraceAcc   = fretData.Ch1.traceMetadata(i).endOfTrace;
        fret = fret(startOfTraceAcc:endOfTraceAcc);
        if fltPrm.fretRangeFilterOn==true
              idxRange  = fret >= fltPrm.absFret - fltPrm.fretRange &...
                          fret <= fltPrm.absFret + fltPrm.fretRange;
             idxRangeFret = sum(idxRange) >= fltPrm.totalDwt; 
             
            traceStat(i).fretInRange = nanmean(fret(idxRange));
            traceStat(i).totalDwell = sum(idxRange);
                         
        else
            idxRangeFret = true;
            traceStat(i).fretInRange = [];
            traceStat(i).totalDwell = [];
        end
        
        
        traceStat(i).flt9       = idxRangeFret;
        
        %------------------------------------------------------------------
        %--- Apply all Filters  
        %------------------------------------------------------------------
        filterVal = idxSnr     &...
                    idxCrt     &...
                    idxIntAcc  &...
                    idxIntDon  &...
                    idxTotal   &...
                    idxAvgFret &...
                    idxFretLt  &...
                    idxRangeFret;
        
        %------------------------------------------------------------------
        %--- ME-Filter: To save time check for Multiple Events only if all 
        %--- other filters are logical true (8)
        %------------------------------------------------------------------
        if filterVal==true
            rawTraces = fretData.Ch1;
            segments  = getSegments(rawTraces,i);
            valFilter = repeatingEventsFilter(segments,...
                                              fltPrm.nEvents,...
                                              fltPrm.trackDist,...
                                              fltPrm.offTimeLim);
            if fltPrm.meFilterOn==true
                  idxMEFAcc  = valFilter.idxFlt;
            else; idxMEFAcc  = true;
            end
            traceStat(i).fltNotME    = idxMEFAcc;
            traceStat(i).MESegm      = segments;
            traceStat(i).MEDist      = valFilter.dist;
            traceStat(i).MESortTimes = valFilter.sortedTimes;
            traceStat(i).MEOffTimes  = valFilter.offTime;
            filterVal = filterVal & idxMEFAcc;
        else
            traceStat(i).fltNotME    = NaN;
            traceStat(i).MESegm      = NaN;
            traceStat(i).MEDist      = NaN;
            traceStat(i).MESortTimes = NaN;
            traceStat(i).MEOffTimes  = NaN;
        end
        
        %------------------------------------------------------------------
        % Save Filter Value
        %------------------------------------------------------------------
        traceStat(i).filterVal = filterVal;
                                 
                        
    end
    close(h);
end


%-------------------------------------------------------------------------%
%                                                                         %
% Corosstalk Correction                                                   %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

function fretData=crosstalkCorr(crtValue, fretData)

    acc = fretData.Ch1.int;
    don = fretData.Ch2.int;
    
    accCorr = acc - crtValue*don;
    fretData.Ch1.int = accCorr;  
end


%-------------------------------------------------------------------------%
%                                                                         %
% Idealize Total Intensity                                                %
%                                                                         %
% Input: N  = total number of experiments                                 %
%        No = number of current experiment                                %
%        fretdata with fields Ch1 and Ch2                                 %
%                                                                         %
%-------------------------------------------------------------------------%
function [idlData,total]=idlTotalIntensity(fretdata)

	% Check the dimensions
    try
        assert(size(fretdata.Ch1.x,1)==size(fretdata.Ch2.x,1))
    catch
        warning('Acceptor and Donor file have different sizes');
    end
    
    % Loop through traces
    [nTraces,nFrames]=size(fretdata.Ch1.x);
    idlData = zeros(nTraces,nFrames);
    total   = zeros(nTraces,nFrames);
    
    h = waitbar(0,'Idealize Total Intensity...');
    
    for i=1:nTraces
        
        % Progress Bar
        waitbar(i / nTraces)
        
        % Use donorRange to define start / end for total Intensity
        startOfTrace =  fretdata.Ch2.traceMetadata(i).startOfTrace;
        endOfTrace   =  fretdata.Ch2.traceMetadata(i).endOfTrace;
        lenBaseline  =  fretdata.Ch2.traceMetadata(i).lenBaseline;
        
        % Idealize total intensity using skmTotal(from Spartan / Daniel Terry)
        donor    = fretdata.Ch2.int(i,startOfTrace:endOfTrace+lenBaseline);
        acceptor = fretdata.Ch1.int(i,startOfTrace:endOfTrace+lenBaseline);   
        tmpTotal = donor + acceptor;
        total(i,startOfTrace:endOfTrace+lenBaseline)=tmpTotal;
        
        % Check if total has NaN Values and idealize data
        if any(isnan(tmpTotal))
            totalNearest = fillmissing(tmpTotal,'linear','EndValues','nearest');
            idlData(i,startOfTrace:endOfTrace+lenBaseline) = skmTotal(totalNearest);
        else
            idlData(i,startOfTrace:endOfTrace+lenBaseline) = skmTotal(tmpTotal);
        end
        
        % Set idl state to zero where the donor is not tracked
        idlData(i,endOfTrace+1:endOfTrace+lenBaseline) = 0;
    end
    close(h)
end


%-------------------------------------------------------------------------%
%                                                                         %
% Apply Filter To Fret Data                                               %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%
function fltTraces=applyFilterToData(fretData,idxFilter)

    fltTraces=[];
    
    % Channel Names (Acc: Ch1, Don: Ch2)  
    chName = {'Ch1','Ch2'};
    
    % Metadata Field Names
    metaFN = fieldnames(fretData.Ch1.traceMetadata);
    traceStatFN = fieldnames(fretData.Fret.traceStat);
         
    % Update Ch1(Acc) and Ch2(Don)
    nTraces=size(fretData.Ch1.x);   
    for c=1:2
        
        % Update data
        Ch=chName{c};
        fltTraces.(Ch).time         = fretData.(Ch).time';
        fltTraces.(Ch).x            = fretData.(Ch).x(idxFilter,:);
        fltTraces.(Ch).y            = fretData.(Ch).y(idxFilter,:);
        if isfield(fretData.(Ch),'xCorr')
            fltTraces.(Ch).xCorr    = fretData.(Ch).xCorr(idxFilter,:);
            fltTraces.(Ch).yCorr    = fretData.(Ch).yCorr(idxFilter,:);
        end
        fltTraces.(Ch).int          = fretData.(Ch).int(idxFilter,:);
        fltTraces.(Ch).snr          = fretData.(Ch).snr(idxFilter,:);

        % Update Metadata for each trace 
        j=0;
        for i=1:nTraces
            
            % Apply Filter
            if idxFilter(i)==0; continue;end
            j=j+1;
            
            %Update Metadata
            for k=1:numel(metaFN)
                fn=metaFN{k};
                if isfield(fretData.(Ch).traceMetadata,fn)
                    fltTraces.(Ch).traceMetadata(j).(fn) =...
                    fretData.(Ch).traceMetadata(i).(fn);
                end
            end
            
        end % End Update Metadata 
    end %End Update Ch1,Ch2
    
    % Update FRET fields
    fltTraces.Fret.idlTotal = fretData.Fret.idlTotal(idxFilter,:);
    fltTraces.Fret.total    = fretData.Fret.total(idxFilter,:);
    
    % Update Fret Trace Stat for each trace 
    m=0;n=0; 
    for i=1:nTraces
        % Apply Filter
        if idxFilter(i)==0
            m=m+1;
            %Update traceStat (rejected)
            for k=1:numel(traceStatFN)
                fn=traceStatFN{k};
                if isfield(fretData.Fret.traceStat,fn)
                    fltTraces.Fret.traceStatRej(m).(fn) =...
                    fretData.Fret.traceStat(i).(fn);
                end
            end
            
        else
            n=n+1;
            %Update traceStat (selected)
            for k=1:numel(traceStatFN)
                fn=traceStatFN{k};
                if isfield(fretData.Fret.traceStat,fn)
                    fltTraces.Fret.traceStatSel(n).(fn) =...
                    fretData.Fret.traceStat(i).(fn);
                end
            end
        end

    end % End Update Metadata 
    
    
    
end %End Function


%-------------------------------------------------------------------------%
%                                                                         %
% Apply Filter To Tracking Matrix                                         %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

 

function trackingPst = applyFilterToTrackingMatrix(tracking,idxFilter)
    
    trackingFN  = fieldnames(tracking);
    tracksFinalFN = fieldnames(tracking.tracksFinal);
    
    
    nTraces=size(tracking.tracksFinal,1);
    assert(nTraces==size(idxFilter,2),...
        'The tracking matrix has a different size than the filter');
    j=0;
    
    % Filter tracksFinal
    for i=1:nTraces

        % Apply Filter
        if idxFilter(i)==0; continue;end
        j=j+1;

        %Update Tracking.mat
        for k=1:numel(tracksFinalFN)
            fn= tracksFinalFN{k};
            if isfield(tracking.tracksFinal,fn)
                trackingPst.tracksFinal(j).(fn) =...
                tracking.tracksFinal(i).(fn);
            end
        end  

    end
    
    % Update Tracking.mat
    for k=1:numel(trackingFN)
        fn = trackingFN{k};
        if strcmp(fn,'tracksFinal'),continue,end
        trackingPst.(fn) = tracking.(fn);
    end
   
end



