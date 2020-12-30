function [transDiffAnalysisRes,errFlag] = basicTransientDiffusionAnalysisv1(tracks,...
    probDim,alphaValues,alphaAsym,minDuration,plotRes,confRadMin,peakAlpha)
%BASICTRANSIENTDIFFUSIONANALYSISV1 detects potential diffusion segments of a track and performs MSS analysis on these segments
%
%SYNOPSIS [transDiffAnalysisRes,errFlag] = basicTransientDiffusionAnalysis(tracks,...
%     probDim,alphaValues,minDuration,plotRes,confRadMin,checkAsym)
%
%INPUT  tracks      : -- EITHER --
%                     Output of trackWithGapClosing (matrix),
%                     -- OR --
%                     Output of trackCloseGapsKalman (structure, possibly
%                     with merges/splits).
%
%       probDim     : Problem dimensionality.
%                     Optional. Default: 2.
%       alphaValues : Alpha-value for classification. Can take the values
%                     0.2, 0.1, 0.05 and 0.01. One can enter one value, in
%                     which case it will be used for both confined and
%                     directed, or two values, where the first will be used
%                     for confined and the second for directed.
%                     Optional. Default: 0.05 for both.
%       minDuration : Minimum amount of rolling windows a track can have to
%                     be segmented 
%                     Optional. Default: 5.
%       plotRes     : 1 to plot results, 0 otherwise.
%                     Optional. Default: 0.
%                     Results can be plotted only if problem is 2D.
%                     color-coding:
%                     *brown: immobile
%                     *blue: confined diffusion.
%                     *cyan: normal diffusion.
%                     *magenta: super diffusion.
%                     *black: unclassified.
%       confRadMin  : 1 to calculate the confinement radius of confined
%                     particles using the minimum positional standard
%                     deviation, 0 to calculate it using the mean
%                     positional standard deviation.
%                     Optional. Default: 0.
%
%
%       peakAlpha   : confidence level for choosing peaks when initially
%       segmenting track. Default : 95
%
%OUTPUT transDiffAnalysisRes : And array (size = number of tracks) with the
%                     field ".segmentClass", which contains the fields:
%           .momentScalingSpectrum: (Number of classification
%                     subparts)-by-(20+probDim) matrix, where each row
%                     contains the information for one classification
%                     subpart, and the columns store the following:
%                     (1) Start frame of subpart.
%                     (2) End frame of subpart.
%                     (3) Classification of subpart: 1 = confined, 2 =
%                         free, 3 = directed.
%                     (4) MSS slope resulting in classification.
%                     (5-11) Generalized diffusion coefficients for moment
%                            orders 0-6.
%                     (12-18) Scaling power for moment orders 0-6.
%                     (19) Normal diffusion coefficient (from the MSD).
%                     (20) Confinement radius or localization precision, if subpart is classified as
%                          confined or immobile (NaN otherwise).
%                     (21/22/23) Center of subpart, if subpart is
%                                classified as confined (NaN otherwise).
%           .momentScalingSpectrum1D: same as above but in 1D for asym
%           tracks, additionally includes preferred direction value
%           .asymmetry, .asymmetryAfterMSS: NOT IMPLEMENTED RIGHT NOW.
%
%       errFlag         : 0 if executed normally, 1 otherwise.
%
%REMARKS While tracks do not have to be linear in order to be asymmetric,
%the last analysis step assumes that tracks are linear.
%
%Khuloud Jaqaman, March 2008
%Updated Tony Vega, July 2016

%% Output
checkAsym = 0; %Needs work so not imposed currently
transDiffAnalysisRes = [];
errFlag = 0;
%% Input

%check whether tracks were input
if nargin < 1
    disp('--basicTransientDiffusionAnalysisv1: Please input at least the tracks to be analyzed!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

if nargin < 3 || isempty(alphaValues)
    alphaValues = -0.05;
end

if nargin < 4 || isempty(alphaAsym)
%     alphaValues = 0.05; 
      alphaAsym = 0.05;
end

if nargin < 5 || isempty(minDuration)
    minDuration = 5;
end

if nargin < 6 || isempty(plotRes)
    plotRes = 1;
elseif plotRes == 1 && probDim ~= 2
    disp('--basicTransientDiffusionAnalysisv1: Cannot plot tracks if problem is not 2D!');
    plotRes = 0;
end

if nargin < 7 || isempty(confRadMin)
    confRadMin = 0;
end


if nargin < 8 || isempty(peakAlpha)
    peakAlpha = 95;
end

if errFlag
    disp('--trackTransientDiffusionAnalysis1: Please fix input variables');
    return
end
    p = mfilename('fullpath');
    load(strcat(p(1:end-33),'positionConfidenceCI.mat'))
    load(strcat(p(1:end-33),'positionConfidenceFC.mat'))

%define window sizes
windowAsym = 5;
windowMSS =11;
windowMSSMin = 20;
halfWindowAsym = (windowAsym - 1) / 2;
halfWindowMSS = (windowMSS - 1) / 2;

%specify MSS analysis moment orders
momentOrders = 0 : 6;

%% Track extraction for analysis

%store input tracks in a new variable
tracksInput = tracks;

%extract segments for analysis if tracks were input as a structure that
%might contain merges and splits
%the point is to reduce compound tracks that contain merges and splits into
%simple separate tracks
%thus this step is not necessary if the tracks were input as a matrix,
%which by definition does not contain unresolved compound tracks.
if isstruct(tracks)

    %get number of input tracks from structure
    numInputTracks = length(tracksInput);

    clear tracks

    [tracks,dummy,compTrackStartRow,numSegments] = ...
                convStruct2MatIgnoreMS(tracksInput);


else

    %get number of input tracks from matrix
    numInputTracks = size(tracksInput,1);

    %indicate rows where tracks start (trivial in this case)
    compTrackStartRow = (1 : numInputTracks)';

    %indicate number of segments in each track (1 for all tracks)
    numSegments = ones(numInputTracks,1);

end

%get number of track segments to be analyzed
numTrackSegments = size(tracks,1);

%get track segment start times, end times and life times
trackSEL = getTrackSEL(tracks);

%find track segments that are long enough for analysis
if checkAsym
    indx4analysis = find(trackSEL(:,3) >= windowAsym);
else
    indx4analysis = find(trackSEL(:,3) >= windowMSSMin);
end


indxNot4analysis = setdiff((1:numTrackSegments)',indx4analysis);

%% Rolling window classification

%reserve memory %This will only have certain traits, probably none of these
trackSegmentClassRes = repmat(struct('asymmetry',NaN(1,3),...
    'momentScalingSpectrum',NaN(1,20+probDim),...
    'momentScalingSpectrum1D',NaN(1,20+probDim),...
    'asymmetryAfterMSS',NaN(1,3)),...
    numTrackSegments,1);

%go over all analyzable track segments
gaussDeriv = cell(numTrackSegments,1);
fullDist = [];
for iTrack = indx4analysis'

    %get track segment start, end and life times
    trackSELCurrent = trackSEL(iTrack,:);

    %% Maximum displacement analysis to divide track segment into parts with potentially different diffusion behavior

    %find the length of each part
    trackPartLength = trackSELCurrent(1,3);

    %get starting point of this part
    trackPartStart = trackSELCurrent(1,1);

    %get number of MSS analysis rolling windows in this part
    numRollWindows = trackPartLength - windowMSS + 1;

    %if number of rolling windows is larger than the minimum
    %required duration of a classification, proceed with rolling
    %window analysis
    if numRollWindows > minDuration(1)

        %initialize max displacement vector
        maxDisplacement = NaN(numRollWindows,1); 
        %go over all windows
        for iWindow = 1 : numRollWindows

            %get window start and end points
            startPoint = trackPartStart + iWindow - 1;
            endPoint = startPoint + windowMSS - 1;

            test = tracks(iTrack,8*(startPoint-1)+1:8*endPoint);
            xTest = test(1:8:end);
            yTest = test(2:8:end);
            X =[xTest',yTest'];
            D = pdist(X,'euclidean');
            maxDisplacement(iWindow) = max(D);



        end %(for iWindow = 1 : numRollWindows)
        
        %% Derivative of Maximum Displacement
        h =1;
        y = maxDisplacement;
        [out] = filterGauss1D(y, 2, 'symmetric'); 
        der = diff(out)/h;
%                 der = gradient(out);
        absDer = abs(der);
        normFactor = ((maxDisplacement(2:end)+maxDisplacement(1:end-1))./2);
        normFactor(normFactor==0) = NaN;
        normDer = absDer./normFactor;
        gaussDeriv{iTrack} = normDer;

    else
        gaussDeriv{iTrack} = [];
    end
end


for iTrack = indx4analysis'
    %% Separate tracks that change and those that don't
    trackFull = gaussDeriv{iTrack};
    if ~isempty(trackFull)
        
        level = prctile(trackFull,peakAlpha);

        %% Asymmetry detection first
        %1. Check track to see if any sections have asymmetry
        %this classification scheme is taken from Huet et al (BJ 2006)
        %it classifies tracks as asymmetric or not, based on the scatter of
        %positions along them
        
        if checkAsym       
            %Divide track into segments using asym minimum
            [segPointsA,peakThresh] = findCloseSegments(trackFull,level,halfWindowMSS,windowAsym,windowMSSMin);
            [segPointsA] = findDiffSegments(segPointsA,windowAsym,peakThresh);% using halfWindowMSS to look at same part of track analyzed above
            n = 1:length(segPointsA)-1;
            difference = segPointsA(n+1)-segPointsA(n);
            trackSELCurrent = trackSEL(iTrack,:);
            partClassAsym = NaN(length(n),3);
            %go over all of these segments and determine if they are
            %asymmetric
            for k = 1:length(segPointsA)-1 
            startPoint = trackSELCurrent(1) +segPointsA(k)-1;
            endPoint  = startPoint+ difference(k)-1;

                %get the particle positions along the track
                coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                %coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';%Not implemented yet
                coordXY = [coordX coordY];

                %determine whether the track is sufficiently asymmetric
                [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);

                %classify track as ...
                %1 = linear, if the asymmetry parameter is larger than the threshold
                %0 = not linear, if the asymmetry parameter is smaller than the
                %threshold
                %otherwise, keep track classification as undetermined
                partClassAsym(k,1) = startPoint;
                partClassAsym(k,2) = endPoint;
                partClassAsym(k,3) = asymFlag;
            end


        %find indices of all tracks classified as asymmetric
        indxAsym = find(partClassAsym(:,3) == 1);
        else
            %Other values which are not defined, given value of NaN or
            %similar
                indxAsym = 0;
                [segPointsA,peakThresh] = findCloseSegments(trackFull,level,halfWindowMSS,6,windowMSSMin);
                [segPointsA,peakThresh] = findDiffSegments(segPointsA,windowMSSMin,peakThresh);% halfWindowMSS,windowMSSMin
        end
        
                        %% New filtering section
                        % Does not work when a section has all NaNs
                       
                if size(segPointsA,1)>=2%<1%
                    pointsTmp = segPointsA;
                    list = [];
                    sel = trackSEL(iTrack,1);
                    for seg = 1:size(segPointsA,1)-2
                            %Code to compare these segments
                            b1 = (pointsTmp(seg))+sel-1;
                            b2 =(pointsTmp(seg+1))+sel-1;
                            e1 = (pointsTmp(seg+1)-1)+sel-1;
                            e2 =(pointsTmp(seg+2)-1)+sel-1;

                            test1 = tracks(iTrack,8*(b1-1)+1:8*e1);
                            xTest1 = test1(1:8:end);
                            yTest1 = test1(2:8:end);
                            X1 =[xTest1',yTest1'];
                            D1n = X1(2:end,:)-X1(1:end-1,:);
                            magVector1 = sqrt(((D1n(:,2)).^2)+((D1n(:,1)).^2));
                            D1 = pdist(X1,'euclidean');

                            test2 = tracks(iTrack,8*(b2-1)+1:8*e2);
                            xTest2 = test2(1:8:end);
                            yTest2 = test2(2:8:end);
                            X2 =[xTest2',yTest2'];
                            D2n = X2(2:end,:)-X2(1:end-1,:);
                            magVector2 = sqrt(((D2n(:,2)).^2)+((D2n(:,1)).^2));
                            D2 = pdist(X2,'euclidean');
                            deg1 = atand(X1(:,2)./X1(:,1));
                            deg2 = atand(X2(:,2)./X2(:,1));
                            %Try to normalize things
                            bins = linspace(0,max(max(D2),max(D1)),20);
%                             [newD1,~] = histcounts(D1,bins);
%                             [newD2,~] = histcounts(D2,bins);
%                             normD1 = newD1/sum(newD1);
%                             normD2 = newD2/sum(newD2);
                            D1S = datasample(D1, 500);
                            D2S = datasample(D2, 500);
                            try
                                [~,p] = kstest2(D1,D2);
                            catch
%                                 warning('Segment Filter: Not enough data in one segment')
                                p = 0;
                            end
                            if p > 0.05
                                list = [list;seg+1];
                            end
                    end
                    if ~isempty(list)
                        segPointsA(list) = [];
                        peakThresh(list-1) = [];
                                if checkAsym
                                    %go over all of these segments and determine if they are
                                    %asymmetric
                                    n = 1:length(segPointsA)-1;
                                    difference = segPointsA(n+1)-segPointsA(n);
                                    for k = 1:length(segPointsA)-1 
                                    startPoint = trackSELCurrent(1) +segPointsA(k)-1;
                                    endPoint  = startPoint+ difference(k)-1;

                                        %get the particle positions along the track
                                        coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                        coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                                        %coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';%Not implemented yet
                                        coordXY = [coordX coordY];

                                        %determine whether the track is sufficiently asymmetric
                                        [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);

                                        %classify track as ...
                                        %1 = linear, if the asymmetry parameter is larger than the threshold
                                        %0 = not linear, if the asymmetry parameter is smaller than the
                                        %threshold
                                        %otherwise, keep track classification as undetermined
                                        partClassAsym(k,1) = startPoint;
                                        partClassAsym(k,2) = endPoint;
                                        partClassAsym(k,3) = asymFlag;
                                    end


                                %find indices of all tracks classified as asymmetric
                                indxAsym = find(partClassAsym(:,3) == 1);
                                end
                    end

                end
                %}
                %% End of new filtering
        %% Back to normal
        %Initialize all variables that will be saved later
        n = 1:length(segPointsA)-1;
        difference = segPointsA(n+1)-segPointsA(n);
        pointClassMSS = NaN(length(n),1);
        mssSlope = pointClassMSS;
        normDiffCoef = pointClassMSS;
        confRadTmp = pointClassMSS;
        centerTmp = NaN(length(n),probDim);
        genDiffCoef = NaN(length(n),length(momentOrders));
        scalingPower = NaN(length(n),length(momentOrders));
        partClassMSS = NaN(length(n),3);
        
        partClassMSS1D = NaN(length(n),3);
        mssSlope1D = pointClassMSS;
        genDiffCoef1D = NaN(length(n),length(momentOrders));
        scalingPower1D = NaN(length(n),length(momentOrders));
        normDiffCoef1D = pointClassMSS;
        confRadius1D = NaN(length(n),2);
        prefDir = NaN(length(n),2);
        trackCenter = NaN(length(n),2);
        %% Diffusion Analysis
        %Now go through segments and get diffusion classification.
        trackSELCurrent = trackSEL(iTrack,:);
        
        for k = 1:length(segPointsA)-1 
        startPoint = trackSELCurrent(1) +segPointsA(k)-1;
        if k == length(segPointsA)-1
           endPoint  = startPoint+ difference(k)-1;  
        else
           endPoint  = startPoint+ difference(k);
        end
            if ismember(k,indxAsym)    

                        %get the positions in this track and their standard deviations

                    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter(k,:),confRadius1D(k,:),prefDir(k,:)] = asymmetricDiffusion(startPoint,endPoint,alphaAsym,probDim,tracks,iTrack);

                    %since not all track segments are linear, put analysis results in their
                    %proper place among all track segment
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;              
                    partClassMSS1D(k,3) = partClassAsym(k,3);
                    partClassMSS(k,1)= startPoint;
                    partClassMSS(k,2)= endPoint;
                    partClassMSS(k,3)= pointClass;
                    mssSlope1D(k) = mssSlopeT;
                    genDiffCoef1D(k,:) = genDiffCoefT;
                    scalingPower1D(k,:) = scalingPowerT;
                    normDiffCoef1D(k) = normDiffCoefT;

            else

                
                        [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                            tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                            probDim,momentOrders,alphaValues(1));


                        if pointClassMSS(k) == 1
                        [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);

                        elseif pointClassMSS(k) == 0
                          [confRadTmp(k)] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                            centerTmp(k,:) = NaN(1,probDim);
                        else
                            confRadTmp(k) = NaN;
                            centerTmp(k,:) = NaN(1,probDim);    
                        end 
                        partClassMSS1D(k,1)= startPoint;
                        partClassMSS1D(k,2)= endPoint;
                        partClassMSS(k,1)= startPoint;
                        partClassMSS(k,2)= endPoint;
                        partClassMSS(k,3)= pointClassMSS(k);
            end
                    
        end
         %% Try to merge unclassified segments
             %Check if there are symmetric unclassified segments
                    partClassTmp = partClassMSS;
                    mssSlopeTmp = mssSlope;
                    unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                    
                    if sum(unclassCheck)>0 && size(partClassMSS,1)>1
                        rnd = find(unclassCheck);
                        check = find(unclassCheck);
                        for un = 1:length(rnd)
                            [partClassMSS,partClassMSS1D,mssSlope] = mergeUnclassSegments(unclassCheck,check(1),partClassMSS,partClassMSS1D,trackFull,halfWindowMSS,mssSlopeTmp);


                            %if parts have been merged, redo diffusion analysis on
                            %these parts. Otherwise continue
                            if size(partClassMSS,1) < size(partClassTmp,1)
        %                         [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
        %                         partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);

                            [partClassMSST,mssSlopeT,normDiffCoefT,confRadTmpT,centerTmpT,genDiffCoefT,scalingPowerT,...
                                partClassMSS1DT,mssSlope1DT,normDiffCoef1DT,confRadius1DT,trackCenterT,genDiffCoef1DT,scalingPower1DT,prefDirT] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
                            %    
                                %% New section to avoid possibly erroneous
                                %classification: If new classification results in
                                %more mobile class, revert to old class
                                changeSeg = find((partClassMSST(:,3)-partClassMSS(:,3)) > 0);
                                if changeSeg
%                                     if size(partClassMSST,1) == length(changeSeg) 
                                        if partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==2 %If confined track became free
                                            oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),2);
                                            
                                        elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==1 %If immobile track became confined
                                            oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),2);
                                        elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==0 %If immobile track became free
                                            oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),2);
                                         %Testing use on switches to lower mobility
                                        elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==1 %If free track became confined
                                            oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),2);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),3);
                                        elseif partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==0 %If confined track became immobile
                                            oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),2);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),3);
                                        elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==0 %If free track became immobile
                                            oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),2);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),3);   
                                        else % We currently don't have confidence for super, so this is just not accepted
                                            newConfidence = 1;
                                            oldConfidence = 0;
                                        end
                                        if newConfidence < oldConfidence %If merge improved confidence, accept.
                                            unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                                            check = find(unclassCheck);
                                            %Accept and make all T's into actual value
                                            partClassMSS = partClassMSST;
                                            partClassTmp = partClassMSS;
                                            mssSlope = mssSlopeT;
                                            mssSlopeTmp = mssSlope;
                                            normDiffCoef = normDiffCoefT;
                                            confRadTmp = confRadTmpT;
                                            centerTmp = centerTmpT;
                                            genDiffCoef = genDiffCoefT;
                                            scalingPower = scalingPowerT;
                                            partClassMSS1D =partClassMSS1DT;
                                            mssSlope1D = mssSlope1DT;
                                            normDiffCoef1D =normDiffCoef1DT;
                                            confRadius1D = confRadius1DT;
                                            trackCenter = trackCenterT;
                                            genDiffCoef1D = genDiffCoef1DT;
                                            scalingPower1D = scalingPower1DT;
                                            prefDir = prefDirT;
                                            clear oldConfidence newConfidence
                                        else %If not, get previous classification and redo classification, CHECK if this is necessary...
                                          partClassMSS = partClassTmp;
                                        [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
                                    partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
                                            clear oldConfidence newConfidence
                                        end
% %                                     else
%                                         %Get segments that will stay new (i.e did not result in positive value) 
%                                         keep = partClassMSS(~ismember(1:size(partClassMSS,1),changeSeg),:);%partClassMSS(find(~(partClassMSST(:,3)-partClassMSS(:,3)) > 0),:);
%                                         %Get segments that will revert (i.e resulted in postive value)
%                                         partClassR = partClassTmp;
%                                         for remove = 1:size(keep,1)
%                                             birth = find(partClassR(:,1)==keep(remove,1));
%                                             death =  find(partClassR(:,2)==keep(remove,2));
%                                             partClassR(birth:death,:) = [];
%                                         end
%                                         %Stich everything together and re-run,
%                                         %technically we have all the information, but
%                                         %this is easiest eay to stich data together
%                                         partClassMSS = sortrows([keep;partClassR]);
%                                         [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
%                                     partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
% %                                     end        

                                else
                                    unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                                    check = find(unclassCheck);
                                    %Accept and make all T's into actual value
                                    partClassMSS = partClassMSST;
                                    partClassTmp = partClassMSS;
                                    mssSlope = mssSlopeT;
                                    mssSlopeTmp = mssSlope;
                                    normDiffCoef = normDiffCoefT;
                                    confRadTmp = confRadTmpT;
                                    centerTmp = centerTmpT;
                                    genDiffCoef = genDiffCoefT;
                                    scalingPower = scalingPowerT;
                                    partClassMSS1D =partClassMSS1DT;
                                    mssSlope1D = mssSlope1DT;
                                    normDiffCoef1D =normDiffCoef1DT;
                                    confRadius1D = confRadius1DT;
                                    trackCenter = trackCenterT;
                                    genDiffCoef1D = genDiffCoef1DT;
                                    scalingPower1D = scalingPower1DT;
                                    prefDir = prefDirT;
                                end
                                %}
                            end
                        end
                    end
        
         %% Reclassify if segments next to each other are the same
            %merge subparts that now have the same classification
            checkSeg = partClassMSS(1:end-1,3)-partClassMSS(2:end,3);
            indicateRev = 0;
            while ~isempty(find(checkSeg==0)) && indicateRev < 2
                            mssSlopeTmp = mssSlope;
                            clear mssSlope
%                             mssSlope(find(checkSeg==0):find(checkSeg==0)+1) = mssSlope(max(find(checkSeg==0):find(checkSeg==0)+1));
%                             mssSlope = unique(mssSlope); %Could be disastorous 
                            partClassTmp = partClassMSS;
                            lifeTimeMSS = (partClassTmp(:,2)-partClassTmp(:,1))+1;
                            pointClassMSSAsymT = partClassMSS1D(:,3);
                            pointClassMSS = partClassMSS(:,3);
                    for iSubpart = 1 : length(pointClassMSS)-1
                        iSubpartPlus1 = iSubpart + 1;
                        while ( (iSubpartPlus1 <= length(pointClassMSS)) && ...
                                ( (pointClassMSS(iSubpart) == pointClassMSS(iSubpartPlus1)) || ...
                                (isnan(pointClassMSS(iSubpart)) && isnan(pointClassMSS(iSubpartPlus1))) ) && ...
                                ( (pointClassMSSAsymT(iSubpart) == pointClassMSSAsymT(iSubpartPlus1)) || ... %same asym
                                (isnan(pointClassMSSAsymT(iSubpart)) && isnan(pointClassMSSAsymT(iSubpartPlus1))) ) )
                            
                            partClassMSS1D(iSubpart,2) = partClassMSS1D(iSubpartPlus1,2);
                            partClassMSS1D(iSubpartPlus1,1) = partClassMSS1D(iSubpart,1);
                            
                            partClassMSS(iSubpart,2) = partClassMSS(iSubpartPlus1,2);
                            partClassMSS(iSubpartPlus1,1) = partClassMSS(iSubpart,1);
                            iSubpartPlus1 = iSubpartPlus1 + 1;
                        end
                    end
                    [~,uniqueParts] = unique(partClassMSS(:,1));
                    for ms = 1:length(uniqueParts)
                        mssSlope(ms) = sum(lifeTimeMSS(partClassMSS(:,2) == partClassMSS(uniqueParts(ms),2)).*mssSlopeTmp(partClassMSS(:,2) == partClassMSS(uniqueParts(ms),2)))./sum(lifeTimeMSS(partClassMSS(:,2) == partClassMSS(uniqueParts(ms),2)));
                    end
                    partClassMSS = partClassMSS(uniqueParts,:);
%                     mssSlope = mssSlope(uniqueParts);
                    partClassMSS1D = partClassMSS1D(uniqueParts,:);
                    
                    %if parts have been merged, redo diffusion analysis on
                    %these parts. Otherwise continue
                    if size(partClassMSS,1) < size(partClassTmp,1)
                        [partClassMSST,mssSlopeT,normDiffCoefT,confRadTmpT,centerTmpT,genDiffCoefT,scalingPowerT,...
                        partClassMSS1DT,mssSlope1DT,normDiffCoef1DT,confRadius1DT,trackCenterT,genDiffCoef1DT,scalingPower1DT,prefDirT] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
                    %    
                        %New section to avoid possibly erroneous
                        %classification: If new classification results in
                        %more mobile class, revert to old class
                        changeSegF = find((partClassMSST(:,3)-partClassMSS(:,3)) > 0);
                        badChange = [];
                        badC=1;
                        if changeSegF
                            for cs= 1:length(changeSegF)
                                    changeSeg = changeSegF(cs);
                                    if partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==2 %If confined track became free
%                                         oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                        newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                        oldConfidenceList = newConfidenceList;
                                        oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),3);
                                        newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),2);

                                    elseif partClassMSS(changeSeg,3) == 0 && partClassMSST(changeSeg,3) ==1 %If immobile track became confined
%                                          oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                        newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                        oldConfidenceList = newConfidenceList;
                                        oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),3);
                                        newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),2); 
                                    elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==0 %If immobile track became free
%                                             oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                            newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                            oldConfidenceList = newConfidenceList;
                                            oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),3);
                                            newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),2);
                                     %Testing use on switches to lower mobility
                                    elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==1 %If free track became confined
%                                         oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                        newConfidenceList = positionConfidenceFC{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                        oldConfidenceList = newConfidenceList;
                                        oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),2);
                                        newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),3);
                                    elseif partClassMSS(changeSeg,3) == 1 && partClassMSST(changeSeg,3) ==0 %If confined track became immobile
%                                         oldConfidenceList = positionConfidenceCI{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                        newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                        oldConfidenceList = newConfidenceList;
                                        oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),2);
                                        newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),3);
                                    elseif partClassMSS(changeSeg,3) == 2 && partClassMSST(changeSeg,3) ==0 %If free track became immobile
%                                         oldConfidenceList = positionConfidenceFC{1,min(length(partClassTmp(ismember(mssSlopeTmp,mssSlope),1):partClassTmp(ismember(mssSlopeTmp,mssSlope),2)),500)};
                                        newConfidenceList = positionConfidenceCI{1,min(length(partClassMSS(changeSeg,1):partClassMSS(changeSeg,2)),500)};
                                        oldConfidenceList = newConfidenceList;
                                        oldConfidence = oldConfidenceList(round(oldConfidenceList(:,1),4)==round(mssSlope(changeSeg),4),2);
                                        newConfidence = newConfidenceList(round(newConfidenceList(:,1),4)==round(mssSlopeT(changeSeg),4),3);    
                                    else % We currently don't have confidence for super, so this is just not accepted
                                        newConfidence = 1;
                                        oldConfidence = 0;
                                    end
                                    if newConfidence > oldConfidence
                                        badChange(badC)=changeSegF(cs);
                                        badC = badC +1;
                                    end
                            end
                            %If there was a bad change, revert it back but
                            %keep good changes
                            if ~isempty(badChange)
                                if size(partClassMSST,1) == length(changeSegF)
                                      partClassMSS = partClassTmp;
                                    [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
                                partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
                                clear oldConfidence newConfidence
                                else
                                    %Get segments that will stay new (i.e did not result in positive value)
                                    keep = partClassMSS(~ismember(1:size(partClassMSS,1),badChange),:);
                                    %Get segments that will revert (i.e resulted in postive value)
                                    segments = [uniqueParts(2:end);length(partClassTmp)+1]-uniqueParts;
                                    for segS = 1:size(partClassMSS,1)
                                        storage.Segments(segS) = {uniqueParts(segS):uniqueParts(segS)+(segments(segS)-1)};
                                    end
                                    revertSeg = cell2mat(storage.Segments(badChange));
                                    %Stich everything together and re-run,
                                    %technically we have all the information, but
                                    %this is easiest eay to stich data together
                                    partClassMSS = sortrows([keep;partClassTmp(revertSeg,:)]);
                                        [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
                                partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
                                clear oldConfidence newConfidence
                                end
                            else
                            %------------------------------------------------------------
%                             if newConfidence < oldConfidence %If merge improved confidence, accept.

                                %Accept and make all T's into actual value
                                partClassMSS = partClassMSST;
                                partClassTmp = partClassMSS;
                                mssSlope = mssSlopeT;
                                normDiffCoef = normDiffCoefT;
                                confRadTmp = confRadTmpT;
                                centerTmp = centerTmpT;
                                genDiffCoef = genDiffCoefT;
                                scalingPower = scalingPowerT;
                                partClassMSS1D =partClassMSS1DT;
                                mssSlope1D = mssSlope1DT;
                                normDiffCoef1D =normDiffCoef1DT;
                                confRadius1D = confRadius1DT;
                                trackCenter = trackCenterT;
                                genDiffCoef1D = genDiffCoef1DT;
                                scalingPower1D = scalingPower1DT;
                                prefDir = prefDirT;
                                clear oldConfidence newConfidence
%                             else %If not, get previous classification and redo classification, CHECK if this is necessary...
%                               partClassMSS = partClassTmp;
%                             [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
%                         partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
%                                 clear oldConfidence newConfidence
                            end                            
                            %-------------------------------------------------
%                             if size(partClassMSST,1) == length(changeSeg)
%                                   partClassMSS = partClassTmp;
%                                 [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
%                             partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
%                            
%                             else
%                                 %Get segments that will stay new (i.e did not result in positive value)
%                                 keep = partClassMSS(~ismember(1:size(partClassMSS,1),changeSeg),:);%partClassMSS(find(~(partClassMSST(:,3)-partClassMSS(:,3)) > 0),:);
%                                 %Get segments that will revert (i.e resulted in postive value)
%                                 segments = [uniqueParts(2:end);length(partClassTmp)+1]-uniqueParts;
%                                 for segS = 1:size(partClassMSS,1)
%                                     storage.Segments(segS) = {uniqueParts(segS):uniqueParts(segS)+(segments(segS)-1)};
%                                 end
%                                 revertSeg = cell2mat(storage.Segments(changeSeg));
%                                 %Stich everything together and re-run,
%                                 %technically we have all the information, but
%                                 %this is easiest eay to stich data together
%                                 partClassMSS = sortrows([keep;partClassTmp(revertSeg,:)]);
%                                 [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
%                             partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym);
%                             end        
                            indicateRev = indicateRev+1;
                        else
                            %Accept and make all T's into actual value
                            partClassMSS = partClassMSST;
                            mssSlope = mssSlopeT;
                            normDiffCoef = normDiffCoefT;
                            confRadTmp = confRadTmpT;
                            centerTmp = centerTmpT;
                            genDiffCoef = genDiffCoefT;
                            scalingPower = scalingPowerT;
                            partClassMSS1D =partClassMSS1DT;
                            mssSlope1D = mssSlope1DT;
                            normDiffCoef1D =normDiffCoef1DT;
                            confRadius1D = confRadius1DT;
                            trackCenter = trackCenterT;
                            genDiffCoef1D = genDiffCoef1DT;
                            scalingPower1D = scalingPower1DT;
                            prefDir = prefDirT;
                        end
                        %}
                        
                    else
                        checkSeg = 1;
                        continue
                    end
                    checkSeg = partClassMSS(1:end-1,3)-partClassMSS(2:end,3);
            end
            

                  
    else
        %analyze whole track 
        trackSELCurrent = trackSEL(iTrack,:);
        
        %Fill in empty 1D data
        mssSlope1D = NaN(1,1);
        mssSlope = NaN(1,1);
        genDiffCoef1D = NaN(1,length(momentOrders));
        genDiffCoef = NaN(1,length(momentOrders));
        scalingPower1D = NaN(1,length(momentOrders));
        scalingPower = NaN(1,length(momentOrders));
        normDiffCoef1D = NaN(1,1);
        normDiffCoef = NaN(1,1);
        confRadius1D = NaN(1,2);
        confRadTmp = NaN;
        centerTmp = NaN(1,probDim);
        
        prefDir = NaN(1,2);
        trackCenter = NaN(1,2);
        partClassMSS1D = [trackSELCurrent(:,1:2) NaN];
        
        if checkAsym ==1

            %get the particle positions along the track
            coordX = (tracks(iTrack,1:8:end))';
            coordY = (tracks(iTrack,2:8:end))';
%             coordZ = tracks(iTrack,3:8:end)';
            coordXY = [coordX coordY];

            %determine whether the track is sufficiently asymmetric
            [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);

            if asymFlag ==1

                [pointClass,mssSlope1D,genDiffCoef1D,scalingPower1D,normDiffCoef1D,trackCenter,confRadius1D,prefDir] = asymmetricDiffusion(trackSELCurrent(:,1),trackSELCurrent(:,2),alphaAsym,probDim,tracks,iTrack);

                partClassMSS1D = [trackSELCurrent(:,1:2) asymFlag];
                partClassMSS = [trackSELCurrent(:,1:2) pointClass];

            else
                
                trackSELCurrent = trackSEL(iTrack,:);
                [pointClassMSS,mssSlope,genDiffCoef,scalingPower,normDiffCoef] = trackMSSAnalysis(...
                    tracks(iTrack,:),...
                    probDim,momentOrders,alphaValues(1));
                %Get confinement info if necessary
                if pointClassMSS == 1
                    [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,:),probDim,confRadMin);
                elseif pointClassMSS == 0
                          [confRadTmp] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                            centerTmp = NaN(1,probDim);
                else
                    confRadTmp = NaN;
                    centerTmp = NaN(1,probDim);
                end 
                partClassMSS =[trackSELCurrent(:,1:2) pointClassMSS]; 
                partClassMSS1D = [trackSELCurrent(:,1:2) asymFlag];

            end
        else
                [pointClassMSS,mssSlope,genDiffCoef,scalingPower,normDiffCoef] = trackMSSAnalysis(...
                    tracks(iTrack,:),...
                    probDim,momentOrders,alphaValues(1));
        %Get confinement info if necessary
                if pointClassMSS == 1
                    [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,:),probDim,confRadMin);
                elseif pointClassMSS == 0
                          [confRadTmp] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                            centerTmp = NaN(1,probDim);
                else
                    confRadTmp = NaN;
                    centerTmp = NaN(1,probDim);
                end 
                partClassMSS =[trackSELCurrent(:,1:2) pointClassMSS]; 
                    
        end

    end
    if size(partClassMSS,1) > 1 && sum(partClassMSS(:,3)) == 0 
        partClassMSS = partClassMSS;
    end
    
    %% Store results
    trackClassMSS = [partClassMSS,mssSlope,genDiffCoef,scalingPower, normDiffCoef,confRadTmp,centerTmp];
    trackClassMSS1D = [partClassMSS1D,mssSlope1D,genDiffCoef1D,scalingPower1D, normDiffCoef1D,confRadius1D,trackCenter,prefDir]; %COMPLETE
    trackSegmentClassRes(iTrack).asymmetry = partClassMSS1D;% [1 352 NaN]
    trackSegmentClassRes(iTrack).momentScalingSpectrum = trackClassMSS;
    trackSegmentClassRes(iTrack).momentScalingSpectrum1D = trackClassMSS1D;
%     trackSegmentClassRes(iTrack).asymmetryAfterMSS = trackClassAsym;% [1 352 NaN] Verify whether this is needed
    
end
%% Store trivial nonclassification information for tracks that are not classifiable
if ~isempty(indxNot4analysis')
for iTrack = indxNot4analysis' 
    trackSELCurrent = trackSEL(iTrack,:);
    trackSegmentClassRes(iTrack).asymmetry(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).momentScalingSpectrum(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).momentScalingSpectrum1D(1:2) = trackSELCurrent(1:2);
%     trackSegmentClassRes(iTrack).asymmetryAfterMSS(1:2) = trackSELCurrent(1:2);
end
end
%% save results in output structure

%reserve memory
segmentClass = struct('asymmetry',[],'momentScalingSpectrum',[],'momentScalingSpectrum1D',[],'asymmetryAfterMSS',[]);
transDiffAnalysisRes = repmat(struct('segmentClass',segmentClass),numInputTracks,1);

%go over all input tracks
for iTrack = 1 : numInputTracks

    %go over the segments of each track
    for iSegment = 1 : numSegments(iTrack)

        %store the segment's classification results
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetry = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetry;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).momentScalingSpectrum = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).momentScalingSpectrum;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).momentScalingSpectrum1D = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).momentScalingSpectrum1D;
%         transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetryAfterMSS = ...
%             trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetryAfterMSS;

    end %(for iSegment = 1 : numSegments(iTrack))

end %(for iTrack = 1 : numInputTracks)

%% plotting

%plot results if requested
if plotRes
    plotTracksTransDiffAnalysis2D(tracksInput,transDiffAnalysisRes,[],1,[],[],checkAsym,[]);
end



end

%% Subfunctions
function [partClassMSS,mssSlope,normDiffCoef,confRadTmp,centerTmp,genDiffCoef,scalingPower,...
partClassMSS1D,mssSlope1D,normDiffCoef1D,confRadius1D,trackCenter,genDiffCoef1D,scalingPower1D,prefDir] = reclassifySegments(partClassMSS,probDim, tracks, iTrack,momentOrders,alphaValues,alphaAsym,confRadMin,checkAsym)                 
    numSeg = size(partClassMSS,1);%This is here bc things have been analyzed, make earlier?
    pointClassMSS = NaN(numSeg,1);
%                             partClassAsym = NaN(numSeg,1);%point to this?
    mssSlope = pointClassMSS;
    normDiffCoef = pointClassMSS;
    confRadTmp = pointClassMSS;
    centerTmp = NaN(numSeg,probDim);
    genDiffCoef = NaN(numSeg,length(momentOrders));
    scalingPower = NaN(numSeg,length(momentOrders));

    partClassMSS1D = NaN(size(partClassMSS,1),3);
    mssSlope1D = pointClassMSS;
    genDiffCoef1D = NaN(numSeg,length(momentOrders));
    scalingPower1D = NaN(numSeg,length(momentOrders));
    normDiffCoef1D = pointClassMSS;
    confRadius1D = NaN(numSeg,2);
    prefDir = NaN(numSeg,2);
    trackCenter = NaN(numSeg,2);
    asymFlag = 0;
    for k = 1:numSeg 
    startPoint = partClassMSS(k,1);
    endPoint  = partClassMSS(k,2);
%     alphaAsym = 0.05;%alphaValues(2);
                if checkAsym
                %get the particle positions along the track
                coordX = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                coordY = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
%                                     coordZ = (tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint))';
                coordXY = [coordX coordY];

                %determine whether the track is sufficiently asymmetric

                    [~,asymFlag] = asymDeterm2D3D(coordXY(:,1:probDim),alphaAsym);
                end
                %classify track as ...
                %1 = linear, if the asymmetry parameter is larger than the threshold
                %0 = not linear, if the asymmetry parameter is smaller than the
                %threshold
                %otherwise, keep track classification as undetermined

                if asymFlag ==1

                    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter(k,:),confRadius1D(k,:),prefDir(k,:)] = asymmetricDiffusion(startPoint,endPoint,alphaAsym,probDim,tracks,iTrack);

                    %since not all track segments are linear, put analysis results in their
                    %proper place among all track segment
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;              
                    partClassMSS1D(k,3) = asymFlag;
                    partClassMSS(k,1)= startPoint;
                    partClassMSS(k,2)= endPoint;
                    partClassMSS(k,3)= pointClass;

                    mssSlope1D(k) = mssSlopeT;
                    genDiffCoef1D(k,:) = genDiffCoefT;
                    scalingPower1D(k,:) = scalingPowerT;
                    normDiffCoef1D(k) = normDiffCoefT;
                else

                    [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                        tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                        probDim,momentOrders,alphaValues(1));
                    
                    partClassMSS1D(k,1)= startPoint;
                    partClassMSS1D(k,2)= endPoint;
                    
                    if pointClassMSS(k) == 1  
                    [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                    elseif pointClassMSS(k) == 0
                          [confRadTmp(k)] = estimLocPre(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim);
                            centerTmp(k,:) = NaN(1,probDim);
                    else
                        confRadTmp(k) = NaN;
                        centerTmp(k,:) = NaN(1,probDim);
                    end
                    partClassMSS(k,3) = pointClassMSS(k);
                end
    end

end

function  [partClassMSS,partClassMSS1D,mssSlope] = mergeUnclassSegments(unclassCheck,check,partClassMSS,partClassMSS1D,trackFull,halfWindowMSS,mssSlopeTmp)%unclassCheck
    %cycle through segments
%     check = find(unclassCheck);
%         [~,check] =sort(difference(difference < windowMSSMin)); 
%         while ~isempty(check) && size(partClassMSS,1)>1
mssSlope =mssSlopeTmp;
        if size(partClassMSS,1)>1
        %If short segment is first, connect to succeeding
            if check(1) == 1 
                partClassMSS(1,2) = partClassMSS(2,2);
                partClassMSS1D(1,2) = partClassMSS1D(2,2);
                partClassMSS(1,3) = partClassMSS(2,3);
                partClassMSS1D(1,3) = partClassMSS1D(2,3);
                partClassMSS(2,:) =[];
                partClassMSS1D(2,:) =[];
                mssSlope(1) = mssSlope(2);
                mssSlope(2) = [];
                unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
% %                 check = find(unclassCheck);
% %                 continue
%             end
             %If segment is last, connect to preceding
            elseif check(1) == length(unclassCheck) %changed from end
                partClassMSS(end-1,2) = partClassMSS(end,2);
                partClassMSS1D(end-1,2) = partClassMSS1D(end,2);
                partClassMSS(end,:) =[];
                partClassMSS1D(end,:) =[];
                mssSlope(end) = mssSlope(end-1);
                mssSlope(end-1) = [];
                unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
% %                 check = find(unclassCheck);
% %                 continue
%             end

            %If segment is somewhere in middle
            elseif ~isempty(check)
                %If preceeding and succeeding segments are the same
                %diffusion, merge everything
                if partClassMSS(check(1)+1,3) == partClassMSS(check(1)-1,3) && isnan(partClassMSS1D(check(1)+1,3)) == isnan(partClassMSS1D(check(1)-1,3))
                    partClassMSS(check(1)-1,2) = partClassMSS(check(1)+1,2);
                    partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1)+1,2);
                    partClassMSS(check(1)+1,:) =[];
                    partClassMSS1D(check(1)+1,:) =[];
                    partClassMSS(check(1),:) =[];
                    partClassMSS1D(check(1),:) =[];
                    mssSlope(check(1)-1) = max(mssSlope(check(1):check(1)+1));
                    mssSlope(check(1)+1) =[];
                    mssSlope(check(1)) =[];
                elseif partClassMSS1D(check(1)+1,3) == partClassMSS1D(check(1)-1,3) && isnan(partClassMSS(check(1)+1,3)) == isnan(partClassMSS(check(1)-1,3))
                    partClassMSS(check(1)-1,2) = partClassMSS(check(1)+1,2);
                    partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1)+1,2);
                    partClassMSS(check(1)+1,:) =[];
                    partClassMSS1D(check(1)+1,:) =[];
                    partClassMSS(check(1),:) =[];
                    partClassMSS1D(check(1),:) =[];  
                    %Need to fix mss Slope here
                 %Otherwise score the segments and merge to weaker score
                else
                    if trackFull(partClassMSS(check(1)+1,1)-(partClassMSS(1,1)+halfWindowMSS)) > trackFull(partClassMSS(check(1),1)-(partClassMSS(1,1)+halfWindowMSS))
                        partClassMSS(check(1)-1,2) = partClassMSS(check(1),2);
                        partClassMSS1D(check(1)-1,2) = partClassMSS1D(check(1),2);
                        partClassMSS(check(1),:) =[];
                        partClassMSS1D(check(1),:) =[];
                        mssSlope(check(1))=[];
                    else
                        partClassMSS(check(1)+1,1) = partClassMSS(check(1),1);
                        partClassMSS1D(check(1)+1,1) = partClassMSS1D(check(1),1);
                        partClassMSS(check(1),:) =[];
                        partClassMSS1D(check(1),:) =[];
                        mssSlope(check(1))=[];
                        
                    end
                end
                       
                unclassCheck  = isnan(partClassMSS(:,3)).*isnan(partClassMSS1D(:,3));
                check = find(unclassCheck);
            end
        end
end
function [pointClass, mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT,trackCenter,confRadius1D,prefDir] = asymmetricDiffusion(startPoint,endPoint,alphaAsym,probDim,tracks, iTrack)
                    trackCoordX = tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint)';
                    deltaCoordX = tracks(iTrack,8*(startPoint-1)+5:8:8*endPoint)';
                    trackCoordY = tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint)';
                    deltaCoordY = tracks(iTrack,8*(startPoint-1)+6:8:8*endPoint)';
                    trackCoordZ = tracks(iTrack,8*(startPoint-1)+3:8:8*endPoint)';
                    deltaCoordZ = tracks(iTrack,8*(startPoint-1)+7:8:8*endPoint)';
                    trackCoord = [trackCoordX trackCoordY trackCoordZ];
                    deltaCoord = [deltaCoordX deltaCoordY deltaCoordZ];
                    trackCoord = trackCoord(:,1:probDim);
                    deltaCoord = deltaCoord(:,1:probDim);

                    %project onto direction of motion
                    [posAlongDir,deltaPosAlongDir] = projectCoordOntoDir(trackCoord,...
                        deltaCoord,[],[]);

                    %construct matrix of linear tracks with projected positions
                    trackCoord2 = [posAlongDir zeros(length(posAlongDir),3) deltaPosAlongDir zeros(length(posAlongDir),3)]';
                    trackCoord2 = trackCoord2(:)';
            %         iAsym = iAsym + 1;
            %         tracksAsym(iAsym,:) = trackCoord2;
                    momentOrders = 0 : 6;
                    [pointClass,mssSlopeT,genDiffCoefT,scalingPowerT,normDiffCoefT] = ...
                    trackMSSAnalysis(trackCoord2,1,momentOrders,alphaAsym);
                %% Confinement of asym 
                    xyCoord = [trackCoordX trackCoordY];

                    %find the eignevalues of the variance-covariance matrix of this track's
                    %positions
                    [eigenVec,eigenVal] = eig(nancov(xyCoord(:,1:probDim)));
                    eigenVal = diag(eigenVal);

                    %calculate the confinement radius along the preferred direction of
                    %motion
                    confRadius1D(1,2) = sqrt( max(eigenVal) * (3) );

                    %calculate the confinement radius perpendicular to the preferred
                    %direction of motion
                    confRadius1D(1,1) = sqrt( mean(eigenVal(eigenVal~=max(eigenVal))) * (probDim + 1) );

                    %calculate the track's center
                    trackCenter = nanmean(xyCoord(:,1:probDim));

                    %store the preferred direction of motion
                    prefDir = eigenVec(:,eigenVal==max(eigenVal))';
end

function [confRadTmp,centerTmp] = estimConfRad(tracks,probDim,confRadMin)

%get subpart's coordinates
xCoord = tracks(1:8:end)';
yCoord = tracks(2:8:end)';
zCoord = tracks(3:8:end)';
xyzCoord = [xCoord yCoord zCoord];

%find the eigenvalues and eigenvectors of the variance-covariance
%matrix of this track's positions
eigenVal = eig(nancov(xyzCoord(:,1:probDim)));

%calculate the track's confinement radius
if confRadMin
    confRadTmp = sqrt( min(eigenVal) * (probDim + 2) );
else
    confRadTmp = sqrt( mean(eigenVal) * (probDim + 2) );
end

%calculate the track's center
centerTmp = nanmean(xyzCoord(:,1:probDim));
end
function [locPreTmp] = estimLocPre(tracks,probDim)

%get subpart's coordinates
xCoord = tracks(1:8:end)';
yCoord = tracks(2:8:end)';
zCoord = tracks(3:8:end)';
xyzCoord = [xCoord yCoord zCoord];

%find the eigenvalues and eigenvectors of the variance-covariance
%matrix of this track's positions
eigenVal = eig(nancov(xyzCoord(:,1:probDim)));

%calculate the track's confinement radius
    locPreTmp = sqrt( mean(eigenVal) );

end

function [segPoints,peakThresh] = findCloseSegments(trackFull,level,halfWindowMSS,windowMSSMin,windowMSSMin2)
        [peaks,locs] = findpeaks(trackFull);
      
        %Get segments and verify that each is a minimum length, if not then
        %join with adjacent segment with smallest peak. May need to repeat
        %after checking that segments are long enough. Make sure entire
        %track is retained
        offLim = [16:19];
        peakThresh = peaks(peaks >= level);
        segPoints = [1;locs(peaks >= level);length(trackFull)];
        segPoints(2:end)=segPoints(2:end)+1+halfWindowMSS; %half of windowMin, 1 for deriv offset
        segPoints(end)=segPoints(end)+1+halfWindowMSS; %half of windowMin, 1 for edn correction
        n = 1:length(segPoints)-1;
        difference = segPoints(n+1)-segPoints(n);
        checkT = find(difference < windowMSSMin2);
        checkM = find(difference >= 16);
        indMin = ismember(checkT,checkM);
        check = checkT(indMin);
    if ~isempty(check) 
        for m = 1:length(check) 
        %If short segment is first, connect to succeeding
            if check(m) ==1 && difference(check(m)) >= 16
                comp1 = difference(check(m));
                if difference(2) -(windowMSSMin2-comp1) >= windowMSSMin && ~ismember(difference(2) -(windowMSSMin2-comp1),offLim)
                    difference(1) = comp1 + (windowMSSMin2-comp1);
                    difference(2) = difference(2) -(windowMSSMin2-comp1);
% %                     checkT = find(difference < windowMSSMin2);
% %                     checkM = find(difference >= 16);
% %                     indMin = ismember(checkT,checkM);
% %                     check = checkT(indMin);
                    % Keep peak but move segPoint
                    segPoints(2,:) = segPoints(2,:) + (windowMSSMin2-comp1);
                    continue
% %                 else
% %                     segPoints(2,:) =[];
% %                     peakThresh(1) =[];
% % 
% %                     n = 1:length(segPoints)-1;
% %                     difference = segPoints(n+1)-segPoints(n);
% %                     checkT = find(difference < windowMSSMin);
% %                     ind = difference(difference < windowMSSMin);
% %                     indMin = difference(difference >= windowMSSMin);
% %                     check = checkT(ind  == min(ind));
% %                     continue
                end
                
            
            
             %If segment is last, connect to preceding
            elseif check(m) == length(difference) && difference(check(m)) >= 16
                comp1 = difference(check(m));
                if difference(check(m)-1) -(windowMSSMin2-comp1) >= windowMSSMin && ~ismember(difference(check(m)-1) -(windowMSSMin2-comp1),offLim)
                    difference(check(m)) = comp1 + (windowMSSMin2-comp1);
                    difference(check(m)-1) = difference(check(m)-1) -(windowMSSMin2-comp1);
% %                     checkT = find(difference < windowMSSMin2);
% %                     checkM = find(difference >= 16);
% %                     indMin = ismember(checkT,checkM);
% %                     check = checkT(indMin);
                    % Keep peak but move segPoint
                    segPoints(end-1,:) = segPoints(end-1,:) - (windowMSSMin2-comp1);
                    continue
% %                 else
% %                     segPoints(end-1,:) = [];
% %                     peakThresh(end) = [];
% %                     n = 1:length(segPoints)-1;
% %                     difference = segPoints(n+1)-segPoints(n);
% %                     checkT = find(difference < windowMSSMin);
% %                     ind = difference(difference < windowMSSMin);
% %                     indMin = difference(difference >= windowMSSMin);
% %                     check = checkT(ind  == min(ind));
% %                     continue
                end
                
            

            %If segment is somewhere in middle, connect to segment with lowest
            %peak
            elseif ~isempty(check) && difference(check(m)) >= 16
                if peakThresh(check(m)) > peakThresh(check(m)-1)
                    %if higher peak is in front, delete peak behind, same
                    %as if it were last segment
                    comp1 = difference(check(m));
                    if difference(check(m)-1) -(windowMSSMin2-comp1) >= windowMSSMin && ~ismember(difference(check(m)-1) -(windowMSSMin2-comp1),offLim)
                        difference(check(m)) = comp1 + (windowMSSMin2-comp1);
                        difference(check(m)-1) = difference(check(m)-1) -(windowMSSMin2-comp1);
% %                         checkT = find(difference < windowMSSMin2);
% %                         checkM = find(difference >= 16);
% %                         indMin = ismember(checkT,checkM);
                        
                        % Keep peak but move segPoint
                        segPoints(check(m),:) = segPoints(check(m),:) - (windowMSSMin2-comp1);
% %                         check = checkT(indMin);
                        continue
% %                     else
% %                         segPoints(end-1,:) = [];
% %                         peakThresh(end) = [];
% %                         n = 1:length(segPoints)-1;
% %                         difference = segPoints(n+1)-segPoints(n);
% %                         checkT = find(difference < windowMSSMin);
% %                         ind = difference(difference < windowMSSMin);
% %                         indMin = difference(difference >= windowMSSMin);
% %                         check = checkT(ind  == min(ind));
% %                         continue
                    end
                else
                    %if higher peak is behind, delete peak in front, same
                    %as if it were first segment
                        comp1 = difference(check(m));
                        if difference(check(m)+1) -(windowMSSMin2-comp1) >= windowMSSMin && ~ismember(difference(check(m)+1) -(windowMSSMin2-comp1),offLim)
                            difference(check(m)) = comp1 + (windowMSSMin2-comp1);
                            difference(check(m)+1) = difference(check(m)+1) -(windowMSSMin2-comp1);
% %                             checkT = find(difference < windowMSSMin2);
% %                             checkM = find(difference >= 16);
% %                             indMin = ismember(checkT,checkM);
                            % Keep peak but move segPoint
                            segPoints(check(m)+1,:) = segPoints(check(m)+1,:) + (windowMSSMin2-comp1);
% %                             check = checkT(indMin);
                            continue
% %                         else
% %                             segPoints(2,:) =[];
% %                             peakThresh(1) =[];
% % 
% %                             n = 1:length(segPoints)-1;
% %                             difference = segPoints(n+1)-segPoints(n);
% %                             checkT = find(difference < windowMSSMin);
% %                             ind = difference(difference < windowMSSMin);
% %                             indMin = difference(difference >= windowMSSMin);
% %                             check = checkT(ind  == min(ind));
% %                             continue
                        end
                end
                       
            end
        end
    end
end

function [segPoints,peakThresh] = findDiffSegments(segPoints,windowMSSMin, peakThresh) %trackFull,level,halfWindowMSS
%         [peaks,locs] = findpeaks(trackFull);%If this slow, try find(trackFull(n) > I(n+1) & trackFull(n) > I(n-1));
% %         idxSwitch = find(peaks >= level); %change to level prctile(trackFull,90)/ level
%         
%         %Get segments and verify that each is a minimum length, if not then
%         %join with adjacent segment with smallest peak. May need to repeat
%         %after checking that segments are long enough. Make sure entire
%         %track is retained 
%         peakThresh = peaks(peaks >= level);
%         segPoints = [1;locs(peaks >= level);length(trackFull)];
%         segPoints(2:end)=segPoints(2:end)+1+halfWindowMSS; %half of windowMin, 1 for deriv offset
%         segPoints(end)=segPoints(end)+1+halfWindowMSS; %half of windowMin, 1 for edn correction
        n = 1:length(segPoints)-1;
        difference = segPoints(n+1)-segPoints(n);
        checkT = find(difference < windowMSSMin);
        indMin = difference(difference >= windowMSSMin); %Is at least one segment analyzable?
        ind = difference(difference < windowMSSMin); %Temporarily
%         allowing short segments, uncomment to revert
        check = checkT(ind  == min(ind));
%         [~,check] =sort(difference(difference < windowMSSMin)); 
        while ~isempty(check) && isempty(indMin)
        %If short segment is first, connect to succeeding
            if check(1) == 1 && difference(check(1)) <16 %And less then a minimum length
                segPoints(2,:) =[];
                peakThresh(1) =[];
                
                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
% %                 difference(1) = difference(1) +1; %offset correction, verify
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                indMin = difference(difference >= windowMSSMin);
                check = checkT(ind  == min(ind));
%                 [~,check] =sort(difference(difference < windowMSSMin));
                continue
            elseif check(1) ==1 && difference(check(1)) >= 16
                comp1 = difference(check(1));
                if difference(2) -(windowMSSMin-comp1) >= windowMSSMin
                    difference(1) = comp1 + (windowMSSMin-comp1);
                    difference(2) = difference(2) -(windowMSSMin-comp1);
                    checkT = find(difference < windowMSSMin);
                    ind = difference(difference < windowMSSMin);
                    indMin = difference(difference >= windowMSSMin);
                    check = checkT(ind  == min(ind));
                    % Keep peak but move segPoint
                    segPoints(2,:) = segPoints(2,:) + (windowMSSMin-comp1);
                    continue
                else
                    segPoints(2,:) =[];
                    peakThresh(1) =[];

                    n = 1:length(segPoints)-1;
                    difference = segPoints(n+1)-segPoints(n);
            % %                 difference(1) = difference(1) +1; %offset correction, verify
                    checkT = find(difference < windowMSSMin);
                    ind = difference(difference < windowMSSMin);
                    indMin = difference(difference >= windowMSSMin);
                    check = checkT(ind  == min(ind));
                    continue
                end
                
            end
             %If segment is last, connect to preceding
            if check(1) == length(difference)  && difference(check(1)) <16
                segPoints(end-1,:) = [];
                peakThresh(end) = [];
                                                
                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
% %                 difference(1) = difference(1) +1; %offset correction, verify
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                indMin = difference(difference >= windowMSSMin);
                check = checkT(ind  == min(ind));
%                 [~,check] =sort(difference(difference < windowMSSMin));
                continue
            elseif check(1) == length(difference) && difference(check(1)) >= 16
                comp1 = difference(check(1));
                if difference(check(1)-1) -(windowMSSMin-comp1) >= windowMSSMin
                    difference(check(1)) = comp1 + (windowMSSMin-comp1);
                    difference(check(1)-1) = difference(check(1)-1) -(windowMSSMin-comp1);
                    checkT = find(difference < windowMSSMin);
                    ind = difference(difference < windowMSSMin);
                    indMin = difference(difference >= windowMSSMin);
                    check = checkT(ind  == min(ind));
                    % Keep peak but move segPoint
                    segPoints(end-1,:) = segPoints(end-1,:) - (windowMSSMin-comp1);
                    continue
                else
                    segPoints(end-1,:) = [];
                    peakThresh(end) = [];
                    n = 1:length(segPoints)-1;
                    difference = segPoints(n+1)-segPoints(n);
                    checkT = find(difference < windowMSSMin);
                    ind = difference(difference < windowMSSMin);
                    indMin = difference(difference >= windowMSSMin);
                    check = checkT(ind  == min(ind));
                    continue
                end
                
            end

            %If segment is somewhere in middle, connect to segment with lowest
            %peak
            if ~isempty(check)
                if peakThresh(check(1)) > peakThresh(check(1)-1)
                    segPoints(check(1),:) =[]; %if higher peak is in front, delete peak behind
                    peakThresh(check(1)-1)=[];
                else
                    segPoints(check(1)+1,:) =[];
                    peakThresh(check(1)) =[];
                end
                       
                n = 1:length(segPoints)-1;
                difference = segPoints(n+1)-segPoints(n);
                checkT = find(difference < windowMSSMin);
                ind = difference(difference < windowMSSMin);
                indMin = difference(difference >= windowMSSMin);
                check = checkT(ind  == min(ind));
%                 [~,check] =sort(difference(difference < windowMSSMin));
            end
        end
end
                    %{
                    %% Lower threshold Diffusion Analysis
                    %If any track has been classified as confined or
                    %immobile, check again with lower threshold
                    recheck = ismember(partClassMSS(:,3),[1,0]);
                    if sum(recheck) == 0
%                         continue
                    else
                        
                        tmpPartClassMSS = partClassMSS(recheck,:);
                        for m = 1:size(tmpPartClassMSS,1)
                          testPoints = segPoints2(segPoints2 >= tmpPartClassMSS(m,1) & segPoints2 <= tmpPartClassMSS(m,2)+1);
                          n = 1:length(testPoints)-1;
                          difference = testPoints(n+1)-testPoints(n);
                          pointClassMSS = NaN(length(n),1);
                          mssSlopeTmp = pointClassMSS;
                          normDiffCoefTmp = pointClassMSS;
                          confRadLow = pointClassMSS;
                          centerLow = NaN(length(n),probDim);
                          genDiffCoefTmp = NaN(length(n),length(momentOrders));
                          scalingPowerTmp = NaN(length(n),length(momentOrders));
                          partClassMSSTmp = NaN(length(n),3);
                          trackSELCurrent = trackSEL(iTrack,:);
                                  
                            for k = 1:length(testPoints)-1 
                            startPoint = trackSELCurrent(1) +testPoints(k)-1;
                            endPoint  = startPoint+ difference(k)-1;

                                        [pointClassMSS(k),mssSlopeTmp(k),genDiffCoefTmp(k,:),scalingPowerTmp(k,:),normDiffCoefTmp(k)] = trackMSSAnalysis(...
                                            tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                                            probDim,momentOrders,alphaValues);

                                        if pointClassMSS(k) <= 1
                                            [confRadLow(k),centerLow(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                                            xCoord = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                            yCoord = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                                            center  = centerLow(k,:);
% %                                             [pointClassMSS(k)]  = immobileDetection(xCoord,yCoord,center);
                                        else
                                            confRadLow(k) = NaN;
                                            centerLow(k,:) = NaN(1,probDim);
                                        end 

                                        partClassMSSTmp(k,1)= startPoint;
                                        partClassMSSTmp(k,2)= endPoint;
                                        partClassMSSTmp(k,3)= pointClassMSS(k);
                                        %Rest of MSS information or not... jsut fill in track

                            end
                            % Connect like segments again
                            %partClassTmp = partClassMSSTmp;
                                for iSubpart = 1 : length(pointClassMSS)-1
                                    iSubpartPlus1 = iSubpart + 1;
                                    
                                    while ( (iSubpartPlus1 <= length(pointClassMSS)) && ...
                                    ( (pointClassMSS(iSubpart) == pointClassMSS(iSubpartPlus1)) || ...
                                    (isnan(pointClassMSS(iSubpart)) && isnan(pointClassMSS(iSubpartPlus1))) ) )
                                    partClassMSSTmp(iSubpart,2) = partClassMSSTmp(iSubpartPlus1,2);
                                    partClassMSSTmp(iSubpartPlus1,1) = partClassMSSTmp(iSubpart,1);
                                    iSubpartPlus1 = iSubpartPlus1 + 1;
                                    end
                                end
                                
                            [~,uniqueParts] = unique(partClassMSSTmp(:,1));
                            partClassMSSTmp = partClassMSSTmp(uniqueParts,:);
                            
                            %create no switch option, nothing changed
                            
                            if sum(ismember(partClassMSSTmp(:,3),2))>0
                                partClassMSSTest =partClassMSS;
%                                 continue
                            elseif ismember(partClassMSSTmp,partClassMSS,'rows')
                                partClassMSSTest =partClassMSS;
%                                 continue
                            else
                                if m ==1
                                    oldClass =partClassMSS(~recheck,:);
%                                     mssSlope0=mssSlope(~recheck,:);
%                                     genDiffCoef0=genDiffCoef(~recheck,:);
%                                     scalingPower0 =scalingPower(~recheck,:);
%                                     normDiffCoef0=normDiffCoef(~recheck,:);
%                                     confRadTmp0=confRadTmp(~recheck,:);
%                                     centerTmp0=centerTmp(~recheck,:);
                                else
                                        oldClass = partClassMSS;

%                                     mssSlope0=mssSlope;
%                                     genDiffCoef0=genDiffCoef;
%                                     scalingPower0 =scalingPower;
%                                     normDiffCoef0=normDiffCoef;
%                                     confRadTmp0=confRadTmp;
%                                     centerTmp0=centerTmp;
                                end
%                                 clear partClassMSS
                                partClassMSSTest = [oldClass;partClassMSSTmp];
%                                 mssSlope=[mssSlope0;mssSlopeTmp];
%                                 genDiffCoef=[genDiffCoef0;genDiffCoefTmp];
%                                 scalingPower =[scalingPower0;scalingPowerTmp];
%                                 normDiffCoef=[normDiffCoef0;normDiffCoefTmp];
%                                 confRadTmp=[confRadTmp0;confRadLow];
%                                 centerTmp=[centerTmp0;centerLow];
                               %Sort out the rest of the variables being
                               %saved
                               %Eventually accomodate for many segments to
                               %be re-analyzed
                            end
                            
                        end
                       partClassMSS = partClassMSSTest; 
                    end
                    
                     %% II Reclassify again if segments next to each other are the same
                     %merge subparts that now have the same classification
                            partClassTmp = partClassMSS;
                            pointClassMSS = partClassMSS(:,3); 
                    for iSubpart = 1 : length(pointClassMSS)-1
                        iSubpartPlus1 = iSubpart + 1;
                        while ( (iSubpartPlus1 <= length(pointClassMSS)) && ...
                                ( (pointClassMSS(iSubpart) == pointClassMSS(iSubpartPlus1)) || ...
                                (isnan(pointClassMSS(iSubpart)) && isnan(pointClassMSS(iSubpartPlus1))) ) )
                            partClassMSS(iSubpart,2) = partClassMSS(iSubpartPlus1,2);
                            partClassMSS(iSubpartPlus1,1) = partClassMSS(iSubpart,1);
                            iSubpartPlus1 = iSubpartPlus1 + 1;
                        end
                    end
                    [~,uniqueParts] = unique(partClassMSS(:,1));
                    partClassMSS = partClassMSS(uniqueParts,:);
                    %Final Diffusion characteristics

                        numSeg = size(partClassMSS,1);
                        pointClassMSS = NaN(numSeg,1);
                        mssSlope = pointClassMSS;
                        normDiffCoef = pointClassMSS;
                        confRadTmp = pointClassMSS;
                        centerTmp = NaN(numSeg,probDim);
                        genDiffCoef = NaN(numSeg,length(momentOrders));
                        scalingPower = NaN(numSeg,length(momentOrders));
                        
                        
                        for k = 1:numSeg 
                        startPoint = partClassMSS(k,1);
                        endPoint  = partClassMSS(k,2);

                                    [pointClassMSS(k),mssSlope(k),genDiffCoef(k,:),scalingPower(k,:),normDiffCoef(k)] = trackMSSAnalysis(...
                                        tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                                        probDim,momentOrders,alphaValues);

                                    if pointClassMSS(k) <= 1
                                        [confRadTmp(k),centerTmp(k,:)] = estimConfRad(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,confRadMin);
                                        xCoord = (tracks(iTrack,8*(startPoint-1)+1:8:8*endPoint))';
                                        yCoord = (tracks(iTrack,8*(startPoint-1)+2:8:8*endPoint))';
                                        center  = centerTmp(k,:);
% %                                         [pointClassMSS(k)]  = immobileDetection(xCoord,yCoord,center);
                                    else
                                        confRadTmp(k) = NaN;
                                        centerTmp(k,:) = NaN(1,probDim);
                                    end 
                            partClassMSS(k,3) = pointClassMSS(k);
                        end
                    
%}