function plotTracksTransDiffAnalysis2D(trackedFeatureInfo,diffAnalysisRes,timeRange,...
    newFigure,image,showConf,asymCheck,bayes)
%PLOTTRACKSTRANSDIFFANALYSIS plots tracks in 2D highlighting the different diffusion segments within each track
%
%SYNOPSIS plotTracksTransDiffAnalysis2D(trackedFeatureInfo,diffAnalysisRes,timeRange,...
%    newFigure,image,showConf)
%
%INPUT  trackedFeatureInfo: -- EITHER -- 
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%                           -- OR -- 
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       diffAnalysisRes   : Structure of diffusion analysis results as
%                           output by the code "trackTransDiffusionAnalysis1".
%       timeRange         : 2-element row vector indicating time range to plot. 
%                           Optional. Default: whole movie.
%       newFigure         : 1 if plot should be made in a new figure
%                           window, 0 otherwise (in which case it will be
%                           plotted in an existing figure window).
%                           Optional. Default: 1.
%       image             : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image.
%       showConf          : 1 to show confinement radii, 0 otherwise.
%                           Optional. Default: 0.
%
%OUTPUT The plot.
%       Color coding:
%       linear according to asymmetry -> red
%       random/unclassified & 2D confined diffusion -> blue
%       random/unclassified & 2D normal diffusion -> cyan
%       random/unclassified & 2D super diffusion -> magenta
%       random & too short to analyze 2D diffusion -> purple
%       too short for any analysis -> black
%
%Khuloud Jaqaman, April 2009
%Edited: Tony Vega 2016

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--plotTracksTransDiffAnalysis2D: Incorrect number of input arguments!');
    return
end

%get number of tracks and number of time points
if isstruct(trackedFeatureInfo) %if tracks are in structure format
    numTracks = length(trackedFeatureInfo);
    tmp = vertcat(trackedFeatureInfo.seqOfEvents);
    numTimePoints = max(tmp(:,1));
    clear tmp
else %if tracks are in matrix format
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;
end

errFlag = 0;

%check whether a time range for plotting was input
if nargin < 3 || isempty(timeRange)
    timeRange = [1 numTimePoints];
else
    if timeRange(1) < 1 || timeRange(2) > numTimePoints
        disp('--plotTracksTransDiffAnalysis2D: Wrong time range for plotting!');
        errFlag = 1;
    end
end

%check whether newFigure was input
if nargin < 4 || isempty(newFigure)
    newFigure = 1;
else
    if newFigure ~= 0 && newFigure ~= 1
        disp('--plotTracksTransDiffAnalysis2D: newFigure should be 0 or 1!');
        errFlag = 1;
    end
end

%check whether user supplied an image
if nargin < 5 || isempty(image)
    image = [];
end

%check whether to plot confinement radii
if nargin < 6 || isempty(showConf)
    showConf = 0;
end

%exit if there are problems in input variables
if errFlag
    disp('--plotTracksTransDiffAnalysis2D: Please fix input data!');
    return
end

%% Pre-processing

if isstruct(trackedFeatureInfo) %if tracks are input in structure format

    %store the input structure as a variable with a different name
    inputStructure = trackedFeatureInfo;
    clear trackedFeatureInfo;
    
    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end

    %if all tracks have only one segment ...
    if max(numSegments) == 1

        %indicate that there are no compound tracks with merging and splitting branches
        mergeSplit = 0;

        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step) 
        %in this case of course every compound track is simply one track
        %without branches
        trackStartRow = (1:numTracks)';

        %store tracks in a matrix
        trackedFeatureInfo = NaN*ones(numTracks,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(i,8*(startTime-1)+1:8*endTime) = inputStructure(i).tracksCoordAmpCG;
        end
        
    else %if some tracks have merging/splitting branches
        
        %indicate that in the variable mergeSplit
        mergeSplit = 1;
        
        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        trackStartRow = ones(numTracks,1);
        for iTrack = 2 : numTracks
            trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);            
        end
        
        %put all tracks together in a matrix
        trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            trackedFeatureInfo(trackStartRow(i):trackStartRow(i)+...
                numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
                inputStructure(i).tracksCoordAmpCG;
        end
        
    end    
    
else %if tracks are not input in structure format

    %indicate that there are no compound tracks with merging and splitting branches
    mergeSplit = 0;
    
    %locate the row of the first track of each compound track in the
    %big matrix of all tracks
    %in this case of course every compound track is simply one track
    %without branches
    trackStartRow = (1:numTracks)';

end

%get the x,y-coordinates of features in all tracks
tracksX = trackedFeatureInfo(:,1:8:end)';
tracksY = trackedFeatureInfo(:,2:8:end)';

%find x-coordinate limits
minXCoord = min(floor(min(tracksX(:))),0);
maxXCoord =  ceil(max(tracksX(:)));

%find y-coordinate limits
minYCoord = min(floor(min(tracksY(:))),0);
maxYCoord =  ceil(max(tracksY(:)));

%get number of track segments to be plotted
numTrackSegments = size(tracksX,2);

%% confinement radius information

% %get track segment center, confinement radii and preferred direction of
% %motion
% trackSegmentCenter = catStruct(1,'diffAnalysisRes.confRadInfo.trackCenter');
% trackSegmentConfRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius');
% trackSegmentPrefDir = catStruct(1,'diffAnalysisRes.confRadInfo.prefDir');
% 
% %determine indices of tracks with one confinement radius
% indxCircle = find( ~isnan(trackSegmentConfRad(:,1)) & isnan(trackSegmentConfRad(:,2)) );
% 
% %determine indices of tracks with 2 confinement radii
% indxRectangle = find( ~isnan(trackSegmentConfRad(:,2)) );

%% Plotting

%if the user wants to plot in a new figure window
if newFigure

    %open new figure window
    figure

    if ~isempty(image) %if user supplied an image
        imshow(image,[]); %plot the image
    else %if user did not supply an image
        imshow(ones(maxYCoord,maxXCoord),[]); %plot an empty image
    end

    %set figure axes limits
    axis([minXCoord maxXCoord minYCoord maxYCoord]);

    %show coordinates on axes
    ah = gca;
    set(ah,'visible','on');

    %label axes
    xlabel('x-coordinate (pixels)');
    ylabel('y-coordinate (pixels)');

end

%hold on figure
hold on

%extract the portion of tracksX and tracksY that is of interest
tracksXP = tracksX(timeRange(1):timeRange(2),:);
tracksYP = tracksY(timeRange(1):timeRange(2),:);

%% New plotting strategy
%Build new matrix which indicates when tracks change diffusion


base = NaN([size(tracksX,1) size(tracksX,2) 4]); %Create 3D matrix to store indices for all diffusion types
asymBase = NaN([size(tracksX,1) size(tracksX,2) 4]);
nanBase = NaN(size(tracksX));
diffTypesF = [];
 if isempty(bayes)
        %get track segment types from diffusion analysis
        trackSegmentType = vertcat(diffAnalysisRes.segmentClass);

        if asymCheck == 1
            for k = 1:length(trackSegmentType)

                for j = 1:size(trackSegmentType(k).momentScalingSpectrum,1)
                base(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k) = ...
                trackSegmentType(k).momentScalingSpectrum(j,3);

                asymBase(trackSegmentType(k).momentScalingSpectrum1D(j,1):trackSegmentType(k).momentScalingSpectrum1D(j,2),k) = ...
                trackSegmentType(k).momentScalingSpectrum1D(j,3);

                end

            end
        else
            for k = 1:length(trackSegmentType)

                for j = 1:size(trackSegmentType(k).momentScalingSpectrum,1)
                    if ~isnan(trackSegmentType(k).momentScalingSpectrum(j,3))
                base(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k,trackSegmentType(k).momentScalingSpectrum(j,3)+1) = ...
                trackSegmentType(k).momentScalingSpectrum(j,3);
            diffTypesF = [diffTypesF ; trackSegmentType(k).momentScalingSpectrum(j,3)];
                    else
                nanBase(trackSegmentType(k).momentScalingSpectrum(j,1):trackSegmentType(k).momentScalingSpectrum(j,2),k) = -1;
                    end
                end

            end
        end
 else
   for m = 1:size(tracksXP,2)
        results = diffAnalysisRes;
        states = results.ML_states{1,m};
        track = results.track{1,m};
        steps = NaN(2,size(track,2));
        steps(:,2:end) = results.steps{1,m};

        track(1,~isnan(steps(1,:))) = states;
        track(2,~isnan(steps(1,:))) = states;
        index = find(isnan(steps(2,:)));
        track(1,index)=track(2,index+1);
        track(2,index)=track(2,index+1);
        base(~isnan(tracksXP(:,m)),m) = track(1,:);
    end  
     
 end

numTimePlot = timeRange(2) - timeRange(1) + 1;


% asymBase(isnan(asymBase))=-1;
% base(isnan(base))=-1;
% base(base==-2)=NaN;
% asymBase(asymBase==-2)=NaN;
diffTypes = unique(diffTypesF);
% diffTypes = unique(base);
copyR = size(tracksXP,1);
copyC = size(tracksXP,2);
removeGapParams.Delimeter = Inf;
removeGapParams.RemoveOtherGaps = true;
removeGapParams.Color = 'k';
removeGapParams.LineStyle = ':';
lineWithGaps(tracksXP,tracksYP,removeGapParams);
for k = 1:length(diffTypes)
 	
    switch diffTypes(k)
        case -1
            ind = find(nanBase(:) ==-1);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
        
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 0 0]);
            end %Allow to change to white if image involved
           %asym linear and unclassified
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
                if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                    for i=1:numTimePlot-1
                        validData=~all(isnan(copyX(i:i+1,:)),1);
                        lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 0 1]);%[1 1 0][0 102/255 51/255]
                    end 
                end
           %symmetric and unclassified
            indA = find(asymBase ==0);
            indAF  = intersect(indA,ind);
                if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                    for i=1:numTimePlot-1
                        validData=~all(isnan(copyX(i:i+1,:)),1);
                        lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0.6 0 1]);
                    end 
                end
        case 0
            ind = find(base(:,:,1) ==0);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
 	
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0.5 0.3 0]);
            end
            %Currently not labeling immobile and linear, doesn't make sense
            %to me
        case 1
 	
            ind = find(base(:,:,2) ==1);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);

            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 0 1]);%Real
% %                 lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color','c');
            end
            %Asym linear and confined
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
            if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                for i=1:numTimePlot-1
                    validData=~all(isnan(copyX(i:i+1,:)),1);
                    lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[1 0.7 0]);
                end 
            end
        case 2
 	
            ind = find(base(:,:,3) ==2);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
 	
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color','c');%Real
% %                 lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 0 1]);
            end
           %asym linear and free
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
            if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                for i=1:numTimePlot-1
                    validData=~all(isnan(copyX(i:i+1,:)),1);
                    lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[1 0 0]);
                end 
            end
        case 3
            ind = find(base(:,:,4) ==3);
            copyX = NaN(copyR,copyC);
            copyY = NaN(copyR,copyC);
            copyX(ind) = tracksXP(ind);
            copyY(ind) = tracksYP(ind);
 	
            for i=1:numTimePlot-1
                validData=~all(isnan(copyX(i:i+1,:)),1);
                lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color','m');
            end
 	           %asym linear and super
            indA = find(asymBase ==1);
            indAF  = intersect(indA,ind);
            if ~isempty(indAF)
                    copyX = NaN(copyR,copyC);
                    copyY = NaN(copyR,copyC);
                    copyX(indAF) = tracksXP(indAF);
                    copyY(indAF) = tracksYP(indAF);

                for i=1:numTimePlot-1
                    validData=~all(isnan(copyX(i:i+1,:)),1);
                    lineWithGaps(copyX(i:i+1,validData),copyY(i:i+1,validData), 'Color',[0 1 0]);
                end 
            end
    end
 	
end
%{
%plot track segments with their appropriate diffusion type colors
for i = 1 : numTrackSegments

    %missing intervals are indicated by a dotted black line
    obsAvail = find(~isnan(tracksXP(:,i)));
    plot(tracksXP(obsAvail,i),tracksYP(obsAvail,i),'k:');
    
    %get the classification of this track segment
    segmentClassMSS  = trackSegmentType(i).momentScalingSpectrum(:,1:3);
    segmentClassAsym = trackSegmentType(i).asymmetryAfterMSS;
    segmentClass = [segmentClassMSS segmentClassAsym(:,3)];
    segmentConfInfo = trackSegmentType(i).momentScalingSpectrum(:,end-2:end);
    
    %remove any parts of the track segment classification before the
    %plotting time interval
    while segmentClass(1,2) < timeRange(1)
        segmentClass = segmentClass(2:end,:);
    end
    segmentClass(1,1) = max(segmentClass(1,1),timeRange(1));
    
    %remove any parts of the track segment classification after the
    %plotting time interval
    while segmentClass(end,1) > timeRange(2)
        segmentClass = segmentClass(1:end-1,:);
    end
    segmentClass(end,2) = min(segmentClass(end,2),timeRange(2));
    
    %shift the track segment classification to make it start at 1
    segmentClass(:,1:2) = segmentClass(:,1:2) - timeRange(1) + 1;
    
    %substract 1 from the starts of track segments after the first to
    %prevent discontinuities in the plotting
    segmentClass(2:end,1) = segmentClass(2:end,1) - 1;

    %plot the different parts with their appropriate colors
    for iPart = 1 : size(segmentClass,1)

        switch segmentClass(iPart,3)

            case -1 %linear based on asymmetry
                plot(tracksXP(segmentClass(iPart,1):segmentClass(iPart,2),i),...
                    tracksYP(segmentClass(iPart,1):segmentClass(iPart,2),i),'r');

            case 1 %confined based on MSS
                plot(tracksXP(segmentClass(iPart,1):segmentClass(iPart,2),i),...
                    tracksYP(segmentClass(iPart,1):segmentClass(iPart,2),i),'b');
                
                %plot confinement radius if requested
                if showConf
                    theta = (0:pi/10:2*pi); %angle
                    xy = [cos(theta') sin(theta')]; %x and y-coordinates
                    circleVal = xy .* segmentConfInfo(iPart,1);
                    plot(segmentConfInfo(iPart,2)+circleVal(:,1),...
                        segmentConfInfo(iPart,3)+circleVal(:,2),'k');
                end

            case 2 %Brownian based on MSS
                plot(tracksXP(segmentClass(iPart,1):segmentClass(iPart,2),i),...
                    tracksYP(segmentClass(iPart,1):segmentClass(iPart,2),i),'c');

            case 3 %directed based on MSS
                plot(tracksXP(segmentClass(iPart,1):segmentClass(iPart,2),i),...
                    tracksYP(segmentClass(iPart,1):segmentClass(iPart,2),i),'m');

            otherwise
                
                switch segmentClass(iPart,4)
                    
                    case 0 %random based on asymmetry but too short for MSS analysis
                        plot(tracksXP(segmentClass(iPart,1):segmentClass(iPart,2),i),...
                            tracksYP(segmentClass(iPart,1):segmentClass(iPart,2),i),'Color',[0.6 0 1]);

                    otherwise %too short fo any analysis
                        plot(tracksXP(segmentClass(iPart,1):segmentClass(iPart,2),i),...
                            tracksYP(segmentClass(iPart,1):segmentClass(iPart,2),i),'k');

                end %(switch segmentClass(iPart,4))

        end %(switch segmentClass(iPart,3))

    end %(for iPart = 1 : size(segmentClass,1))

end %(for i = 1 : numTrackSegments)
%}
%show merges and splits
if mergeSplit

    %go over all tracks
    for iTrack = 1 : numTracks

        %parse sequence of events of this compound track and find merges and
        %splits
        seqOfEvents = inputStructure(iTrack).seqOfEvents;
        indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';
        indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4)) ...
            & seqOfEvents(:,1) > timeRange(1) & seqOfEvents(:,1) <= timeRange(2)))';

        %go over all splits
        for iSplit = indxSplit

            %get time of splitting
            timeSplit = seqOfEvents(iSplit,1);

            %determine row where starting track is located
            rowS = trackStartRow(iTrack) + seqOfEvents(iSplit,3) - 1;

            %determine row where splitting track is located
            rowSp = trackStartRow(iTrack) + seqOfEvents(iSplit,4) - 1;

            %plot split as a dash-dotted line
            plot([tracksX(timeSplit,rowS) tracksX(timeSplit-1,rowSp)], ...
                [tracksY(timeSplit,rowS) tracksY(timeSplit-1,rowSp)],'k-.');

        end

        %go over all merges
        for iMerge = indxMerge

            %get time of merging
            timeMerge = seqOfEvents(iMerge,1);

            %determine row where ending track is located
            rowE = trackStartRow(iTrack) + seqOfEvents(iMerge,3) - 1;

            %determine row where merging track is located
            rowM = trackStartRow(iTrack) + seqOfEvents(iMerge,4) - 1;

            %plot merge as a dashed line
            plot([tracksX(timeMerge-1,rowE) tracksX(timeMerge,rowM)], ...
                [tracksY(timeMerge-1,rowE) tracksY(timeMerge,rowM)],'k--');

        end

    end %(for iTrack = 1 : numTracks)

end %(if mergeSplit)

% %show confinement areas if requested
% if showConf
% 
%     %generate circle to plot
%     theta = (0:pi/10:2*pi); %angle
%     xy = [cos(theta') sin(theta')]; %x and y-coordinates
%     
%     %go over symmetric tracks
%     for iTrack = indxCircle'
%         
%         %plot a circle of radius = confinement radius and centered at the
%         %center of this track
%         circleVal = xy .* trackSegmentConfRad(iTrack,1);
%         plot(trackSegmentCenter(iTrack,1)+circleVal(:,1),...
%             trackSegmentCenter(iTrack,2)+circleVal(:,2),'k');
%         
%     end
%     
%     %go over linear tracks
%     for iTrack = indxRectangle'
%     
%         %get the confinement axes
%         axisPara = trackSegmentPrefDir(iTrack,:);
%         axisPerp = [-axisPara(2) axisPara(1)] * trackSegmentConfRad(iTrack,1);
%         axisPara = axisPara * trackSegmentConfRad(iTrack,2);
%         
%         %find the 4 corners of the confinement rectangle
%         cornerCoord = [-axisPara - axisPerp; -axisPara + axisPerp; ...
%             axisPara + axisPerp; axisPara - axisPerp; -axisPara - axisPerp] ...
%             + repmat(trackSegmentCenter(iTrack,:),5,1);
%         
%         %plot the rectangle
%         plot(cornerCoord(:,1),cornerCoord(:,2),'k');
% 
%     end
%     
% end

%%%%% ~~ the end ~~ %%%%%

