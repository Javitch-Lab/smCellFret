
function donor = findDonorTracks(inFileDonA, inFileDonB)
%--------------------------------------------------------------------------
%
% findDonorTracks.m:
%   The function concatenates mapped donor tracks with real single molecule 
%   tracks in the donor channel.
% 
% Description:
%   To generate smFRET trajectories, we use a locally sensitive 
%   transformation function that maps the acceptor position to the donor 
%   channel at each time point to identify the attenuated donor signal 
%   associated with each acceptor particle during FRET. Therefore, the 
%   donor traces during FRET are mapped traces derived from those of the
%   acceptor. After the loss of the acceptor signal, the mapped donor track
%   is concatenated, if possible, with the nearest matching donor track 
%   available within a given search radius at the position of the acceptor
%   loss to obtain the continuation of the corresponding donor track. This 
%   is accomplished by the findDonorTracks.m function, which uses the 
%   sub-function E = distance(A,B) [1] to calculate the nearest neighbor 
%   distances between mapped and real tracks. 
% 
% Syntax:  
%   donor = findDonorTracks(inFileDonA, inFileDonB)
% 
% Inputs:
%   1. inFileDonA - Files of type 'singleTracesDonA_phc.traces' containing
%                   the fretTraces structure for mapped donor traces only. 
%   2. inFileDonB - Files of type 'singleTracesDonB_phc.traces' containing
%                   the fretTraces structure for all tracked donor traces.
% 
% Outputs:
%   1. donor - The variable is no longer used.
%   2. The function saves the two following files to the #Ch2 working 
%      directory: 
%      a)'singleTracesDon_phc.traces', a file containing the fretTraces 
%        structure of the donor data. In this dataset, mapped donor tracks 
%        are now concatenated after acceptor photobleaching with the 
%        nearest matching neighbor donor track, if possible. 
%      b)'singleTracesDonAOld_phc.traces', a backup file of the source 
%        file 'singleTracesDonA_phc.traces'. 
% 
% Other m-files required: 
%   Subfunctions: distance.m [1]
% 
% See also: 
%   scriptGetFretTraces.m
%
% References:
%   [1] The function E = distance(A,B) was written by  
%       Roland Bunschoten
%       University of Amsterdam
%       Intelligent Autonomous Systems (IAS) group
%       Kruislaan 403  1098 SJ Amsterdam
%       tel.(+31)20-5257524
%       bunschot@wins.uva.nl
%       Last Rev : Oct 29 16:35:48 MET DST 1999
%       Tested   : PC Matlab v5.2 and Solaris Matlab v5.3
%       Thanx    : Nikos Vlassis
%       Copyright notice: You are free to modify, extend and distribute 
%       this code granted that the author of the original code is 
%   	mentioned as the original author of the code.
%
% Authors: 
%   - P.G. Feb 2014
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%
%% Load files
%clear
%currentDir='E:\#BackupMicData\140905_pc5-mCMV-SNAPf-D2s-liveFRET_Clone32\967nMAlexa647_33nMDY549_PcaPcd_140ngml-Tet\test\';
%dataSetCh2='#01Ch2\';
%inFileDonA='singleTracesCh2_phc.traces';
%inFileDonB='singleTracesCh2b_phc.traces';
%tracesDonA=loadTracesCell([currentDir dataSetCh2 inFileDonA]);
%tracesDonB=loadTracesCell([currentDir dataSetCh2 inFileDonB]);

tracesDonA=loadTracesCell(inFileDonA);
tracesDonB=loadTracesCell(inFileDonB);

%% Initialize Parameters
minTraceLen           = 5; % minimum trace length of tracked donor that is 
                           % evaluated
pixCal                = 0.1576; % unit um/pixel 
nnLimit               = 250;     % unit nm 
dt                    = tracesDonA.time(2)- tracesDonA.time(1);
[nTracesD1,nFramesD1] = size(tracesDonA.x);
[nTracesD2,~]         = size(tracesDonB.x);

%% Find Endpoint of Donor Trace1 (mapped donor trace)

for i=1:nTracesD1
    endOfTrace   = tracesDonA.traceMetadata(i).endOfTrace;
    lenBaseline  = tracesDonA.traceMetadata(i).lenBaseline;
    Tbgd = tracesDonA.time(endOfTrace:endOfTrace+lenBaseline);
    if Tbgd(1)-dt > 0  &&  Tbgd(1)+dt <= nFramesD1*dt
       frames = ceil([Tbgd(1)-dt Tbgd(1) Tbgd(1)+dt]/dt);
    else 
       continue
    end
    % End of donorTrace1
    endOfD1 = frames(end);
    
    % Save the interval end of donorTrace1 in a vector
    xDonor1 = tracesDonA.x(i,frames);
    yDonor1 = tracesDonA.y(i,frames);
    if sum(xDonor1)==0 && sum(yDonor1)==0;continue;end
    D1  = [xDonor1; yDonor1];
    
    % Allocate memory for nnTraceIdx matrix
    
    nnTrace=0;
    nnTraceIdx = zeros(10,4);
    nnTraceIdx(:,:) = NaN;
    
    %Compare traceDonA (mapped donor trace) with the part of traceDonB 
    %(tracked donor trace) that contains the actual data (not the bgd data)
    for j=1:nTracesD2
        startOfTrace = tracesDonB.traceMetadata(j).startOfTrace;
        endOfTrace   = tracesDonB.traceMetadata(j).endOfTrace;
        lenBaseline  = tracesDonB.traceMetadata(j).lenBaseline;
        traceLen     = tracesDonB.traceMetadata(j).traceLen;
        
        % Only evaluate traces with a min trace length
        if traceLen < minTraceLen, continue; end
        staOfD2 = startOfTrace;
        endOfD2 = endOfTrace+lenBaseline;
        
        % Make sure that the actual data of tracesDonA and tracesDonB overlap 
        if staOfD2 < endOfD1 && endOfD2 > endOfD1
            xDonor2 = tracesDonB.x(j,frames);
            yDonor2 = tracesDonB.y(j,frames);
            D2  = [xDonor2; yDonor2];
        else
            continue
        end
        
        %Calculate Distance Matrix for the number of frames defined above 
        dist = distance(D1,D2);
        [minDist,Idx] = min(dist(:));
        [row,col] = ind2sub(size(dist),Idx);
        
        % Minimum Distance in nm
        minDist=minDist*pixCal*1000;
       
        % Determine the number NN-traces with a distance < nnLimit
        if minDist < nnLimit && nnTrace <= size(nnTraceIdx,1)
            nnTrace = nnTrace + 1;
            %disp(['Trace :' num2str(j) 'Dist :' num2str(minDist)]);
            nnTraceIdx(nnTrace,1) = j;
            nnTraceIdx(nnTrace,2) = minDist;
            nnTraceIdx(nnTrace,3) = row;
            nnTraceIdx(nnTrace,4) = col;
        end
    end

    % Find the nnTrace with the smallest distance in the nnTraceIdx
    if ~isnan(nnTraceIdx(1,1))
        [~,minInd] = min(nnTraceIdx(:,2));
        j=nnTraceIdx(minInd,1);
        k=nnTraceIdx(minInd,3);
        
        %i
        endOfTraceDonA=tracesDonA.traceMetadata(i).endOfTrace;
        %j
        endOfTraceDonB=tracesDonB.traceMetadata(j).endOfTrace;
        
        %if startOfTrace1 <= startOfTrace2 
            %lenBackset   = tracesDonA.traceMetadata(i).lenBackset;
        
        %if endOfTrace1 >= endOfTrace2
            %endOfTrace = endOfTrace1;
            %endOfTrace   = tracesDonA.traceMetadata(i).endOfTrace;
            %lenBaseline  = tracesDonA.traceMetadata(i).lenBaseline;
            
        if endOfTraceDonA <= endOfTraceDonB
            tracesDonA.x(i,frames(k):end) = tracesDonB.x(j,frames(k):end);
            tracesDonA.y(i,frames(k):end) = tracesDonB.y(j,frames(k):end);
            tracesDonA.int(i,frames(k):end) = tracesDonB.int(j,frames(k):end);
            
            startOfTrace = tracesDonA.traceMetadata(i).startOfTrace;
            endOfTrace   = tracesDonB.traceMetadata(j).endOfTrace;
            lenBaseline  = tracesDonB.traceMetadata(j).lenBaseline;
            
            tracesDonA.traceMetadata(i).ids          = tracesDonB.traceMetadata(j).ids; 
            tracesDonA.traceMetadata(i).endOfTrace   = endOfTrace; 
            tracesDonA.traceMetadata(i).lenBaseline  = lenBaseline;
            tracesDonA.traceMetadata(i).traceLen     = endOfTrace - ...
                                                     startOfTrace + 1;
        else
            tracesDonA.traceMetadata(i).ids          = 'notAssigned'; 
        end
        %elseif startOfTrace1 > startOfTrace2
        %    lenBackset   = tracesDonB.traceMetadata(j).lenBackset;
        %    startOfTrace = tracesDonB.traceMetadata(j).startOfTrace;
        %    endOfTrace   = tracesDonB.traceMetadata(j).endOfTrace;
        %    lenBaseline  = tracesDonB.traceMetadata(j).lenBaseline;
        %end
        %tracesDonA.traceMetadata(i).lenBackset   = lenBackset;
        %tracesDonA.traceMetadata(i).startOfTrace = startOfTrace;
        %tracesDonA.traceMetadata(i).endOfTrace   = endOfTrace; 
        %tracesDonA.traceMetadata(i).lenBaseline  = lenBaseline;
        %tracesDonA.traceMetadata(i).traceLen     = endOfTrace - ...
        %                                             startOfTrace + 1
    else
        tracesDonA.traceMetadata(i).ids          = 'notAssigned';
    end
end
donor = tracesDonA;

%% Save the modified *.phc traces file 

%Rename original Ch2_phc file before overwriting
currentDir=fileparts(inFileDonA);
newFileDonA = 'singleTracesDonAOld_phc.traces';
oldFile = inFileDonA;
newFile = [currentDir '\' newFileDonA];
movefile(oldFile, newFile);

% Save modified Ch2_phc file
currentDir=fileparts(inFileDonA);
outFileDonA = 'singleTracesDon_phc.traces';
saveTracesCell([currentDir '\' outFileDonA],'traces',tracesDonA);


