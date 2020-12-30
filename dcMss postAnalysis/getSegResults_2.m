function [SegResultsFinal,SegStats]=getSegResults_2(results)
%--------------------------------------------------------------------------
%
% getSegResults_2:
%   subfunction of scriptGetSegResults.m
% 
% Description:
%   This function combines all track segment data of a cell into a matrix  
%   called SegResultsFinal and performs a statistical analysis on the four 
%   segment populations: free, confined, directed and immobile diffusion  
%   (i.e. mean lifetime, diffusion coefficient, etc.). 
%   IMPORTANT: All data in the matrix are given in units of frames and 
%   pixels. The user needs to convert the data into appropriate units at a 
%   later time.
%
% Syntax:  
%   [SegResultsFinal,SegStats]=getSegResults_2(results)
% 
% Inputs:
%   results - file results.mat, which is the output of the DC-MSS function
%   'basicTransientDiffusionAnalysisv1.m' [1]
%
%   [1] Vega, A. R., et al. (2018). "Multistep Track Segmentation and 
%   Motion Classification for Transient Mobility Analysis." Biophysical 
%   Journal 114(5): 1018-1025.
% 
% Outputs:
%   1. SegResultsFinal.mat (saved in each selected file directory)            
%      Content of SegResultsFinal is same output as described in
%      basicTransientDiffusionAnalysisv1 (DC-MSS):
%      Column (1) Start frame of segment.
%      Column (2) End frame of segment.
%      Column (3) Classification of segment: 0 immobile, 1 confined, 2 free  
%      diffusion, 3 super diffusion
%      Column (4) MSS slope resulting in classification.
%      Column (5-11) Generalized diffusion coefficients for moment orders 0-6.
%      Column (12-18) Scaling power for moment orders 0-6.
%      Column (19) Normal diffusion coefficient (from the MSD).
%      Column (20) Confinement radius or localization precision, if subpart 
%      is classified as confined or immobile (NaN otherwise).
%      Column (21/22/-) Center of subpart, if subpart is classified as
%      confined (NaN otherwise).
%      Column (23) Trace count is the trace ID number
%      Column (24) track length of each segments in frames
%
%   2. SegStats.mat. Results of the statistical analysis. (saved in each 
%      selected file directory)
%      Row [1-4] the 4 mnotion classes (see Column [1])
%      Column [1]: motion class -> 0 immobile, 1 confined, 2 free diffusion, 
%      3 super diffusion
%      Column [2] nummber of segments in each class
%      Column [3] Facktion of segments out of total number of segments
%      Column [4] Total time spend in each class
%      Column [5] Total fraction of time spend in each class, normalized to
%      total
%      Column [6] Median lifetime
%      Column [7] Mean lifetime 
%      Column [8] Mean diffusion coef
%      Column [9] Std diffusion coef
%      Column [10] Median diffusion coef
%      Column [11] Mean Radii
%      Column [12] std Radii
% 
% 
% See also: 
%   scriptGetSegResults.m
%
% Author: 
%   - Signe Mathiasen May 2018 (origininal function getSeResults.m)
%   - P.G. Aug 2018
%     modified originial function getSeResults.m to process cases where 
%     1. the results matrix is empty: results = [] 
%     2. all track segments <20 frames will be deleted in SegResultsFinal
%   - P.G. March 2019 
%     segResultsFinal contains now all tracks including those with a 
%     length < 20 frames. But the calculation in SegStats does not include 
%     these tracks (tracks < 20 frames) and therefore all results will be 
%     the same as in the previous version.
%     Keeping all the tracks has the advantage the the ID numbers in 
%     column 23 of segResultsFinal are complete. The Ids can be used in
%     trace filters which are based on the diffusion properties. 
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------
%% Load Data
iTrack=numel(results);
SegResults=[];
TraceCount=[];
counter=0;

for k=1:iTrack
        
        counter = counter+1;
        NumSeg  = numel(results(k).segmentClass.momentScalingSpectrum(:,1));            
        CountMatrix = ones(NumSeg,1).*counter;
        TraceCount  = [TraceCount;CountMatrix]; %#ok<AGROW>
        SegResults  = [SegResults;results(k).segmentClass.momentScalingSpectrum]; %#ok<AGROW>

end

%% Creating the final results matrix keeping NaN for the full data

if ~isempty(SegResults)
    NotNans = find(~isnan(SegResults(:,3)));
%%% segLength = (SegResults(NotNans,2)-SegResults(NotNans,1))+1; %add one frame? To be discussed?!
%%% SegResultsFinal = horzcat(SegResults(NotNans,:),TraceCount(NotNans),segLength);

    segLength = (SegResults(:,2)-SegResults(:,1))+1; %add one frame? To be discussed?!
    SegResultsFinal = horzcat(SegResults(:,:),TraceCount(:),segLength);

else
    SegResultsFinal = [];
end

if isempty(SegResultsFinal)
    SegResultsFinal=[];
    SegStats=[0,0,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
              1,0,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
              2,0,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
              3,0,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
   return
end

% Content of SegResultsFinal
% [1]-[22] is same output as described in basicTransientDiffusionAnalysisv1
% [23] Trace count is the trace ID number
% [24] track length of each segments in frames


%% At the end concatenate all diffusion parameters into one final variable
%%% segCount    = numel(SegResultsFinal(:,1));
segCount    = size(NotNans,1);
totalLength = sum(SegResultsFinal(NotNans,24));

% Calculate statistical parameters and save them in SegStats. 
idxClass0 = SegResultsFinal(:,3)==0;
idxClass1 = SegResultsFinal(:,3)==1;
idxClass2 = SegResultsFinal(:,3)==2;
idxClass3 = SegResultsFinal(:,3)==3;

% [1] class, 0 immobile, 1 confined, 2 free diffusion, 3 super diffusion
Class0 = SegResultsFinal(idxClass0,:);
Class1 = SegResultsFinal(idxClass1,:);
Class2 = SegResultsFinal(idxClass2,:);
Class3 = SegResultsFinal(idxClass3,:);


% [2] Nummber of segments in each class
ClassNum=[numel(Class0(:,1));...
          numel(Class1(:,1));...
          numel(Class2(:,1));...
          numel(Class3(:,1))];
% [3] Faction of segments out of total number of segments
ClassFrac=ClassNum./segCount;

% [4] Total time spend in each class
ClassTime=[sum(Class0(:,24));...
           sum(Class1(:,24));...
           sum(Class2(:,24));...
           sum(Class3(:,24))];
% [5] Total fraction of time spend in each class, normalized to total
ClassTimeNorm=ClassTime./totalLength;

% Lifetimes in frames mean and median
l0=median(sort(Class0(:,24)));
l1=median(sort(Class1(:,24)));
l2=median(sort(Class2(:,24)));
l3=median(sort(Class3(:,24))); 

% [6] Median lifetime
LTmedian=[l0;l1;l2;l3];

% [7] Mean lifetime 
LTmean=[mean(Class0(:,24));...
        mean(Class1(:,24));...
        mean(Class2(:,24));...
        mean(Class3(:,24))];

%Diffusion coefficients mean, std and median
% [8] Mean diffusion coef
DiffCoMean=[mean(Class0(:,19));...
            mean(Class1(:,19));...
            mean(Class2(:,19));...
            mean(Class3(:,19))];
% [9] Std diffusion coef
DiffCoStd= [std(Class0(:,19));...
            std(Class1(:,19));...
            std(Class2(:,19));...
            std(Class3(:,19))];
% [10] Median diffusion coef
m0=median(sort(Class0(:,19)));
m1=median(sort(Class1(:,19)));
m2=median(sort(Class2(:,19)));
m3=median(sort(Class3(:,19)));
DiffCoMedian=[m0;m1;m2;m3];

% Confinement radii
ndRadii = Class0(:,20);
cdRadii = Class1(:,20);
% [11] Mean Radii
RadiiMean=[mean(ndRadii);...
           mean(cdRadii);...
           NaN;...
           NaN];
% [12] std Radii
RadiiStd=[std(ndRadii);...
          std(cdRadii);...
          NaN;...
          NaN];

% Generate output variable SegStat      
Class=[0;1;2;3];
SegStats=horzcat(Class,...          % [1] class, 0 immobile, 1 confined, 2 free diffusion, 3 super diffusion
                 ClassNum,...       % [2] nummber of segments in each class
                 ClassFrac,...      % [3] Facktion of segments out of total number of segments
                 ClassTime,...      % [4] Total time spend in each class
                 ClassTimeNorm,...  % [5] Total fraction of time spend in each class, normalized to total
                 LTmedian,...       % [6] Median lifetime
                 LTmean,...         % [7] Mean lifetime 
                 DiffCoMean,...     % [8] Mean diffusion coef
                 DiffCoStd,...      % [9] Std diffusion coef
                 DiffCoMedian,...   % [10] Median diffusion coef
                 RadiiMean,...      % [11] Mean Radii
                 RadiiStd);         % [12] std Radii

%% End



