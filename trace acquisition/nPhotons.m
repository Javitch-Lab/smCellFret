function [nMol, maxIntensity, selList] = nPhotons(filename,QE,...
                                          xTranslation, yTranslation)
%--------------------------------------------------------------------------
%
% nPhotons.m:
%   The function measures the spot intensity at each track position in
%   units of photons per frame.
% 
% Description:
%   The fluorescence intensity of a single molecule spot and the 
%   corresponding background value are calculated over the 5x5 central  
%   region and the adjacent background region as the sum of discrete pixel  
%   intensities. The number of pixels in the 5x5 central region is  
%   consistent with the number of pixels in the adjacent background region. 
% 
% Syntax:  
%   [nMol, maxIntensity, selList] = nPhotons(filename,QE,...
%                                   xTranslation, yTranslation)
% Inputs:
%   1. filename     - traces file with a baseline (all files with the 
%                     extension '_bln.traces')
%   2. QE           - quantum efficiency of the EMCCD camera at the emisson 
%                     wavelength
%   3. xTranslation - correction factor for a coordinate shift in x
%   4. yTranslation - correction factor for a coordinate shift in y
%   5. widthPSF     - the width of the pixel area to capture 80% of the
%                     light from a single molecule PSF. This value is 
%                     defined in 'smCellConstants.m'.
%
% Outputs:
%   1. nMol         - number of molecules that were successfully processed
%   2. maxIntensity - highest photon count detected in the dataset (scalar) 
%   3. selList      - list of indices (not trace Id's) of successfully 
%                     processed molecules (molecules at boundaries are 
%                     removed).
% 
% Other m-files required: 
%   Subfunctions: loadTracesCell.m, saveTracesCell.m
% 
% See also: 
%   scriptGetFretTraces.m, 
%
% Authors: 
%   -P.G. Aug 2011
%   -P.G. Jan 2015 (added calc. of the S/Bgd ratio.)
%   -P.G. May 2017 a)use round instead of floor
%                    b)divde photonsBgd with QE 
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% --- Use as script for test purposes --------------------------------------
% % % %  clear;
% % % %  filename = 'F:\#Data\161114_MH_NT-sf-mGlu2_LiveFRET_VariousPowers_PCAPCD\333nMCy3AC_666nMCy5AC_100ng Tet_PCAPCD-50\3mW\Tracking_FinalOpt Con_WA 82617\MEFfilt_PGsettings_Snr-Crt_Dinfo\#09Ch1\singleTracesAcc_bln.traces';
% % % %  QE = 0.97;
% % % %  xTranslation=0; 
% % % %  yTranslation=0;
%% --------------------------------------------------------------------------

% Set image directory
[path,name,~]=fileparts(filename);
imageDir   = [path filesep 'ImageData\'];

% Load Constants
constants = smCellConstants('');
widthPsf  = constants.widthPSF;
photonConversion = constants.photonConvFactorDelta;

if nargin<2
  QE=0.95; % set Camera Quantum efficiency if not specified
  
elseif nargin<3
  xTranslation=0;  % correction of a xCoordinate shift 
  yTranslation=0;  % correction of a yCoordinate shift
  
end

% Read first image and get image size
list=dir([imageDir '*01.tif']);
imageTmp = imread([imageDir list(1).name]);
[imHeigth,imWidth] = size(imageTmp);
clear imageTmp

%--- Load the trace data
data   = loadTracesCell(filename);
x      = data.x-xTranslation;
y      = data.y-yTranslation;

[nTraces,nFrames]=size(x);
ids = cell(nTraces,1);
for i=1:nTraces
    ids{i}=data.traceMetadata(i).ids;
end

% Allocate memory
% intPsf     = zeros(widthPsf,widthPsf)*NaN;
photons    = zeros(nTraces,nFrames);
photonsBgd = zeros(nTraces,nFrames);
snrPSF     = zeros(nTraces,nFrames);
imageStk   = zeros(imHeigth,imWidth,nFrames);

% Load Image Stk
for i=1:nFrames
    imageStk(:,:,i)=double(imread([imageDir 'fov1_' sprintf('%05u',i) '.tif']));
end

% Settings for the integration over a pixel area 
if     widthPsf==3; step=2; %for 3x3 area 
elseif widthPsf==5; step=3; %for 5x5 area 
elseif widthPsf==7; step=4; %for 7x7 area 
end 
%% --- Calculate number of photons: Loop over every detected molecule time 
%  --- trace

h = waitbar(0,'Please wait...');

selList=[];interrupt=0;
for n=1:nTraces
    xCor=round(x(n,:));
    yCor=round(y(n,:));
    
    % only analyze the part that has non-zero entries
    idx1=find(xCor(1,:)>0,1,'first'); 
    idx2=find(xCor(1,:)>0,1,'last');
    
    % Get the intensity values of a tracked molecule n in time range defined
    % by the index variables: idx1, idx2.
    for t=idx1:idx2
        if xCor(t)==0 && t>1
            photons(n,t)=photons(n,t-1);
            continue; 
        end
        
        %--- Measure the intensity over the pixel area: for test purposes
        %--- use the function 'sumPixelArea.m'
        intPsf=NaN(widthPsf,widthPsf);
        for i=1:widthPsf
            for j=1:widthPsf
                col=yCor(1,t)+i-step;
                row=xCor(1,t)+j-step;
                if col>0 && row >0
                    if col<=imHeigth && row <=imWidth
                        intPsf(i,j) = imageStk(col,row,t);
                    end
                end
            end
        end

        %--- Integrate Intensities in the central region and in the
        %--- boundary region
%         if sum(sum(isnan(intPsf)))>0
%             interrupt=1;
%             disp('found NaN');
%             break
%         else
            %--- Intgrate intensity values over the central pixel area of 
            %--- the psf 
            photons(n,t)    = sum(sum(intPsf(2:widthPsf-1,2:widthPsf-1)))/QE;
            photons(n,t)    = photons(n,t)./photonConversion;

            %--- Integrate background intensity values over the  
            %--- boundary region of the psf
            bgdPsfTmp = [intPsf(1,1:widthPsf);...           % upper row
                intPsf(1:widthPsf,widthPsf)'; ...        % rigth cloumn
                fliplr(intPsf(widthPsf,1:widthPsf));...  % bottom row
                fliplr(intPsf(1:widthPsf,1)')];          % left cloumn
            
            %--- Delete the 1st coloum of the matrix (remove elements that
            %--- were counted twice at the corners of the matrix)             
            bgdPsfTmp(:,1) = [];
            bgdPsfTmp = reshape(bgdPsfTmp,1,numel(bgdPsfTmp));

            %--- Elements of the bgd matrix need to match the number of
            %--- matrix elements of the intPsf. 
            %   size        size of                   number
            % widthPSF    intPsf kernel         surrounding pixels
            %    3x3     1x1 (1 element)   8  elem (choose 1 elem randm)
            %    5x5     3x3 (9 elements)  16 elem (choose 9 elem randm)
            %    7x7     5x5 (25 elements) 24 elem (choose all elem +1 randm)
            
            if widthPsf == 3
                bgdPsfTmp=bgdPsfTmp(5:8);
                bgdPsf=mean(bgdPsfTmp);
                
            elseif widthPsf == 5
                bgdPsfTmp=[bgdPsfTmp(1:4), bgdPsfTmp(9:12)];
                bgdPsf=[bgdPsfTmp,mean(bgdPsfTmp)];
                
            elseif widthPsf == 7
                bgdPsf=[bgdPsfTmp, mean(bgdPsfTmp)];
                
            end
            
            photonsBgd(n,t)  = sum(bgdPsf)/QE;
            photonsBgd(n,t) = photonsBgd(n,t)./photonConversion;

            %--- test if the number of matrix elements are equal for the  
            %--- intPsf and bgdPsf before calculating the S/N ratio
            nIntPsf = numel(intPsf(2:widthPsf-1,2:widthPsf-1));
            nBgdPsf = numel(bgdPsf);
            if nIntPsf == nBgdPsf 
                snrPSF(n,t) = photons(n,t)/photonsBgd(n,t);
            end
%         end
    end
    if interrupt == 0
        selList=[selList n];
    else 
        interrupt = 0;
    end 
    waitbar(n / nTraces)
end
close(h) 
%%
% Calculate number of molecules
nMol=size(selList,2);

% Calculate max photon count in the data set 
maxIntensity=max(max(photons));

% Create an output filename
outfile = [path filesep name '.traces'];
outfile = strrep(outfile,'bln.traces','pht.traces');

% Save traces in the specified in the selection list with new intensity
data.channelNames = {'x','y','int','snr','locBgd'};
data.x   = x(selList',:);
data.y   = y(selList',:);
data.int = photons(selList',:);
data.snr = snrPSF(selList',:);
data.locBgd = photonsBgd(selList',:);
for i=1:numel(selList)
    data.traceMetadata(i).ids =ids{selList(i)};
end
saveTracesCell(outfile,'traces',data);


