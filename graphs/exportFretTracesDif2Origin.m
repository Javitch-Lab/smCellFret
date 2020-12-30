%--------------------------------------------------------------------------
%
% exportFretTracesDif2Origin.m:
%   Export selected trace data from fretTracesDif.mat to 'Origin' 
% 
% Description:
%   The script exports user-selected data of the fretTracesDif.mat file 
%   to the 'Origin' graphics software to generate plots of spatial 
%   Donor/Acceptor (D/A) tracks, D/A diffusion tracks, D/A intensity traces 
%   and FRET traces for publication. D/A intensity traces and FRET traces 
%   are corrected for crosstalk (alpha), direct excitation (delta) and 
%   brightness (gamma).
%
% Syntax:  
%   exportFretTracesDif2Origin
% 
% Inputs:
%   fretTracesDif*.mat - The script will prompt the user to select  
%   e.g. the filefretTracesDif_best.mat, which contains only a subset of  
%   data extracted from the original file fretTracesDif.mat. Use 
%   cellFRETViewtraces to select the subset of data which need to be 
%   plotted for publication (see code section 'Select Trace in 
%   cellFretViewtraces') 
%
%   Parameters - see code section 'Set Parameters'
%   1. Update EMCCD camera pixel size
%   2. Update FRET correction factors. 
%   3. Set the variables doAlphaCorr, doDeltaCorr and doGammaCorr to either 
%   true or false to enable or disable the particular data correction method. 
%
% Outputs:
%   1. track - Workspace variable with the spatial coordinates of the D/A track.
%   2. traces - Workspace variable with the corrected D/A and FRET time traces
%   3. diffusion - Workspace variable with the coordinates of the four 
%   diffusion modes immobile, confined, free and directed. 
% 
% Other m-files required: 
%   Subfunctions: alphaCorrSingleTrace, deltaCorrSingleTrace and 
%   gammaCorrSingleTrace 
% 
% Author: 
%   - P.G. Apr 2019
%   - P.G. Sep 2020, added alpha, delta, gamma correction for FRET
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------
%% Select Trace in cellFretViewtraces
clear;

% 1. load fretTracesDif.mat file in cellFretViewtraces 
% 2. Select the trace that needs to be exported
% 3. Check Bgd correction
% 3. Correct crosstalk in in vietraces OR use the crosstalk correction in 
%    in the script: Set doAlphaCorr = true to enable crosstalk correction
%    in the script. 
% 4. Save the trace e.g as fretTracesDif_bestFret.mat
% 5. Enter the molecule number of the trace saved in fretTracesDif_bestFret.mat 
%    that needs to be exported. If fretTracesDif_bestFret.mat has only one 
%    trace use mol = 1. 

mol = 1;

%% Set Parameters
pixelSize = 0.16; % 0.160 um/pixel

% Crosstalk Factor
alpha = 0.09417; doAlphaCorr = true;

% Direct Excitation Factor
delta = 0.056;   doDeltaCorr = true;

% Gamma Correction Factor
gamma = 0.83; doGammaCorr = true;    

%% Load file

% Get fretTracesDif_bestFret.mat file 
[file,path] = uigetfile;
if isequal(file,0); return; end

% load selection
disp(path);
disp(file);
fretFile = [path filesep file];
load(fretFile);

%% ------------------------------------------------------------------------
%
%                    Export tracks (unit is um)
%
%--------------------------------------------------------------------------

track=struct('xAcc',[],...
             'yAcc',[],...
             'xDon',[],...
             'yDon',[]);

%--- Acceptor
startOfTrace.Ch1 = fretTraces.Ch1.traceMetadata(mol).startOfTrace;
endOfTrace.Ch1   = fretTraces.Ch1.traceMetadata(mol).endOfTrace;
accTracking      = startOfTrace.Ch1:endOfTrace.Ch1;

for i=1:numel(accTracking )
    track(i).xAcc = (fretTraces.Ch1.xCorr(mol,accTracking (i)))*pixelSize;
    track(i).yAcc = (fretTraces.Ch1.yCorr(mol,accTracking (i)))*pixelSize;
end

%--- Donor
startOfTrace.Ch2 = fretTraces.Ch2.traceMetadata(mol).startOfTrace;
endOfTrace.Ch2   = fretTraces.Ch2.traceMetadata(mol).endOfTrace;
donTracking      = startOfTrace.Ch2:endOfTrace.Ch2;

for i=1:numel(donTracking)
    track(i).xDon = (fretTraces.Ch2.x(mol,donTracking(i)))*pixelSize;
    track(i).yDon = (fretTraces.Ch2.y(mol,donTracking(i)))*pixelSize;
end



%% ------------------------------------------------------------------------
%
%                          Export traces
%
%--------------------------------------------------------------------------

trace=struct('dtAcc',[],...
             'intAcc',[],...
             'intDon',[],...
             'total',[],...
             'idl',[],...
             'fret',[],...
             'barAcc',[],...
             'barDon',[]);

%--- Acceptor
lenBackset.Ch1   = fretTraces.Ch1.traceMetadata(mol).lenBackset;
lenBaseline.Ch1  = fretTraces.Ch1.traceMetadata(mol).lenBaseline;
traceStart.Ch1   = startOfTrace.Ch1 - lenBackset.Ch1;
traceEnd.Ch1     = endOfTrace.Ch1 + lenBaseline.Ch1;
totalAccRange         = traceStart.Ch1:traceEnd.Ch1;

for j=1:numel(totalAccRange)
    trace(j).dtAcc  = fretTraces.Ch1.time(totalAccRange(j));
    trace(j).intAcc = fretTraces.Ch1.int(mol,totalAccRange(j)); 
    
    if totalAccRange(j) >= startOfTrace.Ch1 && ...
            totalAccRange(j) <= endOfTrace.Ch1 
        trace(j).barAcc = -100;
    else
        trace(j).barAcc = NaN;
    end

end

% --- Donor
lenBackset.Ch2   = fretTraces.Ch1.traceMetadata(mol).lenBackset;
lenBaseline.Ch2  = fretTraces.Ch2.traceMetadata(mol).lenBaseline;
traceStart.Ch2   = startOfTrace.Ch2 - lenBackset.Ch2;
traceEnd.Ch2     = endOfTrace.Ch2 + lenBaseline.Ch2;
totalDonRange    = traceStart.Ch2:traceEnd.Ch2;
donRange         = startOfTrace.Ch2 : endOfTrace.Ch2;

for j=1:numel(totalDonRange)
    
    trace(j).intDon = fretTraces.Ch2.int(mol,totalDonRange(j));
    trace(j).total  = fretTraces.Fret.total(mol,totalDonRange(j));
    trace(j).idl    = fretTraces.Fret.idlTotal(mol,totalDonRange(j));
    fret = (trace(j).intAcc./trace(j).total)*trace(j).idl; 

    if ~isempty(trace(j).intAcc)
        if isnan(fret) && totalDonRange(j) < startOfTrace.Ch2 
            fret=0;
        elseif isnan(fret) && totalDonRange(j) > endOfTrace.Ch2 
            fret=0;
        end
    else
        fret=NaN;
    end
    
    trace(j).fret   = fret; 
    if totalDonRange(j) >= startOfTrace.Ch2  && ...
            totalDonRange(j) <= endOfTrace.Ch2 
        trace(j).barDon = -150;
    else
         trace(j).barDon = NaN;
    end
   
end

%% data for correction
time     = cell2mat({trace.dtAcc}.');
intAcc   = cell2mat({trace.intAcc}.'); % total acc range
intDon   = cell2mat({trace.intDon}.'); % total donor range
idlTotal = cell2mat({trace.idl}.');    % total donor range

accRange = cell2mat({trace.barAcc}.')<0; % idl acc range in workspace var traces
donRange = cell2mat({trace.barDon}.')<0; % idl don range in workspace var traces


%% crosstalk correction 
if doAlphaCorr == true
    [intAcc,intDon,acorr,fret] = alphaCorrSingleTrace(...
        alpha, time, intAcc, intDon, totalAccRange, idlTotal);

    % update workspace var traces
    for j=1:numel(intAcc)
        trace(j).intAcc = intAcc(j);
        trace(j).fret   = fret(j);
        trace(j).total  = intAcc(j)+intDon(j);
    end
end

%% direct excitation correction
if doDeltaCorr == true
    % Initialize the idealization matrix with logical zeros
    [nTraces, nFrames] = size(fretTraces.Ch1.x);
    idl = false(1,nFrames);

    % set idl-matrix TRUE over the donor range
    idl(startOfTrace.Ch2:endOfTrace.Ch2) = true; % idl don range

    % Calc mean total intensity per trace over donor range 
    % Subsitute zeros with NaN to get correct mean total intensities 
    intTotal = cell2mat({trace.total}.'); % use updated total for calculation
    total_idl2D = intTotal.*donRange;
    total_idl2D(total_idl2D==false)=NaN;
    meanTotal = nanmean(total_idl2D); 

    % do the correction
    [intAcc, intDon, dcorr, fret] = deltaCorrSingleTrace(...
        delta, time, intAcc, intDon, accRange, totalAccRange,idlTotal, meanTotal );

    % update workspace var traces
    for j=1:numel(intAcc)
        trace(j).intAcc = intAcc(j);
        trace(j).fret   = fret(j);
        trace(j).total  = intAcc(j)+intDon(j);
    end
end

%% gamma correction
if doGammaCorr == true
    [intAcc, intDon, gcorr, fret] = gammaCorrSingleTrace(...
        gamma, time, intAcc, intDon, totalAccRange, idlTotal);

    % update workspace var traces
    for j=1:numel(intAcc)
        trace(j).intDon = intDon(j);
        trace(j).fret   = fret(j);
        trace(j).total  = intAcc(j)+intDon(j);
    end
end

%% ------------------------------------------------------------------------
%
%                          Export Diffusion
%
%--------------------------------------------------------------------------
         
diffusion=struct('xImb',[],...
                 'yImb',[],...
                 'xCnf',[],...
                 'yCnf',[],...
                 'xFre',[],...
                 'yFre',[],...
                 'xSpr',[],...
                 'ySpr',[]);


for k=1:numel(accTracking)
diffusion(k).xImb = fretTraces.Diff.Ch1.idlImb(mol,accTracking(k))*track(k).xAcc;
diffusion(k).yImb = fretTraces.Diff.Ch1.idlImb(mol,accTracking(k))*track(k).yAcc;

diffusion(k).xCnf = fretTraces.Diff.Ch1.idlCnf(mol,accTracking(k))*track(k).xAcc;
diffusion(k).yCnf = fretTraces.Diff.Ch1.idlCnf(mol,accTracking(k))*track(k).yAcc;

diffusion(k).xFre = fretTraces.Diff.Ch1.idlFre(mol,accTracking(k))*track(k).xAcc;
diffusion(k).yFre = fretTraces.Diff.Ch1.idlFre(mol,accTracking(k))*track(k).yAcc;

diffusion(k).xSpr = fretTraces.Diff.Ch1.idlSpr(mol,accTracking(k))*track(k).xAcc;
diffusion(k).ySpr = fretTraces.Diff.Ch1.idlSpr(mol,accTracking(k))*track(k).yAcc;


end
% Replace zeros with NaN
[diffusion([diffusion.xImb]==0).xImb]=deal(NaN);
[diffusion([diffusion.yImb]==0).yImb]=deal(NaN);

[diffusion([diffusion.xCnf]==0).xCnf]=deal(NaN);
[diffusion([diffusion.yCnf]==0).yCnf]=deal(NaN);

[diffusion([diffusion.xFre]==0).xFre]=deal(NaN);
[diffusion([diffusion.yFre]==0).yFre]=deal(NaN);

[diffusion([diffusion.xSpr]==0).xSpr]=deal(NaN);
[diffusion([diffusion.ySpr]==0).ySpr]=deal(NaN);



%% End
clearvars -except track trace diffusion
disp('finished');

%%
