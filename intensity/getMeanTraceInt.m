%--------------------------------------------------------------------------
%
% getMeanTraceInt:
%   Script that calculates the average total trace intensity and the 
%   average acceptor trace intensity using manually measured trace 
%   lifetimes. 
% 
% Syntax:  
%   getMeanTraceInt
% 
% Inputs:
%   1. The script propmts the user to select one or multiple 
%      fretTraces*.mat files (e.g. fretTracesImob_flt_allFret.mat) 
%   2. The script expects a similar file fretTracesImob_flt_allFret_lt.txt
%      in the same directory. This text file must contain the 
%      comma-separated lifetime data of the total intensity time trace and 
%      the acceptor time trace. Note: these are manually measured lifetime 
%      data.
% 
% Outputs:
%   dataAll - Workspace variable, with the following four columns:
%    Col 1: ltTotal, lifetime of the total intensity trace
%    Col 2: meanTotalInt, mean total intensity
%    Col 3: ltAcc, lifetime of the acceptor time trace 
%    Col 4:  meanAccIn, mean acceptor intensity 
%   
% Example: 
%  dataAll Matrix: 
%    ltTotal  | meanTotalInt     | ltAcc   | meanAccInt 
%    [23,      727.029563962151,    23,     183.435444160910;
%     60,      492.379196644783,    31,     152.179032819194;
%    329,      583.175760631459,   185,     267.527664174673;
%   1336,      700.244737610564,  1171,     183.220958046794]
% 
% Authors: 
%   - P.G. 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
% Get files from user
% -------------------------------------------------------------------------

% Get FILES
clear;
matfiles = getFiles('*.mat','Select *.mat files');

% Check if files array is empty
if any(cellfun('isempty',matfiles))== true; return;end

% Show selection
for z=1:numel(matfiles)
    [path,name,ext]=fileparts(matfiles{z});
    disp(path);
    disp([name ext]);
end

%% -------------------------------------------------------------------------
% Filter data
% -------------------------------------------------------------------------

% Total Intensity Filter Parameters 
fixedMean = 457.91;
fixedStdv = 218.144/2;
lwLim = fixedMean-2*fixedStdv;
upLim = fixedMean+3*fixedStdv;

% 
nCells = numel(matfiles);
data    =struct([]);
h = waitbar(0,'Load Data ...');

for j=1:nCells 
    clear tmp;
    waitbar(j / nCells)
    
    % load traces file
    tracesFile = matfiles{j};
    manEstimated_lt = strrep(matfiles{j},'.mat','_lt.txt');
    load(tracesFile);

    % Load manually modified LT (b/c of segementation in immobile traces)
    estLt=readmatrix(manEstimated_lt);
    totalInt = fretTraces.Fret.total;
    accInt   = fretTraces.Ch1.int;
    [nTraces,nFrames]=size(totalInt);
    
    % save filename
    data(j).file = tracesFile;

    % logical Idl Array 
     
    idl = fretTraces.Fret.idlTotal; 
    idlDon = false(nTraces,nFrames);
    idlAcc = false(nTraces,nFrames);
    for k=1:nTraces
        
        % Donor Idl
        sotDon=fretTraces.Ch2.traceMetadata(k).startOfTrace;
        eotDon=sotDon+estLt(k,1)-1;
%         disp([sotDon, eotDon]);
        idlDon(k,sotDon:eotDon)=true;
        
        % Acceptor Idl
        sotAcc=fretTraces.Ch1.traceMetadata(k).startOfTrace;
        eotAcc=sotAcc+estLt(k,2)-1;
%         disp([sotAcc, eotAcc]);
        idlAcc(k,sotAcc:eotAcc)=true;
    end
    
    % Calc Total intensity where trace is idealized 
    totalInt_idl2D = totalInt.*idlDon;
    totalInt_idl2D(totalInt_idl2D==false)=NaN;
    meanTotalPerTrace = nanmean(totalInt_idl2D,2);
    data(j).total = meanTotalPerTrace;
    
    % Calc Total intensity where trace is idealized 
    accInt_idl2D = accInt.*idlAcc;
    accInt_idl2D(accInt_idl2D==false)=NaN;
    meanAccPerTrace = nanmean(accInt_idl2D,2);
    data(j).acc = meanAccPerTrace;
    
    %output
    dataAll=[estLt(:,1),data.total,estLt(:,2),data.acc];
    
%     % Get Filter Idx
%     idx = meanTotalPerTrace > lwLim & meanTotalPerTrace < upLim;
%     data(j).flt = idx;
%     data(j).nMol = sum(idx);
    
%     % Apply Filter to data
%     fretTraces  = applyFilterToData(fretTraces,idx);
%     data(j).fltTraces = fretTraces;
%     
%     % Save Traces to File
%     [pathName, fileName, ext]= fileparts(tracesFile);
%     newExt = strrep(ext,'.mat','_flt.mat');
%     save([pathName filesep fileName newExt],'fretTraces')
%     
    
end
close(h)



