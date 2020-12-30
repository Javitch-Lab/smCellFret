%--------------------------------------------------------------------------
%
% scriptContourPlot_flt:
%   Batch processing script that generates data for a 1d FRET histogram and 
%   a 2d FRET pseudo color contour plot. 
%
% Description:
%   Running the script will open the file explorer and prompt the user to 
%   select several *.traces files (located in #Ch2 directories). The script 
%   uses the SPARTAN [1] function makecplot.m and generates two ASCII data
%   files which will be saved under the same file name as the input file 
%   with the additional extension *._2dHst.txt and *._1dHst.txt. ASCII 
%   text files can be imported into a data analysis software (e.g. Origin) 
%   to generate 2dim-FRET contour plots and 1dim-FRET histogram 
%   distributions.
% 
% Syntax:  
%   scriptContourPlot_flt
% 
% Inputs:
%   1. *_sync*.traces files - Fret-traces files that have been 
%      post-synchronized using 'scriptPostSyncTraces.m'. The script prompts 
%      the user to select multiple traces-files.
%   2. 'integrWindwow' parameter (see code section 'Set initial parameters') 
%      - An integer value n that defines the time interval [0, n] (unit of  
%      n is frames) of the 2d-FRET histogram (pseudo-color contour plot).   
%      The 1d-FRET histogram is calculated by integration over the interval  
%      [0, n] of the 2d-FRET histogram.  
%   3. applyTotalIntensityFilter is a binary variable (see code section 'Set
%      initial parameters') - The variable should have the value 'false'. 
%
%      Todo: Remove the total intensity filter option from the script. 
%      Filtering by total intensity is handled by a separate script 
%      (scriptTotalFilter_traces.m).           
% 
% Outputs:
%   1. ASCII data txt file with extensions *_1dHst.txt - The file contains 
%      the data for the 1d-Fret histogram. 
%   2. ASCII data txt file with extensions *_2dHst.txt - The file contains 
%      the data for the 2d-Fret contour.
%   NOTE: Import the data, for example, into the "Origin" graphics software 
%   to generate the histogram and contour plots.  
% 
% Other m-files required: 
%   Subfunctions: SPARTAN functions makecplot.m and cascadeConstants.m [1] 
% 
% See also: 
%   scriptContourPlot.m, scriptPostSyncTraces.m , scriptTotalFilter_traces.m
%   
% References 
%   [1] Juette, M. F., et al. (2016). "Single-molecule imaging of 
%   non-equilibrium molecular ensembles on the millisecond timescale." 
%   Nature Methods 13(4): 341-344.
%
% Author:
%   - P.G. 2018, scriptContourPlot
%   - P.G. 2019, added batch processing and a total intensity filter option
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Set initial parameters
clear;

integrWindwow = 125; %unit is frames
applyTotalIntensityFilter = false; 


%% -------------------------------------------------------------------------
% Get files from user
% -------------------------------------------------------------------------

% Get FILES
tfiles = getFiles('*.traces','Select *.traces files');

% Check if files array is empty
if any(cellfun('isempty',tfiles))== true; return;end

% Show selection
for z=1:numel(tfiles)
    [path,name,ext]=fileparts(tfiles{z});
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
nCells = numel(tfiles);
data    =struct([]);
h = waitbar(0,'Load Data ...');
for j=1:nCells
    clear tmp;
    waitbar(j / nCells)
    
    % load traces file
    tracesFile = tfiles{j};
    tmpData=loadTraces(tracesFile);
    totalInt = tmpData.total;
    [nTraces,nFrames]=size(totalInt);
    
    % save filename
    data(j).file = tracesFile;
    
    if applyTotalIntensityFilter==true

        % logical Idl Array 
        idl = false(nTraces,nFrames);
        traceLenDon = [tmpData.traceMetadata.traceLenDon].';
        for i=1:nTraces
            %dwtTmp = dwells{2,i}(2);
            dwtTmp = traceLenDon(i); %unit is frames
            idl(i,1:dwtTmp)=true;
        end

        % Calc Total intensity where trace is idealized 
        totalInt_idl2D = totalInt.*idl;
        totalInt_idl2D(totalInt_idl2D==false)=NaN;
        meanTotalPerTrace = nanmean(totalInt_idl2D,2);
        data(j).total = meanTotalPerTrace;

        % Get Filter Idx
        idx = meanTotalPerTrace > lwLim & meanTotalPerTrace < upLim;
        data(j).flt = idx;
        data(j).nMol = sum(idx);

        % Apply Filter to data
        tracesData  = applyFilterToData(tmpData,idx);
        data(j).fltTraces = tracesData;

        % Save Traces to File
        [pathName, fileName, ext]= fileparts(tracesFile);
        newExt = strrep(ext,'.traces','_flt.traces');
        saveTraces([pathName filesep fileName newExt],tracesData)
    
    else
        
        % Get Filename
        [pathName, fileName, ext]= fileparts(tracesFile);
        
        data(j).flt = [];
        data(j).nMol = size(tmpData.donor,1);
        data(j).fltTraces = tmpData;
        tracesData  = tmpData;
    end
    
    %---------------- Make Histogram -------------------------------------- 
    try
        constants = cascadeConstants;% cascade constants from spartan v5.3.1 
        options   = constants.defaultMakeplotsOptions;

        % normalize contour plot
        contourLength = integrWindwow; %unit is frames

        % normalize to the number of traces
        options.cplot_normalize_to_max = false;
        cHst2D = makecplot(tracesData,options); % makecplot from spartan v5.3.1 
                                               % cHst2D is normalized by number of traces
                                               % the sum of every colum in cHst2D=1 

        % Calculate Contour Profile and nTraces 
        % In this normalization the sum of cHst1D gives 1
        popHist_contourLength = nansum(cHst2D(2:end,2:contourLength+1),2);
        cHst1D  = [ cHst2D(2:end,1) 100*(popHist_contourLength/contourLength)];

        data(j).cHst2D = cHst2D;
        data(j).cHst1D = cHst1D;

        % % % normalize to the number of frames
        % % % options.cplot_normalize_to_max = 4000;
        % % % cHst2D = makecplot(filename,options); %makecplot from spartan v5.3.1 
        % % % popHist_nTraces = sum(cHst2D(2:end,2:end),2);
        % % % cHst1D  = [ cHst2D(2:end,1) popHist_nTraces];

        % Cntrl Step: Integrate the Histogram to get the number of traces or ...
        % frames (contour length)
        sum1DHst = sum(cHst1D(1:end,2));
        disp(['sum1DHst :' num2str(sum1DHst)]);
        
        % Save Hist Data
        newExt = strrep(ext,'.traces','_2dHst.txt');
        writematrix(cHst2D,[pathName filesep fileName newExt]);
        newExt = strrep(ext,'.traces','_1dHst.txt');
        writematrix(cHst1D,[pathName filesep fileName newExt]);
        
    catch
        disp(['data set ' num2str(j) ' is empty'])
    end
    %----------------------------------------------------------------------
    
    
end
close(h)

%% Disp Number of Molecules
nMol = sum([data.nMol].');
disp(['nMol = ' num2str(nMol)]);


%% =======================================================================%
%                                                                         %
%                        Local Callback Functions                         %
%                                                                         %
%=========================================================================%


%% filter matrix array

function fltTraces=applyFilterToData(fretData,idxFilter)

    
    % Update data
    [nTraces,nFrames]=size(fretData.donor); 
    
    % Create traces object
    fltTraces = TracesFret(nTraces,nFrames);
    
    % Metadata Field Names
    metaFN = fieldnames(fretData.traceMetadata);
             
    % Update data
    fltTraces.donor            = fretData.donor(idxFilter,:);
    fltTraces.acceptor         = fretData.acceptor(idxFilter,:);
    fltTraces.fret             = fretData.fret(idxFilter,:);
    fltTraces.time             = fretData.time;
    
    fltTraces.traceMetadata    = [];

    % Update traceMetadata for each trace 
    j=0;
    for i=1:nTraces

        % Apply Filter
        if idxFilter(i)==0; continue;end
        j=j+1;

        %Update Metadata
        for k=1:numel(metaFN)
            fn=metaFN{k};
            if isfield(fretData.traceMetadata,fn)
                fltTraces.traceMetadata(j).(fn) =...
                fretData.traceMetadata(i).(fn);
            end
        end

    end % End Update Metadata
end
%% End































