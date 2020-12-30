function varargout = gammaCorrectCell(input, gammaFactor)
% gammacorrect scale acceptor intensity by estimated gamma values.
%
%   gammacorrect( FNAME ) Corrects for unequal sensitivity and/or quantum yield
%   ('apparent brightness') of donor and acceptor fluorophores in two-color
%   FRET traces by scaling the donor intensity. This scaling factor (gamma) is
%   estimated for each trace as the magnitude increase in donor intensity
%   divided by the magnitude decrease in acceptor intensity upon acceptor
%   photobleaching. A single fixed gamma value is used to make the correction,
%   obtained from the median of apparent trace gamma values, excluding outliers.
%   FNAME is the path of the .traces file to correct. The corrected traces are
%   automatically saved in a new file with the extension '_gcorr.traces'.
%
%     NOTE: for best results, use this function after crosstalk and other
%     corrections and trace selection in autotrace.
%   
%     See Ha et al (1999) PNAS 95, p.893 for original algorithm description.
%
%   gammacorrect( FILES ) corrects each file in the cell array FILES,
%   calculating a distinct gamma value for each file.
%
%   [DATA,GAMMA] = gammacorrect( FNAME ) returns the corrected data and
%   gamma value used for correction instead of saving to file.
%
%   DATA = gammacorrect( DATA ) corrects the TracesFret object DATA directly.
%   (NOTE: input is modified in place. To preserve the original, create a copy.)
%
%   [...] = crosstalkcorrect( ..., GAMMA ) uses a user-supplied gamma value
%   instead of calculating one.
% 
%   See also: calc_gamma, scaleacceptor, crosstalkcorrect, correctTraces.
%
%   Copyright 2014-2017 Cornell University All Rights Reserved.
%
% =========================================================================
%
%  gammaCorrectCell is a modified version of the SPARTAN function
%  gammacorrect
%
% =========================================================================
%
% gammaCorrectCell:
%   Performs a gamma correction
% 
% Description:
%   The function gammaCorrectCell.m scales the donor intensities by a 
%   correction factor gamma. Gamma correction is necessary to adjust for
%   differences in the donor - acceptor fluorescence quantum yields and
%   fluorescence detection efficiencies. 
% 
% Syntax:  
%   varargout = gammaCorrectCell(input, gammaFactor)
% 
% Inputs:
%   input - traces file with the extension *_a_d.traces
%   gammaFactor - empirically determined gamma value
%   
%   When executing the function, without input arguments the user is
%   prompted to: 
%   1. select a list of trace files that are already corrected for spectral
%      crosstalk and direct excitation (files with the extension *_a_d.traces). 
%   2. either confirm the default gamma factor or enter a new value. 
% 
% Outputs:
%   Corrected traces are saved in the same directory as the input file with 
%   the additional extension g for gamma (e.g. *sync_flt_a_d_g.traces). 
% 
% Other m-files required: 
%   Subfunctions: correctDonCell.m
% 
% See also: 
%   alphaCorrectCell.m,  deltaCorrectCell.m
%
% Author: 
%   P.G. Sep 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% Check range of input arguments
narginchk(0,2);
nargoutchk(0,2);

%% Get input file from user
if nargin<1 || isempty(input)
    % If no files given, obtain a list from the user.
    filter = {'*.*traces','Traces files (*.traces,*.rawtraces)';
              '*.*','All files (*.*)'};
    input = getFiles(filter);
    if isempty(input), return; end  %user hit cancel.
    
elseif iscell(input)
    if ~all(cellfun(@ischar,input))
        error('Input must be a Traces object, filename string, cell array of filenames');
    end
    if nargout>0
        error('Output only supported for scalar input')
    end
end

%% Set gamma factor by user
if nargin<2 || isempty(gammaFactor)
    
    gammaDefault = '0.83';
    gamma  = inputdlg({'Enter new Value:'},...           %prompt
                     'Set Global Direct Exitation Factor',... %dlg titel
                     [1 55],...                       %width of line
                     {gammaDefault});                 %default answer

    if isempty(gamma) || isnan(str2double(gamma{1,1}))
        % user hit cancel or input field was empty
        return
    else
        gammaFactor = str2double(gamma{1,1});
    end
    
end

%% Correct traces and save/return result
if isa(input,'Traces') || ischar(input)
    [varargout{1:nargout}] = gcorr(input, gammaFactor);

elseif iscell(input)
    for i=1:numel(input)
        gcorr(input{i}, gammaFactor);
    end
end


end %function gammacorrect





%% Actually correct for gamma.
function varargout = gcorr(input, gammaFactor)

% Load fluorescence data
if ischar(input)
    data = loadTraces( input );
elseif isa(input,'TracesFret')
    data = input;
else
    error('Invalid input. Must be TracesFret object or path to .traces file');
end

if isChannel(data,'acceptor2')
    disp('Warning: this function may not work correctly with multi-color FRET');
end

% To do: Estimate the mean gamma value across all traces in the file

% Keep track of adjustments in trace metadata.
% Change Donor scale factor 
idxA1 = strcmpi(data.channelNames,'donor');
[data.traceMetadata.scaleFluor] = to_col(data.traceMetadata.scaleFluor);
% Get scaleFluor matrix with the old value
scaleFluor = cat(2,data.traceMetadata.scaleFluor);

% Update the scaleFluor matrix with the new gamma Factor. The old
% gamma factor is still saved in metadata.
scaleFluor(idxA1,:) = gammaFactor; 

% correct traces with new gamma
data = correctDonCell(data,scaleFluor,[]);
% update the Fret values
data.recalculateFret();

% Save resulting data
if ischar(input) && nargout==0
    [p,f,e] = fileparts( input );
    %reset filename if it has already extension '_g'
    if ~contains([f e],'_g.')
        outFilename = fullfile(p, [f '_g' e]);
    else
        outFilename = fullfile(p, [f e]);
    end
    saveTraces( outFilename, data );
end

output = {data,gammaFactor};
[varargout{1:nargout}] = output{1:nargout};

end




