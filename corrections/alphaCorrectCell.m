function varargout = alphaCorrectCell(input, crtFactor)
% -------------------------------------------------------------------------
%
% alphaCorrectCell:
%  Performs a spectral crosstalk correction.
% 
% Description:
%   The function alphaCorrectCell multiplies the donor intensities by a 
%   correction factor alpha and subtracts these adjusted intensities from 
%   the acceptor intensity. The alpha factor for spectral crosstalk takes 
%   into account the amount of donor fluorescence leakage into the acceptor 
%   emission channel. 
%
% Syntax:  
%   varargout = alphaCorrectCell(input, crtFactor)
% 
% Inputs:
%   1. input     -  traces file with the extension *.traces
%   2. crtFactor -  empirically determined crosstalk value. 
%   
%   When executing the function, without input arguments the user is prompted to 
%   1. select a list of trace files (e.g. *sync_flt.traces) to be corrected.
%   2. confirm the default alpha factor or enter a new value
% 
% Outputs:
%   Corrected traces are saved in the same directory as the input file
%   withg
%   the additional extension 'a' for alpha (e.g. *sync_flt_a.traces).   
% 
% Other m-files required: 
%   Subfunctions: correctCrtCell.m
% 
% See also: 
%   deltaCorrectCell.m,  gammaCorrectCell.m
%
% Permissions: 
%   alphaCorrectCell is a revised version of SPARTAN's [1] crosstalkorrect.m'
%   function. The function 'crosstalkorrect.m' has been modified with  
%   permission from the Blanchard Laboratory. 
%
% References:
%   [1] Juette, M. F., et al. (2016). "Single-molecule imaging of 
%   non-equilibrium molecular ensembles on the millisecond timescale." 
%   Nature Methods 13(4): 341-344.
%
% Authors: 
%   -P.G. Sep 2020
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
    input = getFiles('*.traces;*.rawtraces');
    if isempty(input), return; end  %user hit cancel
    
elseif iscell(input)
    if ~all(cellfun(@ischar,input))
        error('Input must be a Traces object, filename string, cell array of filenames');
    end
    if nargout>0
        error('Output only supported for scalar input')
    end
end

%% Set crosstalk value by user
if nargin<2 || isempty(crtFactor)
    
    crtDefault = '0.09417';
    crt  = inputdlg({'Enter new Value:'},...        %prompt
                     'Set Global Crosstalk Factor',... %dlg titel
                     [1 50],...                     %width of line
                     {crtDefault});                 %default answer

    if isempty(crt) || isnan(str2double(crt{1,1}))
        % user hit cancel or input field was empty
        return
    else
        crtFactor = str2double(crt{1,1});
    end
    
end

%% Correct traces and save/return result
if isa(input,'Traces') || ischar(input)
   [varargout{1:nargout}] = crtCorr(input, crtFactor);

elseif iscell(input)
    % no output when using file explorer to load traces file(s). 
    for i=1:numel(input)
        crtCorr(input{i}, crtFactor);
    end %for each file
end


end %function crosstalkcorrect



%% Actually correct the crosstalk
function varargout = crtCorr(input, ctval)

if ischar(input)
    data = loadTraces(input);
elseif isa(input,'Traces')
    data = input;
else
    error('Invalid input argument. Must be a filename or Traces object');
end

if isChannel(data,'acceptor2')
    disp('Warning: this function may not work correctly with multi-color FRET');
end

% To do: Estimate the crosstalk if not given by the user.

% Make the crosstalk correction and recalculate FRET
data = correctCrtCell(data, ctval,[]); %crosstalk is new value
data.recalculateFret();

% Save resulting data
if ischar(input) && nargout==0
    [p,f,e] = fileparts( input );
    %reset filename if it has already extension '_a'
    if ~contains([f e],'_a.')
        outFilename = fullfile(p, [f '_a' e]);
    else
        outFilename = fullfile(p, [f e]);
    end
    saveTraces( outFilename, data );
end

output = {data,ctval};
[varargout{1:nargout}] = output{1:nargout};

end %FUNCTION CrtCorr


