function varargout = deltaCorrectCell(input, deltaFactor)
%--------------------------------------------------------------------------
% 
% deltaCorrectCell:
%   Performs a direct excitation correction.
% 
% Description:
%   The function deltaCorrectCell.m multiplies the mean total trace 
%   intensities by a correction factor delta and subtracts these adjusted 
%   intensities from the acceptor intensities. This correction is based on 
%   the assumption that the direct acceptor excitation intensity can be
%   approximated by the total intensity. Here the delta factor is a measure 
%   of the percentage of direct acceptor excitation with the donor laser.
% 
% Syntax:  
%   varargout = deltaCorrectCell(input, deltaFactor)
% 
% Inputs:
%   input       - traces file with the extension *_a.traces
%   deltaFactor - empirically determined delta factor
%   
%   When executing the function, without input arguments the user is prompted to 
%   1. select a list of traces files that are already crosstalk corrected 
%      (files with extension *_a.traces)
%   2. confirm the default delta factor or enter a new value
% 
% Outputs:
%   Corrected traces are saved in the same directory as the input file with
%   the additional extension d for delta (e.g. *sync_flt_a.traces)
% 
% Other m-files required: 
%   Subfunctions: correctDExtCell
% 
% See also: 
%   alphaCorrectCell.m,  gammaCorrectCell.m
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
nargoutchk(0,4);

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

%% Set factor for direct excitation by user
if nargin<2 || isempty(deltaFactor)
    
    deltaDefault = '0.056';
    delta  = inputdlg({'Enter new Value:'},...           %prompt
                     'Set Global Direct Exitation Factor',... %dlg titel
                     [1 55],...                       %width of line
                     {deltaDefault});                 %default answer

    if isempty(delta) || isnan(str2double(delta{1,1}))
        % user hit cancel or input field was empty
        return
    else
        deltaFactor = str2double(delta{1,1});
    end
    
end

%% Correct traces and save/return result
if isa(input,'Traces') || ischar(input)
    [varargout{1:nargout}] = dcorr(input, deltaFactor);

elseif iscell(input)
    % no output when using file explorer to load traces file(s). 
    for i=1:numel(input)
        dcorr(input{i}, deltaFactor);
    end %for each file
end


end %function delta correct



%% Actually correct for direct acceptor excitation
function varargout = dcorr(input, deltaFactor)

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

% Add field for direct Excitation if it does not exist
% directExtn field with entries for
% position 1: deltaFactor = 0 
% position 2: calc mean total intensities per trace
% position 3: actual correction factor -> deltaFactor*meanTotalInt

if ~isfield(data.traceMetadata,'directExtn')
    [data.traceMetadata.directExtn] = deal( zeros(1, 3) );
    for i=1:data.nTraces
        
        % Initialize the idealization matrix with logical zeros
        idl = false(1,data.nFrames);
        
        % set idl-matrix TRUE over the donor range
        donSOT= data.traceMetadata(i).startOfTraceDon; %unit is frames
        donEOT= data.traceMetadata(i).endOfTraceDon; %unit is frames
        idl(donSOT:donEOT) = true;
        
        % Calc mean total intensity per trace over donor range 
        % Subsitute zeros with NaN to get correct mean total intensities 
        total_idl2D = data.total(i,:).*idl;
        total_idl2D(total_idl2D==false)=NaN;
        data.traceMetadata(i).directExtn(1,2)= nanmean(total_idl2D,2);
        
        
        
    end        
end

% Create a new variable directExtn with the new delta factor and total
% intensities (from prev calculation) and pass this on for new directExtn 
% correction
indexes = 1:data.nTraces;
directExtn(1,indexes) = deltaFactor; % pos1 -> delta factor
directExtn(2,indexes) = arrayfun( @(x) x.directExtn(2), data.traceMetadata(indexes) )'; % pos2 -> mean total per trace 
directExtn(3,indexes) = bsxfun(@times, directExtn(1,indexes), directExtn(2,indexes))'; % pos3 ->  pos1 * pos2

% Make the correction and recalculate FRET
data = correctDExtCell(data, directExtn,[]);

% Recalculate total intensities and FRET based on corrected datacorrectD
thresh=repmat(500,data.nTraces,1);
data.recalculateFret(thresh);

% Check distributions 
meanTotalPerTrace = arrayfun( @(x)x.directExtn(2), data.traceMetadata(indexes) )';
extnCorrAccInt = arrayfun( @(x)x.directExtn(3), data.traceMetadata(indexes) )';

% Save resulting data with extension d (for deltaFactor correction)
if ischar(input) && nargout==0
    [p,f,e] = fileparts( input );
    if ~contains([f e],'_d.')
        outFilename = fullfile(p, [f '_d' e]);
    else
        outFilename = fullfile(p, [f e]);
    end
    saveTraces( outFilename, data );
end

output = {data,deltaFactor,meanTotalPerTrace,extnCorrAccInt};
[varargout{1:nargout}] = output{1:nargout};

end %FUNCTION DCORR


