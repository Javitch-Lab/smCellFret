function data = correctDonCell(data, scaleFluor, indexes)
% correctTraces  Apply crosstalk and scaling corrections to fluorescence traces.
%
%    DATA = correctTraces(DATA, CROSSTALK, SCALING) modifies the Traces object DATA
%    in place (left hand argument is optional), applying the given CROSSTALK
%    and channel SCALING corrections in that order. Use data.recalculateFret
%    to update FRET values from corrected fluorescence traces.
%      CROSSTALK(source channel number, destination channel number, trace ID)
%      SCALING(channel number, trace ID)
%
%    ... = correctTraces(..., INDEXES) only adjustes the traces listed in the
%    vector INDEXES. CROSSTALK and SCALING include values for all traces.
%
%    If the final dimension of either CROSSTALK or SCALING is unity, the
%    value will be applied to all traces.
%
%    This function does not alter FRET traces. Use the TracesFret(4) method
%    recalculateFret to update FRET traces.
%
%    See also crosstalkcorrect, gammacorrect, scaleAcceptor, bgsub.
%
%   Copyright 2007-2015 Cornell University All Rights Reserved.
%
% =========================================================================
%
%  correctDonCell is a modified version of the SPARTAN function
%  correctTraces 
%
% =========================================================================
%
% correctDonCell:
%   correctDonCell is a sub-function of gammaCorrectCell. 
%
% Description:
%   The function applies the gamma correction to the actual data.
%
% Syntax:  
%   data = correctDonCell(data, scaleFluor, indexes)
% 
% Inputs:
%   data       - traces object
%   scaleFluor - gamma factor matrix
%   indexes    - only adjust the traces listed in the vector indexes
% 
% Outputs:
%   data      - traces object with corrected donor intensities
%
% See also: 
%   gammaCorrectCell.m
%
% Author: 
%   P.G. Sep 2020
%
% Copyright:
%   2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

narginchk(2,5);
nargoutchk(1,1);

nFluor = numel(data.idxFluor);


% If no traceMetadata (called from gettraces), use defaults.
if ~isfield(data.traceMetadata,'scaleFluor')
    [data.traceMetadata.scaleFluor] = deal( ones(nFluor,1) );
end

% If no updated scaling values given, use current value (no change).
if isempty(scaleFluor)
    % Convert vectors to the same, correct, orientation.
    [data.traceMetadata.scaleFluor] = to_col(data.traceMetadata.scaleFluor);
    scaleFluor = cat(2,data.traceMetadata.scaleFluor);
end

% % If a single value is supplied, apply it to all traces.
% if size(scaleFluor,2)==1 % only one element
%     scaleFluor = repmat(scaleFluor, [1 data.nTraces]);
% end

% Verify input argument sizes match
if isempty(indexes)
    indexes = 1:data.nTraces;
end

if islogical(indexes), indexes=find(indexes); end

if ~(max(indexes)<=data.nTraces && max(indexes)<=size(scaleFluor,2))
    error('Input argument size mismatch');
end


%% Undo previous corrections, reverting to state as acquired.
chNames = data.channelNames(data.idxFluor);  %increasing wavelength order.

% Undo previous scaling by dividing fluorescence by the factor in metadata.
for ch=1:nFluor
    chscale = arrayfun( @(x)x.scaleFluor(ch), data.traceMetadata(indexes) )';
    data.(chNames{ch})(indexes,:) = bsxfun(@rdivide, data.(chNames{ch})(indexes,:), chscale);  %indexes
end

%% Apply new corrections
% Scale fluorescence to correct for unequal brightness/sensitivity
for ch=1:nFluor
    data.(chNames{ch})(indexes,:) = bsxfun(@times, ...
                      data.(chNames{ch})(indexes,:), scaleFluor(ch,indexes)');
end

% Save new correction parameters in metadata
for i=to_row(indexes)
    data.traceMetadata(i).scaleFluor = scaleFluor(:,i);
end



end %function correctTraces


