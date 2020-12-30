function data = correctCrtCell(data, crosstalk, indexes)
% correctTraces  Apply crosstalk to fluorescence traces.
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
%  correctCrtCell is a modified version of the SPARTAN function
%  correctTraces 
%
% =========================================================================
%
% correctCrtCell:
%   correctCrtCell is a sub-function of alphaCorrectCell. 
%
% Description:
%   The function applies the alpha correction to the actual data.
% 
% Syntax:  
%   data = correctCrtCell(data, crosstalk, indexes)
% 
% Inputs:
%   data      - traces object
%   crosstalk - the crosstalk value
%   indexes   - only adjust the traces listed in the vector indexes
% 
% Outputs:
%   data      - traces object with corrected acceptor intensities
%
% See also: 
%   alphaCorrectCell.m
%
% Author: 
%   P.G. Sep 2020
%
% Copyright:
%   2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

narginchk(2,3);
nargoutchk(1,1);

nFluor = numel(data.idxFluor);


% If no traceMetadata (called from gettraces), use defaults.
if ~isfield(data.traceMetadata,'crosstalk')
    [data.traceMetadata.crosstalk]  = deal( zeros(nFluor,nFluor) );
end


% If no updated crosstalk value is given, use current value (no change).
if isempty(crosstalk)
    % get current value from metadata
    crosstalk = cat(3,data.traceMetadata.crosstalk);
end


% If a single value is supplied, apply it to all traces.
if size(crosstalk,3)==1 % only one element
    crosstalk = repmat(crosstalk, [1 1 data.nTraces]);
end

% Verify input argument sizes match
if isempty(indexes)
    indexes = 1:data.nTraces;
end

if islogical(indexes), indexes=find(indexes); end

if ~(max(indexes)<=data.nTraces && max(indexes)<=size(crosstalk,3))
    error('Input argument size mismatch');
end

%% Undo previous corrections, reverting to state as acquired.

ch1 = 'donor';
ch2 = 'acceptor';

% Get old crt value from metadata
ct = arrayfun( @(x)x.crosstalk(1,2), data.traceMetadata(indexes) )';
oldCrtInt =  bsxfun(@times, data.(ch1)(indexes,:), ct);

% Calculate raw data
data.(ch2)(indexes,:) = data.(ch2)(indexes,:) + oldCrtInt;


%% Apply new corrections

% Subtract new crosstalk
ch1 = 'donor';
ch2 = 'acceptor';
ct = squeeze( crosstalk(1,1,indexes) );
newCrtInt = bsxfun(@times, data.(ch1)(indexes,:), ct);
% Subtract over the whole data range
data.(ch2)(indexes,:) = data.(ch2)(indexes,:) - newCrtInt;

% Save new correction parameters in metadata
for i=to_row(indexes)
    data.traceMetadata(i).crosstalk(1,2)  = crosstalk(:,:,i);
end

end %function correctCrtCell


