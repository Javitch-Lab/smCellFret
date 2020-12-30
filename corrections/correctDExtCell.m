function data = correctDExtCell(data, directExtn, indexes )
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
%  correctDExtCell is a modified version of the SPARTAN function
%  correctTraces 
%
% =========================================================================
%
% correctDExtCell:
%   correctDExtCell is a sub-function of deltaCorrectCell. 
%
% Description:
%   The function applies the delta correction to the actual data.
% 
% Syntax:  
%   data = correctDExtCell(data, directExtn, indexes)
% 
% Inputs:
%   data       - traces object
%   directExtn - directExtn matrix
%   indexes    - only adjust the traces listed in the vector indexes
% 
% Outputs:
%   data      - traces object with corrected acceptor intensities
%
% See also: 
%   deltaCorrectCell.m
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
if ~isfield(data.traceMetadata,'directExtn')
    [data.traceMetadata.directExtn]  = deal( zeros(nFluor,nFluor) );
end

% If no updated values for direct excitation given, use current value (no change).
if isempty(directExtn)
    directExtn = cat(1,data.traceMetadata.directExtn)';
end

% Verify input argument sizes match
if isempty(indexes)
    indexes = 1:data.nTraces;
end

if islogical(indexes), indexes=find(indexes); end

if ~(max(indexes)<=data.nTraces && max(indexes)<=size(directExtn,2))
    error('Input argument size mismatch');
end

%% Undo previous corrections, reverting to state as acquired.
% Undo previous direct exitation corrections by loading previous correction
% factor from metadata.

ch2 = 'acceptor';
oldAccIntCorr    = arrayfun( @(x)x.directExtn(3), data.traceMetadata(indexes) )';
data.(ch2)(indexes,:) = data.(ch2)(indexes,:) + oldAccIntCorr;

%% Apply new corrections

% Subtract directExtn
ch2 = 'acceptor';
accIntCorr = directExtn(3,indexes)';
% Get data range
for i=indexes 
    sot = data.traceMetadata(i).startOfTraceAcc;
    eot = data.traceMetadata(i).endOfTraceAcc;
    notZero = data.(ch2)(i,sot:eot)>0;
    data.(ch2)(i,sot:eot) = ((data.(ch2)(i,sot:eot) - accIntCorr(i))).*notZero;
end
%notZero = data.acceptor(:,:)>0;
%data.(ch2)(indexes,:) = (data.(ch2)(indexes,:) - accIntCorr).*notZero;

% Save new correction parameters in metadata
for i=to_row(indexes)
    data.traceMetadata(i).directExtn  = directExtn(:,i)';
end

end %function correctTraces


