function [plotWindow] = rtdPlots(varargin)
% RTDPLOTS() Real-time delivery (RTD) plotting
%
% RTDPLOTS()
% Prompts user to select trace files and recreates plots in  the format 
% generated by rtdTool().
%
% RTDPLOTS(options)
% The struct 'options' accepts the same user-defined settings as rtdTool(),
% but only handles display-related settings.
%
% plotWindow = RTDPLOTS()
% plotWindow = RTDPLOTS(options)
% Returns a handle, plotWindow, to the generated plot window.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


plotWindow = 0;

% Handle any user-specified settings.
if nargin
    opt = varargin{1};
end

% Set plot-only mode for rtdTool().
opt.plotOnly = 1; 

% Prompt user for trace files if not specified.
if ~isfield(opt, 'fileList')
    filter = {'*.traces','Binary Traces Files (*.traces)'; ...
              '*.rawtraces','Raw Traces Files (*.rawtraces)'; ...
              '*.*','All Files (*.*)'};
    prompt = 'Select trace files. Hit cancel when finished.';
    opt.fileList = getFiles(filter, prompt);
    
    if numel(opt.fileList)==0, return; end %no files selected
end

% Run rtdTool().
plotWindow = rtdTool(opt);

end
