%--------------------------------------------------------------------------
%
% scriptDiff2Frettraces.m:
%   Assignment of DC-MSS diffusion data to FRET data.
% 
% Description:
%   The script links track segments classified as either free, confined, 
%   directional or immobile diffusion to the FRET data. Updated FRET 
%   data, which contain motion dynamics can then be displayed in the
%   cellFretViewtraces program for further analysis.      
% 
% Syntax:  
%   scriptDiff2Frettraces.m
% 
% Inputs: 
%   The script prompts the user to provide Ch1 AND Ch2 folder locations.   
%   1. Ch1 folders must contain files of type 'SegResultsFinal.mat' (output 
%   of scriptGetSegResults.m)
%   2. Ch2 folders must contain files of type 'fretTracesPst.mat' (output of
%   scriptFretFilter.m)
%   3. OPTIONAL: If the Channel Ch2 folder also contains files of type 
%   'SegResultsFinal.mat', scriptDiff2Frettraces will also process this data.
% 
% Outputs:
%   1. fretTracesDif.mat - Structure variable with the additional field 
%   Diff, which contains information about the duration of the four
%   different diffusion modes. Diff has the following subfields: 
%   a) segResFinal: contains the results of the DC-MSS motion analysis    
%   b) idlIm: binary n x m matrix with n traces and m frames. Thus each row 
%      represents a time trace, which takes the value 1 if the particle is 
%      immobile, otherwise 0.  
%   c) idlCnf: binary matrix which takes the value 1 if the particle is 
%      confined, otherwise 0.   
%   d) idlFre: binary matrix which takes the value 1 if the particle is 
%      freely diffusing, otherwise 0.  
%   e) idlSpr: binary matrix which takes the value 1 if the particle is 
%      in the mode of directed diffusing, otherwise 0.  
%   Note: fretTracesDif.mat can be loaded into cellFretViewtraces for
%   manual inspection.
% 
% See also: 
%   scriptFretFilter.m, scriptGetSegResults.m, cellFretViewtraces.   
%
%
% Author: 
%   - P.G. Mar 2019
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Get Directories (Ch1 and Ch2 directories required)

clear
dirArray = getDirs('Highlight and Open ACCEPTOR (Ch1) AND DONOR (Ch2) Folders (HIT CANCEL TO END SELECTION)');
nExperiments = length(dirArray(1,:));
if nExperiments==0,return,end

%% Loop through all experiments

% Create waitbar
h = waitbar(0,'Adding Diffusion Data ....');

% Set current Directopry
for j=1:nExperiments 
    
    waitbar(j / (nExperiments))
    
    % Set current Directopry
    basepath = dirArray{2,j};
    dataSet  = dirArray{1,j};
    currentDir = [basepath filesep dataSet filesep];
    
    disp(dataSet);
    
    % Load SegResultsFinal (in Dir Ch1) 
    if contains(dataSet,'Ch1') && ...
            exist([currentDir 'SegResultsFinal.mat'],'file')
        diffCh1=load([currentDir 'SegResultsFinal.mat']);

    % Load post filtered fretTraces & SegResultsFinal (in Dir Ch2) 
    elseif contains(dataSet,'Ch2')
        
        load([currentDir 'fretTracesPst.mat']);
   
        %------------------------------------------------------------------
        %--- Add Diffusion as field to fretTracesPst 
        %------------------------------------------------------------------
        disp('add diffusion data...');
        
        if exist([currentDir 'SegResultsFinal.mat'],'file')
            diffCh2=load([currentDir 'SegResultsFinal.mat']);
            try
                fretTraces.Diff.Ch2.segResFinal=diffCh2.segResultsFinal;
            catch 
                %old format
                fretTraces.Diff.Ch2.segResFinal=diffCh2.s;
            end
        else
            warning('no SegResultsFinal.mat in Ch2') 
        end
        
        try
            fretTraces.Diff.Ch1.segResFinal=diffCh1.segResultsFinal;
        catch
            % old format
            fretTraces.Diff.Ch1.segResFinal=diffCh1.s;
        end
        
        
        %------------------------------------------------------------------
        %--- Add logical matrix for diffuson segment analysis  
        %------------------------------------------------------------------
        disp('add logical idealization data for segment analysis...');
        
        nTraces=size(fretTraces.Ch1.x,1);
        
        % load idl FRET Data (donorRange)
        idl     = fretTraces.Fret.idlTotal;
        
        % Initialize new logical arrays
        idlImob = false(size(fretTraces.Fret.idlTotal));
        idlConf = false(size(fretTraces.Fret.idlTotal));
        idlFree = false(size(fretTraces.Fret.idlTotal));
        idlSupr = false(size(fretTraces.Fret.idlTotal));
        
        for k=1:nTraces   
            % Get number of trace segments
            ids = fretTraces.Diff.Ch1.segResFinal(:,23); % the ids vector (column)
            % Get number of segments in a track
            lia = ismember(ids,k); %lia is a logical vector (column)
            % Find the idx where lia is logical true 
            idx = lia == 1;
            % Get Segments
            seg    = fretTraces.Diff.Ch1.segResFinal(idx,1:2);
            state  = fretTraces.Diff.Ch1.segResFinal(idx,3);
            % Plot Fret traces
            nSegs = sum(lia);
            for l = 1:nSegs 
                range=seg(l,1):seg(l,2);
                switch state(l) % Chk state of segment l in trace k
                    case 0 % Immobile
                        idlImob(k,range)=idl(k,range);
                    case 1 % Confined
                        idlConf(k,range)=idl(k,range);
                    case 2 % Free
                        idlFree(k,range)=idl(k,range);
                    case 3 % Super
                        idlSupr(k,range)=idl(k,range);
                end % Switch Cmd
            end % Segments
        end % Traces
        
        % Save data
        fretTraces.Diff.Ch1.idlImb = idlImob;
        fretTraces.Diff.Ch1.idlCnf = idlConf;
        fretTraces.Diff.Ch1.idlFre = idlFree;
        fretTraces.Diff.Ch1.idlSpr = idlSupr;
        
        %------------------------------------------------------------------
        %--- Save new fretTraces matrix
        %------------------------------------------------------------------
        outfile=strrep('fretTracesPst.mat','Pst','Dif');
        save([currentDir outfile],'fretTraces');
        
    else 
        warning('file does not exist')
        continue
       
    end

    
end % End Experiment Loop
close(h);
disp('Finished'); 

%% End