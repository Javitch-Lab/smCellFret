
function [ resultsCumFit,tData,yData,data ] = cumFit_lt( ltValues,dt,nExp )
%--------------------------------------------------------------------------
%
% cumFit_lt.m:
%   The function performs a cumulative fitting of up to 4 exponential
%   functions to a data set. 
% 
% Syntax:  
%   [ resultsCumFit,tData,yData,data ] = cumFit_lt( ltValues,dt,nExp )
% 
% Inputs:
%   1. ltValues - Column vector with n elements, which contains the
%      measured track length values (unit is frames)
%   2. dt - Variable that keeps the time resolution (unit is s)
%   3. nExp - Number of exponential functions to be fitted to the data
% 
% Outputs:
%   1. cumFitResults - The workspace variable contains the lifetimes, 
%      amplitudes, BIC and AIC states and a Goodnes of Fit Statistics. 
%   2. tData - Time axis data (raw data)
%   3. yData - Cumulative lifetime data (raw data)
%   4. data - This structure variable has the following four fields:  
%      a) data.t = same as 2
%      b) data.y = same as 3
%      c) data.tFit = Time axis data of the fitting function.
%      d) data.yFit = Cumulative lifetime data of the fitting function.
%   5. Fitting reults are also shown in the matlab command window (see 
%      example below).    
%
% Example: 
%   Output in the matlab command window: 
%     Fitting Parameters
%          General model:
%          fcdf(x) = a1*(1-exp((-x+x0)/t1))+(1-a1)*(1-exp((-x+x0)/t2))
%          Coefficients (with 95% confidence bounds):
%            a1 =      0.7072  (0.6903, 0.7241)
%            x0 =      0.7456  (0.7348, 0.7564)
%            t1 =       1.369  (1.324, 1.415)
%            t2 =       8.786  (8.234, 9.337)
%
%     Fitting Errors (Goodness of Fit Statistics)
%                sse: 0.0083
%            rsquare: 0.9991
%                dfe: 144
%         adjrsquare: 0.9991
%               rmse: 0.0076
% 
%     Information Criterion AIC / BIC:
%         aic: 1.6391e+03
%         bic: 1.6504e+03
% 
% See also: 
%   scriptCumFit.m
%
% Authors: 
%   - P.G.
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


%% Load data for testing
% clear;
% load('combinedData.mat');
% ltValues = [cellData.traceMetadata.traceLen]';
% ltValues = ltValues - 1;
% dt       = 0.04;
% nExp     = 2;

%% Load Data

nSeg           = size(ltValues,1);
lt             = sort(ltValues*dt);

% Compute the empirical cumulative distribution function from 
[y,t]          = ecdf(lt);

% Prepare data inputs for curve fitting using 'fit' function
[tData, yData] = prepareCurveData( t, y );

%% Initialize Fitting Models

% Choose model for lifetime fit
expModel=nExp;

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

% Display option in the command window
opts.Display = 'final';

% Robust linear least-squares fitting method
%opts.Robust = 'Off';
opts.Robust = 'LAR'; % specifies the least absolute residual method.
%opts.Robust = 'Bisquare'; % specifies the bisquare weights method.

% Algorithm to use for the fitting procedure
% opts.Algorithm ='Trust-Region'; %to be used with bound constraints 
opts.Algorithm ='Levenberg-Marquardt';

opts.DiffMinChange = 1.0e-8;
opts.DiffMaxChange = 0.1;
opts.MaxFunEvals   = 600;
opts.MaxIter       = 400;
opts.TolFun        = 1.0e-6;
opts.TolX          = 1.0e-6;
opts.Exclude       = [];%1:15; % indices of t values that are exclude 


switch expModel
    
    % 1 exp fit
    case 1
        ft   = fittype( '1*(1-exp((-x+x0)/t1))',...
                        'independent', 'x',...
                        'dependent', 'y',...
                        'coefficients',{'x0','t1'});
        % Set lower and upper bounds  
          opts.Lower      = [0.8 0.000001]; % [x0 t1]
          opts.Upper      = [Inf Inf];    % [x0 t1]
         opts.StartPoint = [tData(1) 1];   % [x0 t1]
    
    % 2 exp fit
    case 2
        ft   = fittype( 'a1*(1-exp((-x+x0)/t1))+(1-a1)*(1-exp((-x+x0)/t2))',...
                        'independent', 'x',...
                        'dependent', 'y',...
                        'coefficients',{'a1','x0','t1','t2'});
        %Set lower and upper bounds (only for Trust Region)  
        opts.Lower      = [0 0.08 0.000001 1]; % [a1 x0 t1 t2]
        opts.Upper      = [Inf Inf 5 50];       % [a1 x0 t1 t2]
         opts.StartPoint = [0.5 tData(1) 0.1 5];       % [a1 x0 t1 t2]
        
    % 3 exp fit    
    case 3
        ft   = fittype( 'a1*(1-exp((-x+x0)/t1))+a2*(1-exp((-x+x0)/t2))+(1-a1-a2)*(1-exp((-x+x0)/t3)) ',...
                        'independent', 'x',...
                        'dependent', 'y',...
                        'coefficients',{'a1','a2','x0','t1','t2','t3'});
        % Set lower and upper bounds  
%         opts.Lower      = [0 0 0.000001 0.000001 13]; % [a1 a2 t1 t2 t3]
%         opts.Upper      = [Inf Inf Inf Inf 17];            % [a1 a2 t1 t2 t3]
         opts.StartPoint = [0.33 0.33 tData(1) 0.1 1 10];     % [a1 a2 t1 t2 t3]
       
    % 4 exp fit    
    case 4
        ft   = fittype( 'a1*(1-exp((-x+x0)/t1))+a2*(1-exp((-x+x0)/t2))+(1-a1-a2)*(1-exp((-x+x0)/t3))+(1-a1-a2-a3)*(1-exp((-x+x0)/t4)) ',...
                        'independent', 'x',...
                        'dependent', 'y',...
                        'coefficients',{'a1','a2','a3','x0','t1','t2','t3','t4'});
        % Set lower and upper bounds  
%         opts.Lower      = [0 0 0 tData(1) 0.000001 0.000001 0.000001 0.000001]; % [a1 a2 a3 x0 t1 t2 t3 t4]
%         opts.Upper      = [Inf Inf Inf tData(1) Inf Inf Inf Inf];               % [a1 a2 a3 x0 t1 t2 t3 t4]
         opts.StartPoint = [0.25 0.25 0.25 tData(1) 0.01 0.1 1 10];                % [a1 a2 a3 x0 t1 t2 t3 t4]
end

%% Cumulative Histogram

fig1=figure('Position',[100 100 500 850]); 

% Fit Distribution
clear results
a1=[];a2=[];a3=[];a4=[];
t1=[];t2=[];t3=[];t4=[];

switch expModel
    case 1
        % Fit 1 exp model to data.
         try   
            % Plot Graph 
            ax1=subplot(3,1,1);
            
            [fcdf, gof] = fit( tData, yData, ft, opts );
            f1=stairs(t, y,'Marker','.','LineStyle','-','Color','b');hold on
            data.y=y;
            data.t=t;
            
            tRegion  = tData(~ismember(1:numel(tData),opts.Exclude));
            yRegion  = yData(~ismember(1:numel(tData),opts.Exclude));
            t=[0;tRegion];
            f2=line([0;t], fcdf([0;t]),'Color','red');hold on
            disp([' nSeg = ' num2str(nSeg)])
            legend([f1 f2 ],...
                   [' nSeg = ' num2str(nSeg)],...
                   [num2str(expModel) ' Exp Fit'],...
                   'Location', 'SouthEast');
            ylabel('CDF');ylim([0,1.1]);
            grid on;grid minor 
            
            data.tFit=t;
            data.yFit=fcdf(t);
            
            % Plot residuals.
            ax2=subplot(3,1,2);
            g = plot( fcdf, tRegion, yRegion, 'residuals' );
            legend( g, 'Residuals', 'Zero Line', 'Location', 'SouthEast' );
            ylabel Residuals
            grid on;grid minor
            
            % Define Probability Density Function
            ax3=subplot(3,1,3);
            a1=1;
            t1=fcdf.t1; 
            %t1 = 1/1.26562; %Bic=23497;
            cumPdf = @(x) a1*(1/t1)*exp(-x/t1);
                  
        catch
            title('fit does not converge');
            warning('Choose Different Start Point');
            resultsCumFit=[];
            return
        end
        
    case 2
        % Fit 2 exp model to data.
        try
            % Plot Graph
            ax1=subplot(3,1,1);
            [fcdf, gof] = fit( tData, yData, ft, opts );
            f1=stairs(t, y,'Marker','.','LineStyle','-','Color','b');
            data.y=y;
            data.t=t;

            tRegion  = tData(~ismember(1:numel(tData),opts.Exclude));
            yRegion  = yData(~ismember(1:numel(yData),opts.Exclude));
            t=[0;tRegion];
            f2=line(t, fcdf(t),'Color','red');hold on
            legend([f1 f2 ],...
                   [' nSeg = ' num2str(nSeg)],...
                   [num2str(expModel) ' Exp Fit'],...
                   'Location', 'SouthEast');
            ylabel('CDF');ylim([0,1.1]);
            grid on;grid minor 
            
            data.tFit=[0;t];
            data.yFit=fcdf([0;t]);
            
            % Plot residuals.
            ax2=subplot(3,1,2);
            g = plot( fcdf, tRegion, yRegion, 'residuals' );
            legend( g, 'Residuals', 'Zero Line', 'Location', 'SouthEast' );
            ylabel Residuals
            grid on; grid minor
            
            % Define Probability Density Function
            ax3=subplot(3,1,3);
            a1=fcdf.a1;a2=1-a1;
            t1=fcdf.t1;t2=fcdf.t2;
            cumPdf = @(x) a1*(1/t1)*exp(-x/t1) + a2*(1/t2)*exp(-x/t2);
           
        catch
            title('fit does not converge');
            warning('Choose Different Start Point');
            resultsCumFit=[];
            return
        end

    case 3
        % Fit 3 exp model to data.
        try
            % Plot Graph
            ax1=subplot(3,1,1);
            [fcdf, gof] = fit( tData, yData, ft, opts );
            f1=stairs(t, y,'Marker','.','LineStyle','-','Color','b');
            data.y=y;
            data.t=t;
            
            tRegion  = tData(~ismember(1:numel(tData),opts.Exclude));
            yRegion  = yData(~ismember(1:numel(yData),opts.Exclude));
            t=[0;tRegion];
            f2=line(t, fcdf(t),'Color','red');hold on
            legend([f1 f2 ],...
                   [' nSeg = ' num2str(nSeg)],...
                   [num2str(expModel) ' Exp Fit'],...
                   'Location', 'SouthEast');
            ylabel('CDF');ylim([0,1.1]);
            grid on;grid minor 
            
            data.tFit=[0;t];
            data.yFit=fcdf([0;t]);
            
            % Plot Residuals
            ax2=subplot(3,1,2);
            g = plot( fcdf, tRegion, yRegion, 'residuals' );
            legend( g, 'Residuals', 'Zero Line', 'Location', 'SouthEast' );
            ylabel Residuals
            grid on; grid minor
            
            % Define Probability Density Function
            ax3=subplot(3,1,3);
            a1=fcdf.a1;a2=fcdf.a2;a3=(1-a1-a2);
            t1=fcdf.t1;t2=fcdf.t2;t3=fcdf.t3;   
            cumPdf = @(x) a1*(1/t1)*exp(-x/t1) + a2*(1/t2)*exp(-x/t2) + a3*(1/t3)*exp(-x/t3);
   
        catch
          title('fit does not converge');
          warning('Choose Different Start Point');
          return
        end
        
        case 4
        % Fit 4 exp model to data.
        try 
            % Plot Graph
            ax1=subplot(3,1,1);
            [fcdf, gof] = fit( tData, yData, ft, opts );
            f1=stairs(t, y,'Marker','.','LineStyle','-','Color','b');
            data.y=y;
            data.t=t;
            
            tRegion  = tData(~ismember(1:numel(tData),opts.Exclude));
            yRegion  = yData(~ismember(1:numel(yData),opts.Exclude));
            t=[0;tRegion];
            f2=line(t, fcdf(t),'Color','red');hold on
            legend([f1 f2 ],...
                   [' nSeg = ' num2str(nSeg)],...
                   [num2str(expModel) ' Exp Fit'],...
                   'Location', 'SouthEast');
            ylabel('CDF');ylim([0,1.1]);
            grid on;grid minor 
            
            data.tFit=[0;t];
            data.yFit=fcdf([0;t]);
            
            % Plot Residuals
            ax2=subplot(3,1,2);
            g = plot( fcdf, tRegion, yRegion, 'residuals' );
            legend( g, 'Residuals', 'Zero Line', 'Location', 'SouthEast' );
            ylabel Residuals
            grid on; grid minor
            
            % Define Probability Density Function
            ax3=subplot(3,1,3);
            a1=fcdf.a1;a2=fcdf.a2;a3=fcdf.a3;a4=1-a1-a2-a3;
            t1=fcdf.t1;t2=fcdf.t2;t3=fcdf.t3;t4=fcdf.t4;   
            cumPdf = @(x) a1*(1/t1)*exp(-x/t1) + a2*(1/t2)*exp(-x/t2) + a3*(1/t3)*exp(-x/t3) + a4*(1/t4)*exp(-x/t4);
   
        catch
          title('fit does not converge');
          warning('Choose Different Start Point');
          resultsCumFit=[];
          return
        end
end

%% Calculate AIC / BIC
logLvector = [lt fcdf(lt) log(cumPdf(lt))];
logL       = sum(logLvector(:,3));
numParam   = 2*expModel-1;
numObs     = nSeg;
[aic,bic]= aicbic(logL,numParam,numObs);
criterion.aic = aic;
criterion.bic = bic;

%--- Display and Plot Fitting Results
disp('Fitting Parameters');
disp(fcdf);
disp('Fitting Errors');
disp(gof);
disp('Information Criterion');
disp(criterion);
h1=line(t, cumPdf(t),'Color','red');
h2=line(t, cumPdf(t),'Color','red');hold on
legend([h1 h2],...
        ['AIC: ' num2str(aic)],...
        ['BIC: ' num2str(bic)],...
        'Location', 'SouthEast' );
xlabel('Time (s)'); ylabel PDF
grid on; grid minor
linkaxes([ax1,ax2,ax3],'x');

% Save Fitting Results

ltVar={'Lt1','Lt2','Lt3','Lt4'};
ampVar={'A1','A2','A3','A4'};

% Initialize result structure
for i=1:numel(ltVar)
    fn=ampVar{i};
    resultsCumFit(1).(fn) = []; 
    fn=ltVar{i};
    resultsCumFit(1).(fn) = []; 
end

% Sort and save results
resultLt  = [t1,t2,t3,t4];
resultAmp = [a1,a2,a3,a4];
[~,idx]=sort(resultLt);
for i=1:numel(idx)
    fn=ampVar{i};
    resultsCumFit(1).(fn) = resultAmp(idx(i));   
    fn=ltVar{i};
    resultsCumFit(1).(fn) = resultLt(idx(i));
end

% Save information criteria
resultsCumFit(1).BIC = bic;
resultsCumFit(1).AIC = aic;

% Save fitting errors (goodness of fit statistics)
% Reference: Wiki coefficient of determination
%            Wiki Root-mean-square deviation
% sse:          Sum of squares due to error
% rsquare:      R-squared (see Wiki: coefficient of determination)
% dfe:          Degrees of freedom in the error
% adjrsquare:   Degree-of-freedom adjusted coefficient of determination
% rmse          Root mean squared error (standard error)
errVar={'sse','rsquare','dfe','adjrsquare','rmse'};
for j=1:numel(errVar)
    fn=errVar{j};
    resultsCumFit(1).(fn) = gof.(fn);
end

%% End













