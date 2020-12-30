function [time,velocity]=particleVelocity(dt,xCor,yCor)

%--------------------------------------------------------------------------
%
% particleVelocity:
%   Calculates the speed from position and time
% 
% Syntax:  
%   [time,velocity]=particleVelocity(dt,xCor,yCor)
% 
% Inputs:
%   1. dt - frame rate in seconds 
%   2. xCor -  x coordinate vector
%   3. yCor - y coordinate vector
% 
% Outputs:
%   1. time - time vector
%   2. velocity - 1D velocity field vector
% 
% See also: 
%   smCellViewtraces
%
% Authors: 
%   - P.G. Jun 2011
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

% Determine dimesion of the trace
ind=xCor(1,:)>0;
dimX=size(ind(ind>0),2);

% Allocate Memory
len=dimX;
velMat=zeros(5,len-1);

% Calculate velocity of a single trace
if len>=2
    for h=2:len
        velMat(1,h-1) = (h-1)*dt;                                %vector dt
        velMat(2,h-1) = xCor(1,h)-xCor(1,h-1);                   %vector dx
        velMat(3,h-1) = yCor(1,h)-yCor(1,h-1);                   %vector dy
        velMat(4,h-1) = (velMat(2,h-1)^2 + velMat(3,h-1)^2)^0.5; %vector dr
        velMat(5,h-1) = velMat(4,h-1)/dt;                        %vector dv
    end
    %dr       = [velMat(1,:)' velMat(2,:)' velMat(3,:)' velMat(4,:)' velMat(5,:)' ] 
    time     = velMat(1,:);
    velocity = velMat(5,:);
else
    time     = dt;
    velocity = 0;
end    


end

