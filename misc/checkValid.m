function [ output_args ] = checkValid( data )
%--------------------------------------------------------------------------
%
% checkValid:
%   A function to return a warning if the input argument is empty
% 
% Syntax:  
%   [ output_args ] = checkValid( data )
% 
% Inputs:
%   data  - a scalar, a vector, a matrix or a multidimensional array. 
% 
% Outputs:
%   output_args - a warning if data is empty
%
% Authors: 
%   - Jozsef Meszaros 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------


if isempty(data)
   warning('empty data')
end

end

