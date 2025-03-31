% Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-18-02

function [] = input_test(x, m, t, varargin)
% Conducts basic tests on the input values and throws error messages.
%
% Inputs:
%   x: Time series.
%   m: Embedding dimension.
%   t: Time delay.
%   epsilon: (Optional) Threshold parameter. Default no value.
%% input parser
p = inputParser;

addRequired(p, 'x', @isnumeric);
addRequired(p, 'm', @isnumeric);
addRequired(p, 't', @isnumeric);


addOptional(p, 'epsilon',[], @isnumeric);

parse(p, x, m, t, varargin{:});

%% initialise variables
x = p.Results.x;
m = p.Results.m;  
t = p.Results.t;
epsilon = p.Results.epsilon;

%% input tests 
if isscalar(x) 
    error('Input Error: Time series is length 1.');
else
end

if length(x) < m 
    error('Input error: Embedding dimension is larger than time series.');
else
end

if ~isnumeric(t) || ~isnumeric(m) || ~isscalar(t) || ~isscalar(m) || t < 1 || m < 1 || mod(t,1) ~= 0 || mod(m,1) ~= 0
        error('Input error: Both t and m must be integers larger than 1.');
else
end

if length(x) < t
    error('Input error: Time delay is larger than time series.');
else
end

if isempty(epsilon)
elseif epsilon < 0 || epsilon > 1
    error ('Input error: Epsilon should be between 1 and 0')
else
end

if length(x)<(m-1)*t
    error ('Input error: m and t are to large (m-1)*t must be smaller than the time series length')
else
end

if min(x) <= 0
    error ('Input error: Time series has 0 values')
end

end