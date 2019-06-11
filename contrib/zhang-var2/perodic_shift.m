function y =  perodic_shift(x, shift)
% function y = perodic_shift(x, shift)
%   y = perodic_shift(x, shift) interpolates to find 
%   the signal values at shifted locations
%   
%   y = perodic_shift(x, shift) assumes x is a perodic 
%   signal. Use linear interpolation
% 
%   IN: 
%       x       length N vector
%       shift   scalar, unit in samples
%   OUT:
%       y       N by 1 vector
%
%   Copyright 2005-7-6, Yingying Zhang, The University of Michigan


% Check dimensions to ensure y orientation same as x.
[m,n] = size(x);      
if m == 1, x = x.'; end % x is now column vector

% Check sign of shift and separate the integer and fraction shifts
if shift >= 0
    shift_int = floor(shift);
    shift_frac = shift - shift_int; % it's tau
    ind = shift_int+1;
    xpad = [x(ind:end); x(1:ind)];
elseif shift <= 0
    shift_int = floor(shift);
    shift_frac = shift - shift_int; % it's tau
    ind = shift_int+1;
    xpad = [x((end+ind):end); x(1:(end+ind))];
end

% use conv for fractional shift: h(n) = [tau 1-tau]
h = [shift_frac; 1 - shift_frac];
yy = conv(xpad, h);

y = yy(2:end-1); % column

if m == 1, y = y.'; end


    