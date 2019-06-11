 function y = mod0(x,b)
%function y = mod0(x,b)
%	the usual mod function returns values in the range [0,b)
%	whereas this mod function returns values in the range [-b/2,b/2)
%	which frequently arises in signal processing problems!
%	Jeff Fessler

y = mod(x + b/2, b) - b/2;
