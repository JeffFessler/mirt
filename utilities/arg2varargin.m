function [varargin] = arg2varargin(arg)
% converts a variable input arguments structure, arg,  into a cell of 
% fieldnames and inputs to be passed along to another function as a 
% varargin cell, i.e. [out] = myfunc(in1, in2, varargin{:})
%
% IMPORTANT NOTE! when passing the output of this function, varargin, to
% the next fucntion, the "{:}" at the end is necessary
%
% Inputs:
%   arg     struct of variable inputs
%
% Outputs:
%   vargin  cell of fieldname strings and inputs
%
% Written 2021-06-03, Melissa Haskell, University of Michigan

varargin = [fieldnames(arg), struct2cell(arg)]';


end

