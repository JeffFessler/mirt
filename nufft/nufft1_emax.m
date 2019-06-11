 function [emax, err] = nufft1_emax(varargin)
%function [emax, err] = nufft1_emax(varargin)
%
% Compute exponential approximation error for each input frequency for 1D NUFFT
% and return worst-case.
%
% in
%	everything is passed to nufft1_build()
% out
%	emax	[1,1]	worst-case error over all frequencies specified
%
% Copyright 2006-4-17, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

st = nufft1_build(varargin{:});
%pr sum(st.sn)-175.3399
err = st.err;
emax = max(err);
