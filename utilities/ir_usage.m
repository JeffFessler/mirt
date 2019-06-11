 function ir_usage(mfile_name)
%function ir_usage(mfile_name)
%|
%| Show "usage" of a given mfile in a screen-friendly way then fail with error.
%| If mfile_name is not provided, then infer from "caller_name.m"
%|
%| 2016-08-07 Jeff Fessler, University of Michigan

if nargin < 1
	mfile_name = caller_name;
	if isempty(mfile_name)
		help(mfilename), error(mfilename)
		return
	end
end

more on
help(mfile_name)
more off

%evalin('caller', 'error(''test'')')
error([mfile_name ' usage'])
%fail(mfile_name)
