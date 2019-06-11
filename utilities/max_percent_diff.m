  function d = max_percent_diff(s1, s2, varargin)
%|function d = max_percent_diff(s1, s2, [options])
%|
%| compute the "maximum percent difference" between two signals: s1, s2
%| options
%|	1		use both arguments as the normalizer
%|	string		print this
%|
%| Copyright 2000-9-16, Jeff Fessler, University of Michigan

if nargin == 1 && streq(s1, 'test'), max_percent_diff_test, return, end
if nargin < 2, ir_usage, end

base = '';

if ischar(s1)
	t1 = s1;
	s1 = evalin('caller', t1);
	base = [caller_name ': '];
else
	t1 = inputname(1);
end
if ischar(s2)
	t2 = s2;
	s2 = evalin('caller', t2);
else
	t2 = inputname(2);
end

use_both = 0;
doprint = (nargout == 0);
while (length(varargin))
	arg1 = varargin{1};
	if isnumeric(varargin{1})
		use_both = 1;
		varargin = {varargin{2:end}};
		continue
	end

	if ischar(arg1)
		doprint = 1;
		base = [arg1 ': '];
		varargin = {varargin{2:end}};
		continue
	end
end


% first check that we have comparable signals!
if ~isequal(size(s1), size(s2))
	printm(['size(%s) = %s'], inputname(1), mat2str(size(s1)))
	printm(['size(%s) = %s'], inputname(2), mat2str(size(s2)))
	error 'dimension mismatch'
end

if any(isnan(s1(:)))
	warn([mfilename ': NaN values in %s in %s!?'], t1, caller_name)
end
if any(isnan(s2(:)))
	warn([mfilename ': NaN values in %s in %s!?'], t2, caller_name)
end
s1 = doubles(s1);
s2 = doubles(s2);

if use_both
	denom = max(abs([s1(:); s2(:)]));
	if ~denom
		d = 0;
	else
		d = max(abs(s1(:)-s2(:))) / denom;
	end
else
	denom = max(abs(s1(:)));
	if ~denom
		denom = max(abs(s2(:)));
	end
	if ~denom
		d = 0;
	else
		d = max(abs(s1(:)-s2(:))) / denom;
	end
end

d = full(d) * 100;

if doprint
	printf([base 'max_percent_diff(%s, %s) = %g%%'], t1, t2, d)
	if ~nargout
		clear d
	end
end


%
% max_percent_diff_test
%
function max_percent_diff_test
v1 = [0 1000];
v2 = [0 1001];
max_percent_diff([0 1000], [1 1000])
max_percent_diff(v1, v2, 'numeric')
max_percent_diff '[0 1000]' '[1 1000]'
