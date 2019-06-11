 function y = stackpick(x, ii)
%function y = stackpick(x, ii)
%| pick one (or more) of the last "slices" of a multi-dimensional array
%| think x(:,:,...,:,ii)
%| This is useful in conjunction with stackup().
%| Copyright 2005-6-18, Jeff Fessler, University of Michiga

if nargin == 1 && streq(x, 'test'), stackpick_test, return, end
if nargin < 2, ir_usage, end

s = size(x);
if any(ii < 1) || any(ii > s(end))
	pr ii
	fail('bad index! last dimension is %d', s(end))
end

% an alternative way to do this would be
% index = repmat({':'}, 1, ndims(x));
% index{end} = ii;
% y = x(index{:});

%if length(s) == 2 && s(2) == 1 % [M,1] case; ambiguous but treat as [M]
%	if ii > 1, fail('[%d,1] input with ii=%d', size(x,1), ii), end
%	y = x;
%	warn 'ambiguous, resolving as if 1D'
%	todo: add an option not to warn?
%else
	x = reshapee(x, [], s(end));
	y = x(:,ii);
	y = reshape(y, [s(1:end-1) length(ii)]);
%end

% stackpick_test
function stackpick_test
x = reshape(1:(5*3*4), [5 3 4]);
y = stackpick(x, 2);
jf_equal(y, x(:,:,2))

% 1d
jf_equal(stackpick([1:9], 2), 2) % [1,M]
jf_equal(stackpick([1:9]', 1), [1:9]') % [M,1]

% 2d
x = reshape(1:(4*3), [4 3]);
y = stackpick(x, 2);
jf_equal(y, [5:8]')

% multiple index
y = stackpick(x, [1 3]);
jf_equal(y, [[1:4]' [9:12]'])
