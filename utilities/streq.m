  function tf = streq(a, b, n)
%|function tf = streq(a, b [,n])
%|
%| return 1 if two strings "a" and "b" are equal
%| (optionally checking only up to 1st n chars)
%| caution: whereas strcmp allows comparisons of cell arrays,
%| this routine allows only two strings.
%|
%| Jeff Fessler

if nargin == 1 && strcmp(a, 'test'), streq_test, return, end
if nargin < 2, ir_usage, end

if ~ischar(a) || ~ischar(b), tf = false; return, end

if nargin == 2
	tf = strcmp(a,b);
elseif nargin == 3
	tf = strncmp(a,b,n);
else
	error(mfilename)
end

function streq_test
jf_equal(true, streq('ok', 'okay', 2))
jf_equal(false, streq('yes', 'no'))
jf_equal(false, streq(7, ' '))
