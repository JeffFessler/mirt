 function A = ir_convmtx1(h, n, mtype)
%function A = ir_convmtx1(h, n, mtype)
%|
%| Version of convmtx that supports 'same' option.
%| 2017-03-13, Jeff Fessler, University of Michigan

if nargin == 1 && streq(h, 'test'), ir_convmtx1_test, return, end
if nargin < 2, ir_usage, end

if nargin < 3, mtype = ''; end

if ndims(h) > 2 || min(size(h)) ~= 1
	fail('1D h')
end

switch mtype
case ''
	A = convmtx(h, n);

case 'convn'
	A = convn(eye(n), h(:)'); % trick

case 'same'
	A = convn(eye(n), h(:)', 'same'); % trick

case 'old_way'

	h = h(:); % 1D
	%hc = [h; zeros(n-length(h), 1, class(h))];
	hc = [h; zeros(n-1, 1, class(h))];
	%h = [zeros(length(h)-1, 1, class(h)); h];
	hr = [h(1); zeros(n-1, 1, class(h))];
	A = toeplitz(hc, hr)';

%{ for 'same' option something like this:
	trim1 = floor(length(h)/2);
	trim2 = floor((length(h)-1)/2);
	A = A(:, (trim1+1):(end-trim2));
%}

otherwise
	fail('bad mtype "%s"', mtype)
end


function ir_convmtx1_test
nlist = {8, 9}; % check odd,even
hlist = {1:4, 1:5};
for ii = 1:numel(nlist)
	n = nlist{ii};
	for jj = 1:numel(hlist)
		h = hlist{jj};

		A0 = convmtx(h, n);
		A1 = ir_convmtx1(h, n);
		jf_equal(A0, A1)
		A3 = ir_convmtx1(h, n, 'convn');
		jf_equal(A0, A3)

		A2 = ir_convmtx1(h, n, 'same');
	%	size(A2)
	end
end
