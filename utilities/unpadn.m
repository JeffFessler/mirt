  function out = unpadn(mat, newdim)
%|function out = unpadn(mat, newdim)
%| pad old to newdim, preserving 'center point'
if nargin == 1 && streq(mat, 'test'), unpadn_test, return, end
if nargin < 2, ir_usage, end

olddim = size(mat);
if any(newdim > olddim), error('must be smaller'), end

idx = cell(1, length(newdim));
for ii=1:length(newdim)
	pad = olddim(ii) - newdim(ii);
	if ~rem(pad,2) % even
		offset = pad/2;
	else
		if rem(olddim(ii),2) % odd
			offset = floor(pad/2);
		else
			offset = ceil(pad/2);
		end
	end
	idx{ii} = [1:newdim(ii)] + offset;
end
out = mat(idx{:});

function unpadn_test

jf_equal(unpadn([0 1 2 1 0], [1 3]), [1 2 1])
jf_equal(unpadn([0 1 2 1 0], [1 4]), [0 1 2 1])
jf_equal(unpadn([0 1 2 1], [1 3]), [1 2 1])
jf_equal(unpadn([0 0 1 2 1 0], [1 4]), [0 1 2 1])
