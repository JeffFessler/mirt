 function out = ndgrid_jf(otype, varargin)
%function out = ndgrid_jf('cell', varargin)
%function out = ndgrid_jf('mat', varargin)
%| version of ndgrid where output is cell array (default) or large matrix
%| think:
%| [out{1} ... out{M}] = ndgrid_cell(in{1}, ..., in{M});
%| also supported is ndgrid_jf('mat'|'cell', {in1, in2, ...})
%|
%| fix: is there a simple way to do this with nargout?
%| Copyright 2007, Jeff Fessler, University of Michigan

if nargin == 1 && streq(otype, 'test'), ndgrid_jf_test, return, end
if nargin < 2, ir_usage, end

switch otype
case 'cell'
	is_cell = 1;
case 'mat'
	is_cell = 0;
otherwise
	error 'unknown output type'
end

% handle a convenient special case (this avoids strum.arg{:})
if length(varargin) == 1 && iscell(varargin{1})
	out = ndgrid_jf(otype, varargin{1}{:});
return
end

if length(varargin) == 1
	if is_cell
		out = {varargin{1}};
	else
		out = varargin{1}(:);
	end
return
end

nn = length(varargin);

for ii=nn:-1:1
	varargin{ii} = full(varargin{ii});
	siz(ii) = numel(varargin{ii});
end

out = cell(1,nn);
for ii=1:nn
	x = varargin{ii}(:);
	s = siz; s(ii) = []; % remove i-th dimension
	x = reshape(x(:,ones(1,prod(s))), [length(x) s]); % expand x
	out{ii} = permute(x, [2:ii 1 ii+1:nn]); % permute to i'th dimension
end

if ~is_cell
	out = stackup(out{:});
end

%
% ndgrid_jf_test
%
function ndgrid_jf_test
x = 1:2;
y = 1:3;
z = 1:4;
out = ndgrid_jf('cell', x, y, z);
[xx yy zz] = ndgrid(x, y, z);
if ~isequal(xx, out{1}) || ~isequal(yy, out{2}) || ~isequal(zz, out{3})
	error 'bug'
end

out = ndgrid_jf('mat', x, y, z);
if ~isequal(xx, out(:,:,:,1)) ...
	|| ~isequal(yy, out(:,:,:,2)) ...
	|| ~isequal(zz, out(:,:,:,3))
	error 'bug'
end

out2 = ndgrid_jf('mat', {x, y, z});
jf_equal(out, out2)

out = ndgrid_jf('cell', x);
out = ndgrid_jf('mat', x);
