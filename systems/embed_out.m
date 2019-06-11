  function [y, eo_fun] = embed_out(y, Mdim)
%|function [y, eo_fun] = embed_out(y, Mdim)
%|
%| Dual or adjoint of embed_in().  See embed_in() for details.
%| To be used around _back() routines that expect [(M) *L]
%| input and produce [(N) *L] output that then must be reshaped
%| to be either [(N) (L)] or [np *L]
%|
%| in
%|	y	[*M (L)] or [(M) (L)]	input data array(s), possibly as columns
%|	Mdim				(M)
%| out
%|	y	[(M) *L]		output data arrays
%|	eo_fun	strum object with methods:
%|	x = eo_fun.shape(x, mask, np)	reshape x from [(N) *L] to be either
%|					[np (L)] or [(N) (L)], depending on y
%|
%| Copyright 2006-12-9, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(y, 'test'), embed_out_test, clear y, return, end
if nargin < 2, Mdim = size(y,1); end % trick: try to guess...

if size(y,1) == prod(Mdim) % convert [*M (L)] to [(M) *L]
	state.column = true;
	dimy = size(y);
	state.diml = dimy(2:end); % [(L)]
	y = reshape(y, [Mdim prod(state.diml)]); % [(M) *L]

else % convert [(M) (L)] to [(M) *L]
	state.column = false;
	dimi = [size(y) 1];
	jf_equal(dimi(1:length(Mdim)), Mdim)
	state.diml = dimi(length(Mdim)+1:end); % (L) (possibly empty)
	y = reshape(y, [Mdim prod(state.diml) 1]); % [(N) *L]
end

eo_fun = strum(state, {'shape', @embed_out_shape});


%
% embed_out_shape()
%
function x = embed_out_shape(state, x, mask, np)

if nargin < 4, np = sum(mask(:)); end

diml = state.diml;

if state.column % column in yields column out, i.e., [(N) *L] to [np (L)]
	if any(diml > 1)
		x = reshape(x, numel(mask), []); % [*N *L]
		x = x(mask,:); % [np *L]
		x = reshape(x, [np diml]); % [np (L)]
	else
		x = x(mask); % [np,1]
	end

else % [(N) *L] to [(N) (L)]
	if any(diml > 1)
		x = reshape(x, [size(mask) diml]); % [(N) (L)]
	end
end


%
% embed_out_test
%
function embed_out_test
ig = image_geom('nx', 10, 'ny', 8, 'dx', 1);
ig.mask = ig.circ > 0;

Mdim = [6 7];
y0 = ones(Mdim); % [(M)]
dl = [2 3];

[y2 eo] = embed_out(y0(:), Mdim); % single column
x2 = eo.shape(ig.ones, ig.mask);
jf_equal(size(x2), [ig.np 1])

[y2 eo] = embed_out(y0, Mdim); % usual 2d
x2 = eo.shape(ig.ones, ig.mask);
jf_equal(size(x2), [ig.nx ig.ny])

yr = repmat(y0, [1 1 dl]); % [(M) (L)]
[y2 eo] = embed_out(yr, Mdim);
x2 = eo.shape(ones([ig.nx ig.ny prod(dl)]), ig.mask);
jf_equal(size(x2), [ig.nx ig.ny dl])

yr = repmat(y0(:), [1 dl]); % [*M (L)]
[y2 eo] = embed_out(yr, Mdim);
x2 = eo.shape(ones([ig.nx ig.ny prod(dl)]), ig.mask);
jf_equal(size(x2), [ig.np dl])
