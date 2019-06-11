  function [mat sp] = aspire_buff2mat(buff)
%|function [mat sp] = aspire_buff2mat(buff)
%|
%| convert 'buffer' read from .wtf to matlab sparse matrix
%| primarily for testing
%|
%| in
%|	buff	byte	read using wtfmex()
%|
%| out
%|	mat	sparse	matlab sparse matrix, full sized
%|	sp	struct	header information etc.
%|
%| Copyright 2008-9-23, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(buff, 'test'), aspire_buff2mat_test, return, end

sp = aspire_buff_parse(buff); % parse buffer
mat = aspire_buff_sp2mat(sp); % convert to matrix


% aspire_buff_sp2mat()
% make sparse matrix from parsed buffer
function mat = aspire_buff_sp2mat(sp)

switch sp.group_by
case 0 % by_row
	ii = [];
	for kk=1:sp.ngroup
		ii = [ii; repmat(kk, double(sp.length(kk)), 1)];
	end
	jj = 1 + sp.index;
case 1 % by_col
	ii = 1 + sp.index;
	jj = [];
	for kk=1:sp.ngroup
		jj = [jj; repmat(kk, double(sp.length(kk)), 1)];
	end
otherwise
	fail 'bug'
end

% matlab stupidly insists on 'double' arguments!
mat = sparse(double(ii), double(jj), double(sp.value), ...
	double(sp.nb*sp.na), double(sp.nx*sp.ny), double(sp.nwt));



%
% aspire_buff2mat_test
%
function aspire_buff2mat_test

ig = image_geom('nx', 22, 'ny', 20, 'dx', 2);
sg = sino_geom('par', 'nb', 24, 'na', 18, 'dr', 1.8);
%ig.mask = ig.circ(ig.fov/2) > 0;

sw = sg.dr * 2;
As = Gtomo2_strip(sg, ig, 'strip_width', sw);

pair = aspire_pair(sg, ig, 'support', 'array', 'strip_width', sw);
types = {'col', 'row'};
for it=1:length(types)
	buff = wtfmex('asp:gensys', pair', types{it}, uint8(ig.mask), int32(0));

	[mat sp] = aspire_buff2mat(buff);
	jf_equal(sp.mask, ig.mask, 'accept_logical_eq_uint8', true) % trick

%	tmp = As.arg.G; % Gsparse
	tmp = As.arg.matrix; % Gmatrix
	equivs(mat, tmp)
end

if 0
	pr sp
	im pl 2 2
	im(1, mat)
	im(2, tmp)
	im(3, mat - tmp)
	keyboard
end
