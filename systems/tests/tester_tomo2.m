 function tester_tomo2(A1, mask, varargin)
%function tester_tomo2(A1, mask, [options])
%|
%| Test suite for a Fatrix-type 2D or 3D system object, including Gblock tests.
%|
%| in
%|	A1	Fatrix or fatrix2
%|	mask	logical [nx ny] or [nx ny nz]
%|
%| option
%|	'A2'	Fatrix	optional 2nd object for testing and comparisons
%|			(it is also put through all the same tests)
%|	'multi'	0|1	see Fatrix_test_basic()
%|	'halt'	0|1	see Fatrix_test_basic()
%|	'equiv_thresh'	instead of jf_equal, use equivs() with this threshold
%|			default: [], which means use jf_equal
%|	'nblock'	# of blocks for testing block version (default 2)
%|			set to 0 to prevent block test
%|
%| Copyright 2005-8-2, Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

arg.multi = true;
arg.A2 = [];
arg.halt = true;
arg.equiv_thresh = [];
arg.nblock = 2; % for block test
arg = vararg_pair(arg, varargin, 'subs', {'G2', 'A2'});

Fatrix_test_basic(A1, mask, 'multi', arg.multi, ...
	'halt', arg.halt, ...
	'name', inputname(1), 'caller', caller_name)

switch ndims(mask)
case 2
	[nx, ny] = size(mask);
	ig = image_geom('nx', nx, 'ny', ny, 'dx', 1);
	x = ellipse_im(ig, []);

case 3
	[nx, ny, nz] = size(mask);
	ig = image_geom('nx', nx, 'ny', ny, 'nz', nz, 'dx', 1);
	x = ellipsoid_im(ig, '');
	x = x .* mask;
end

if arg.nblock
	tester_tomo2_block(A1, mask, x, arg.nblock, arg.multi, arg.equiv_thresh)
end

if ~isempty(arg.A2)
	Fatrix_test_basic(arg.A2, mask, 'multi', arg.multi, ...
		'name', inputname(4), 'caller', caller_name)
	tester_tomo2_compare(A1, arg.A2, x)
end


%
% tester_tomo2_compare()
%
function tester_tomo2_compare(A1, A2, x)

y1 = A1 * x;
y2 = A2 * x;
my_compare(y1, y2, 'A*x')

x1 = A1' * y1;
x2 = A2' * y1;
equivs(x1, x2)

j = round(size(A1,2) / 2); % roughly a middle pixel
y1 = A1(:,[j j+1]);
y2 = A2(:,[j j+1]);
y1 = [y1(:,1) y1(:,2)]; % for fatrix2
y2 = [y2(:,1) y2(:,2)]; % for fatrix2
my_compare(y1, y2, 'A(:,j)')

% check A(:,:)
if 0 && size(x,1) < 100
	t1 = A1(:,:);
	t2 = A2(:,:);
	my_compare(t1, t2, '(:,:)');
%	mpd = max_percent_diff(t1,t2);
%	printf('A(:,:)	mpd %g', mpd)
%	if mpd/100 > 1e-6, error 'A(:,:)', end
end

printm('passed %s', A1.caller)


%
% tester_tomo2_block()
% now block version
%
function tester_tomo2_block(A1, mask, x, nblock, multi, equiv_thresh)

B1 = Gblock(A1, nblock);

% B * x
y1 = A1 * x;
y2 = B1 * x;
my_compare(y1, y2, 'B*x')

y0 = y1;
odim = A1.arg.odim;
na = odim(end);
switch length(odim)
case 3
	getp = @(x,i) x(:,:,i);
case 2
	getp = @(x,i) x(:,i);
otherwise
	fail 'not done'
end
% caution: handle single view 3D projectors [ns nt 1 nrep]
catp = @(a,b) cat(length(odim)+1, a, b);


% B' * y
x1 = A1' * y0;
x2 = B1' * y0;
my_compare(x1, x2, 'B''*y')

comp = @(t1, t2, str) my_compare2(t1, t2, str, equiv_thresh);

%
% block operations
%
for k=1:nblock
	ia = k:nblock:na;
	str = sprintf('B{%d}', k);

	% check B{k} * x
	t1 = A1 * x;
	t1 = getp(t1,ia);
	t2 = B1{k} * x;
	comp(t1, t2, 'B{k} * x')

	% B{k} * [x x]
	if multi
		t2 = B1{k} * stackup(x,x);
		comp(catp(t1, t1), t2, [str ' * [x x]'])
	end

	% check B{k} * x(mask)
	t2 = B1{k} * x(mask);
	t2 = reshape(t2, size(t1));
	comp(t1, t2, [str ' * x()'])

	% check B{k} * [x(mask) x(mask)]
	if multi
		t2 = B1{k} * [x(mask) x(mask)];
		comp([t1(:) t1(:)], t2, [str ' * [x() x()]'])
	end

	% check B{k}' * y()
	tmp = block_insert(ia, odim, getp(y0,ia));
	t1 = A1' * tmp(:);
	t2 = B1{k}' * col(getp(y0,ia));
	if isempty(equiv_thresh)
%		my_compare(t1, t2, [str ''' * y'])
		mpd = max_percent_diff(t1,t2);
		if mpd
			printf([str '''*y mpd %g'], mpd)
			if mpd/100 > 1e-6
				fail 'B{k}'' * y()'
%				keyboard
			end
		end
	else
		comp(t1, t2, [str ''' * y'])
	end
end


% my_compare()
function my_compare(t1, t2, arg)
%max_percent_diff(t1, t2, arg)
%try
	jf_equal(t1, t2)
%catch
%	warn('%s failed', arg)
%	keyboard
%end


function my_compare2(t1, t2, arg, equiv_thresh)
if ~isempty(equiv_thresh)
	equivs(t1, t2, 'thresh', equiv_thresh)
else
	try
		my_compare(t1, t2, arg)
	catch
	%	keyboard
		fail('testing %s', arg)
	end
end


% block_insert()
% for block or OS methods, make data array that is all zeros but at "ia"
function out = block_insert(ia, dims, data)
out = zeros(dims);
switch length(dims)
case 2
	out(:,ia) = data;
case 3
	out(:,:,ia) = data;
otherwise
	fail 'not done'
end
