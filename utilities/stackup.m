 function ss = stackup(varargin)
%function ss = stackup(x1, x2, ...)
%|
%| Stack up 2D arrays to make 3D array, or stack up 3D arrays to make 4D array.
%|
%| This generalizes how [a; b] "stacks up" 1D vectors to make a 2D array.
%| This is useful in conjunction with stackpick().
%| It is akin to cat(ndims(x1)+1, x1, x2, ...)
%|
%| Copyright 2005-6-18, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if nargin == 1 && streq(varargin{1}, 'test'), stackup_test, return, end

arg1 = varargin{1};
switch ndims(arg1)
case 2
	% special usage: stackup(x, 'n3', n3)
%|function ss = stackup(x1, 'n3', n3) is like "op rep"
	if length(varargin) == 3 && streq(varargin{2}, 'n3')
		fail 'no longer supported'
		% fix: redo with repmat?
		n3 = varargin{3};
		ss = zeros([size(arg1) n3], class(arg1));
		for i3=1:n3
			ss(:,:,i3) = arg1;
		end
	return
	end

	% 2d stackup, allowing some of the others to be 3d
	% fix: refine
	nz = 0;
	for ii=1:length(varargin)
		varargin{ii} = squeeze(varargin{ii});
		nz = nz + size(varargin{ii},3);
	end
	ss = zeros([size(arg1) nz], class(arg1));
	iz = 0;
	for ii=1:length(varargin)
		nz = size(varargin{ii},3);
		ss(:,:,iz+[1:nz]) = varargin{ii};
		iz = iz + nz;
	end

case 3

	ss = zeros([size(arg1) length(varargin)], class(varargin{1}));
	for ii=1:length(varargin)
		ss(:,:,:,ii) = varargin{ii};
	end

otherwise
	fail 'only stacking 2d -> 3d and 3d -> 4d done'
end


% stackup_test
function stackup_test
x1 = ones(5,3);
x2 = 2*x1;
y1 = stackup(x1, x2);
jf_equal(y1, cat(3, x1, x2))

%y1 = stackup(x1, 'n3', 5);
%jf_equal(y1, repmat(x1, [1 1 5]))

x1 = ones(5,4,3);
x2 = 2*x1;
y1 = stackup(x1, x2);
jf_equal(y1, cat(4, x1, x2))
