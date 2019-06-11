 function ss = stack(varargin)
%function ss = stack(x1, x2, ...)
%function ss = stack(x1, 'n3', n3) is like "op rep"
% stack up 2D arrays to make 3D array
% this generalizes how [a; b] "stacks up" 1D vectors to make a 2D array.

if nargin < 1, ir_usage, end

arg1 = varargin{1};
if ndims(arg1) == 2
	% special usage: stack(x, 'n3', n3)
	if length(varargin) == 3 && streq(varargin{2}, 'n3')
		n3 = varargin{3};
		ss = zeros([size(arg1) n3]);
		for i3=1:n3
			ss(:,:,i3) = arg1;
		end
	return
	end

	ss = zeros([size(arg1) length(varargin)]);
	for ii=1:length(varargin)
		ss(:,:,ii) = varargin{ii};
	end
else
	error 'only stacking 2D to make 3D done'
end
