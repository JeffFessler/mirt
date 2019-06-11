 function y = ir_imfilter_many(x, psf, varargin)
%function y = ir_imfilter_many(x, psf, varargin)
%|
%| as of 2015-04-12, octave's imfilter can handle only 2D input
%| whereas matlab can handle 3D inputs
%| this is a kludge to deal with that.
%|
%| 2015-04-12, Jeff Fessler, University of Michigan

if nargin < 2, ir_usage, end

if ~ir_is_octave || ndims(x) == ndims(psf)
	y = imfilter(x, psf, varargin{:}); % matlab or usual sizes
return
end

% hereafter for octave with unequal sizes

if ndims(x) == 3 && ndims(psf) == 2
	tmp = imfilter(x(:,:,1), psf, varargin{:});
	nz = size(x,3);
	y = zeros([size(tmp) nz]);
	y(:,:,1) = tmp;
	for iz = 2:nz
		y(:,:,iz) = imfilter(x(:,:,iz), psf, varargin{:});
	end
return
end

pr size(x)
pr size(psf)
fail('imfilter in octave does not support such sizes, nor does my kludge')
