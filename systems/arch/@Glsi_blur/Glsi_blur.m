 function ob = Glsi_blur(nx, ny, psf, mask, chat)
%function ob = Glsi_blur(nx, ny, psf, mask, chat)
%	Construct Glsi_blur object, which does 2d convolution
%	See Glsi_blur_test.m for example usage.
%	Basically, you create a system matrix object by calling:
%		G = Glsi_blur(nx, ny, psf)
%	and then you can use it thereafter by typing commands like
%		y = G * x;
%	which will auto-magically perform the convolution.
%	So you can represent convolution by a matrix multiplication
%	in programs, rather than calling a subroutine.
%	Why bother?  Because it is easier to debug, and more general,
%	to write the routines using matrices in the first place.
%
%	Besides simple utilities like display, there are the following
%	capabilities of this object:
%	y = G * x		convolution
%	? x = G' * y		"transposed" convolution
%	? y = Gt(:,ii)' * x	partial convolution
%	? x = Gt(:,ii) * y	partial convolution
%

if ~isvar('chat') | isempty(chat)
	chat = 0;
end

%
%	default object
%
ob.nx = 0;	% input image size
ob.ny = 0;
ob.nb = 0;	% output image size
ob.na = 0;
ob.psf = [];	% the psf
ob.chat = chat;
ob.power = 1;
ob.mask = [];
%ob.ia_start = 0;	% change for ordered subsets if Gt(:,ii)
%ob.ia_inc = 1;		% change for ordered subsets if Gt(:,ii)
ob.is.empty	= true;
ob.is.transpose = false;
ob.is.transpose_after_sub = false;
ob.is.subref = false;
ob.is.masked = false;	% set to 1 if G(:,mask(:))

if nargin < 3
	warning 'Glsi_blur called with too few arguments'
	help Glsi_blur
	ob = class(ob, 'Glsi_blur');
	return
end

if nargin > 5
	help Glsi_blur
	error nargin
end


	%
	%	input arguments
	%
	ob.nx = nx;
	ob.ny = ny;
	ob.psf = psf;
	if isvar('mask') & ~isempty(mask)
		if any(size(mask) ~= [nx ny]), error 'mask size', end
		ob.mask = mask;
	else
		ob.mask = true(nx,ny);
	end
	ob.nb = nx;
	ob.na = ny;

	ob.is.empty = false;
	ob = class(ob, 'Glsi_blur');
