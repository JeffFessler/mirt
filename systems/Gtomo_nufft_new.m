  function ob = Gtomo_nufft_new(sg, ig, varargin)
%|function ob = Gtomo_nufft_new(sg, ig, options)
%|
%| Construct Gtomo_nufft object, for Fourier-based forward and backprojection.
%| Currently handles 2D parallel and fan-beam cases.
%|
%| See Gtomo_nufft_test.m for example usage.
%|
%| Basically, you create an object calling:
%|		A = Gtomo_nufft(...)
%| and then you can use it thereafter by typing commands like
%|		y = A * x;
%| which will auto-magically do the multiplication by calling a mex file.
%|
%| in
%|	sg	strum	see sino_geom()
%|	ig	strum	see image_geom()
%|
%| optional arguments for which a value must follow:
%|	'chat'		verbose printing of debug messages
%|	'xscale'	use -1 to flip in x direction
%|	'yscale'	use -1 to flip in y direction
%|
%|	'beam'		structure for beam PSF (see below)
%|	'basis'		structure for image-domain basis (see below)
%|
%|	'is.dsft1'	use full 1D DSFT rather than 1D NUFFT
%|	'is.dsft2'	use full 2D DSFT rather than 2D NUFFT
%|			(these are expensive but useful for debugging)
%|	'is.shift0'	apply shifts so center of rotation is at 0 (is default)
%|	'is.complex'	return complex projection and adjoint (unusual)
%|
%| interpolation specification
%|
%|	'J'		neighborhood size (default is 6)
%|	'Kd'		# of FFT points along each dimension 2*[nx ny]
%|	'Krho'		number of radial samples in polar space, default is nb.
%|
%|	The default (recommended) interpolator uses minmax with KB scaling.
%|	These others below are provided for making comparisons.
%|	'is.kaiser'	use minmax tuned kaiser-bessel interpolator
%|	'interp'	{ cell array of arguments to pass to nufft_init() }
%|			(This allows complete user control of interpolator.)
%|	Example: straw-man bilinear interpolation is:
%|		..., 'interp', {'linear'}, 'J', 2, ...
%|	Example: for straw-man uniform scaling factors in NUFFT:
%|		..., 'interp', {'minmax:user', {{1}, {1}}, {{0.5}, {0.5}}}
%|	Example: for table-based NUFFT for large-scale problems:
%|		..., 'interp', {'table', 2^11, 'minmax:kb'}
%|
%| out
%|	ob [nd,np]	np = sum(mask(:)), so it is already "masked"
%|
%| For more help, see Gtomo_nufft_test.m
%|
%| Based on this paper:
%| Yingying Zhang-O'Connor, J A Fessler, IEEE T-MI 25(5):582-9, May 2006
%| "Fourier-based forward and back-projectors in iterative fan-beam
%| tomographic image reconstruction"
%| doi 10.1109/TMI.2006.872139
%|
%| Copyright 2005-1-19, Jeff Fessler, University of Michigan
%| Fan-beam contributed by Yingying Zhang.  Basis etc. by Samuel Matej.

if nargin == 1 && streq(sg, 'test'), run_mfile_local('Gtomo_nufft_test'), return, end
if nargin < 1, help(mfilename), error(mfilename), end

if streq(sg.type, 'fan')
	if sg.dfs ~= 0
		fail 'only arc detector implemented!'
	end

	fanarg = { ...
		'dis_src_det', sg.dsd, ...
		'dis_iso_det', sg.dod, ...
		'source_offset', sg.source_offset, ...
		'offset_s', sg.offset_s, ...
		'ds', sg.ds
	};
else
	fanarg = { ...
		'ds', sg.dr, ...
		'offset_s', sg.offset_r
	};
end

ob = Gtomo_nufft(ig.mask, [sg.nb sg.na], ...
	'dx', ig.dx, ...
	'yscale', -ig.dy / ig.dx, ...
	'orbit', sg.orbit, ...
	'orbit_start', sg.orbit_start, ...
	'strip_width', sg.strip_width, ...
	fanarg{:}, ...
	varargin{:});
