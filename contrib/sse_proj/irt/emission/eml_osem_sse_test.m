% eml_osem_sse_test.m
% compare aspire and matlab versions of E-ML-OSEM
% Copyright 2001-07-23, Jeff Fessler, The University of Michigan

%
% generate data
%
if ~isvar('xfbp'), printm 'data'
	if has_aspire
		f.dir	= test_dir;
		f.wtf	= [f.dir 't,g.wtf'];
		f.wtr	= strrep(f.wtf, 'wtf', 'wtr');
		f.yi	= [f.dir 'yi.fld'];
		f.ci	= [f.dir 'ci.fld'];
		f.ri	= [f.dir 'ri.fld'];
		f.mask	= [f.dir 'mask.fld'];
	end
	em_test_setup; f.is3b = 0;
%	em3_test_setup; f.is3b = 1;
prompt
end


%
% block object for block iterative
%
if ~isvar('Gb'), printm 'Gb'
	f.nblock = 2;
	Gb = Gblock(G, f.nblock, 0);
end


%
% matlab iterations
%
if ~isvar('xmat'), printm 'matlab E-ML-OSEM'
	f.niter = 8;
	f.pixmax = 6;
	f.precon = 'fast';
	f.precon = 'classic';

	if f.is3b % fix: needs work!
		shaper = sprintf('reshape(proj, [%d %d])', n.n1*n.n2, n.a);
		shaper = inline(shaper, 'proj');
	else
		shaper = inline('x', 'x');
	end
	xinit = ig.ones;
	[xmat precon] = eml_osem_sse(xinit(ig.mask), Gb, ...
			shaper(yi), shaper(ci), shaper(ri), ...
			'niter', f.niter, 'pixmax', f.pixmax, ...
			'precon', f.precon);
	im clf, im(121, xmat, 'Matlab E-ML-OSEM iterations'), cbar horiz
prompt
end
