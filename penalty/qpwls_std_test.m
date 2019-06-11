% qpwls_std_test
% examine FFT-based prediction of QPWLS estimate pixel variance
% particularly w.r.t. effect of "mask"

if ~isvar('G2')
	down = 1;
	ig1 = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
	ig2 = ig1;
	ig2.mask = ig2.circ(220) > 0;
	sg = sino_geom('ge1', 'down', down);

	% now a double sized version
	ig3 = image_geom('nx', ig1.nx*2, 'ny', ig1.nx*2, 'fov', ig1.fov*2);

	G1 = Gtomo_dd(sg, ig1);
	G2 = Gtomo_dd(sg, ig2);
	G3 = Gtomo_dd(sg, ig3);
	im pl 2 1
	im(1, ig1.mask, 'mask1') % square
	im(2, ig2.mask, 'mask1') % circle
prompt
end

if ~isvar('std2')
	l2bs = [5:14];

	for ii=1:length(l2bs)
		l2b = l2bs(ii);
		R1 = Robject(ig1.mask);
		R2 = Robject(ig2.mask);
		R3 = Robject(ig3.mask);

		[psf var fw1(ii)] = qpwls_psf(G1, R1, 2^l2b, ig1.mask);
		std1(ii) = sqrt(var);
		[psf var fw2(ii)] = qpwls_psf(G2, R2, 2^l2b, ig2.mask);
		std2(ii) = sqrt(var);
		[psf var fw3(ii)] = qpwls_psf(G3, R3, 2^l2b, ig3.mask);
		std3(ii) = sqrt(var);
	end
end

if im
	clf
	plot(l2bs, std2, 'o', l2bs, std1, 's', l2bs, std3, '+')
	xlabel 'log_2(\beta)'
	ylabel 'relative standard deviation'
	legend('circle mask', 'full mask', 'double size')
	ir_savefig fig_qpwls_std_fan
end
