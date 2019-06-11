% Gtomo_nufft_test_adj.m
% Test the adjoint of the Gtomo_nufft object

% create tomo class object
if ~isvar('G'), printm 'Gtomo_nufft_test_adj'
	nx = 8;
	ny = 6;
	nb = 10;
	na = 4;
	ktype = {'dsft'};
	ktype = {'interp', {'kaiser'}};

	fan_arg = {}; % parallel
	fan_arg = { % fan beam code
		'dis_src_det', 10^9, ...
		'channel_offset', 0, ...
		'dis_iso_det', 10^3, ...
	};
	mask = true(nx, ny);
	mask(1) = 0;
	G = Gtomo_nufft(mask, [nb na], ...
		fan_arg{:}, ...
		'orbit', 360, ...
		'chat', 0, ktype{:});
%	im(real(G.dsft))
end

if 1
	y = zeros(nb, na);
	y(nb/2+1,2) = 1;
	x = G' * y(:);
	x = embed(x, mask);
	if im
		clf, im(221, abs(x), 'abs'), cbar
		im(222, real(x), 'real'), cbar
		im(223, imag(x), 'imag'), cbar
	end
prompt
end

if 1
	Fatrix_test_basic(G, mask)
%	tester_tomo2(G, mask)

	[A Aa dif] = test_adjoint(G);
	im(224, dif, 'A-Aa'''), cbar
	printf('max diff = %g%%', max_percent_diff(A, Aa'))
prompt
end
