% Gblur_test.m
% Test the Gblur object

if 0 % test small even-sized psf
	types = {'fft,same', 'conv,same', 'conv,per'};
	types = {'imfilter,mirror'};

	for ii=1:numel(types)
		B = Gblur(true(6,2), 'psf', [1 -2]', 'type', types{ii});
		im plc 1 3
		im(1, full(B)')
		im(2, full(B'))
		im(3, full(B)' - full(B')), drawnow
		fatrix2_tests(B)
		test_adjoint(B, 'complex', 1);
	end
return
end


if 1 % test 'conv,valid' option
	psf = [0 1 2 1 0; 1 2 4 3 1; 0 2 3 1 0];
	idim = [16 20];
	mask = true(idim);
	Acv = Gblur(mask, 'psf', psf, 'type', 'conv,valid');
	x = shepplogan(idim(1), idim(2), 1);
	im plc 2 2
	im(1, x)
	y1 = Acv * x;
	y2 = conv2(x, psf, 'valid');
	jf_equal(y1, y2)
	y3 = Acv * x(:);
	x3 = Acv' * y3;
	x3 = reshape(x3, idim);
	fatrix2_tests(Acv)
%	test_adjoint(Acv, 'big', 1e-11, 'complex', 1)
	test_adjoint(Acv);
%return
end


if ~isvar('A'), printm 'setup Gblur_test'
	psf = [0 1 2 1 0; 1 2 4 3 1; 0 2 3 1 0];
	psf = psf / sum(psf(:));
	idim = [64 70];

	mask = true(idim);

	types = {'conv,same', 'conv,per', 'conv,valid', 'fft,same'};
	if 2 == exist('imfilter')
		types{end+1} = 'imfilter,same';
		types{end+1} = 'imfilter,circ';
%		types{end+1} = 'imfilter,mirror'; % adjoint fails
	end

	for ii=1:numel(types)
		A{ii} = Gblur(mask, 'psf', psf, 'type', types{ii});

		if isa(A, 'Fatrix')
			Fatrix_test_basic(A{ii}, mask, 'halt', 0) % paces
		else
			fatrix2_tests(A{ii})
		end
		test_adjoint(A{ii}, 'big', 1e-11, 'complex', 1)
	end

	im plc 3 3
	im(1, A{1}.arg.psf, 'psf'), cbar
	im(2, A{1}.arg.mask, 'mask'), cbar

	% compare conv and fft versions, for phantom that avoids image edges
	if 1
		x = shepplogan(idim(1), idim(2), 1);
		y{1} = A{1} * x;
		b{1} = A{1}' * y{1};
		im(3, x, 'x')
		im(4, y{1}, 'A * x')
		im(5, b{1}, 'A'' * y')
		for ii=2:numel(types)
			y{ii} = A{ii} * x;
			if ~streq(types{ii}, 'conv,valid')
				b{ii} = A{ii}' * y{1};
				equivs(y{ii}, y{1})
				equivs(b{ii}, b{1})
			end
		end
	prompt
	end
end


if 1 % time Fatrix vs fatrix2
	idim = [2^9 2^10];
	mask = true(idim);
	mask(1) = false; % test
	x = repmat(7*single(mask), [1 1 2^4]);
%	type = 'conv,same'; % same speed in this case
%	type = 'fft,same'; % fatrix2 faster in this case
	for ii=1:numel(types)
		type = types{ii};
		if streq(type, 'conv,valid')
			continue % no Fatrix version
		end
		A1 = Gblur(mask, 'psf', psf, 'class', 'Fatrix', 'type', type);
		A2 = Gblur(mask, 'psf', psf, 'class', 'fatrix2', 'type', type);
		A1(:,1); % warm up
		A2(:,1);
%profile on
		cpu etic
		y1 = A1 * x;
		cpu('etoc', ['Fatrix ' type ' time:'])
%profile report
%prompt
%profile on
		cpu etic
		y2 = A2 * x;
		cpu('etoc', ['fatrix2 ' type ' time: '])
		jf_equal(y1, y2)
	end
%profile report
end
