% Gnufft_test.m
% Test the Gnufft object (vs exact Gdsft)

%% create Gnufft class object
if 1 || ~isvar('A'), printm 'setup Gnufft_test'
	im plc 3 3
	if 1 % 2d
		N = [32 30];
		J = [6 7];
	%	N = [32 1]; J = [5 1];
		omega = linspace(0, 10*2*pi, 201)'; % crude spiral:
		omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);
		if im, im subplot 1, plot(omega(:,1), omega(:,2), '.'), end
	else % 3d
		N = [16 12 14];
		J = [6 7 5];
		tmp = linspace(0, 10*2*pi, 201)';
		omega = pi*[cos(tmp) sin(tmp)].*tmp(:,[1 1])/max(tmp);
		omega(:,3) = pi * tmp / max(tmp); % spiral cone
		if im, im subplot 1, plot3(omega(:,1), omega(:,2), omega(:,3), '.'), end
	end

	K = 2*N;
	if N(2) == 1
		args = {omega(:,1), N(1), J(1), K(1)}; omega(:,2) = 0;
	else
		args = {omega, N, J, K};
	end

%	mask = true(N); mask(1,1) = false;
% todo: mask!

	args = {args{:}, 'table', 2^10, 'minmax:kb'}; % test with table

	A = Gnufft(args);

	Ad = Gdsft(omega, N);
end


%% test save/load
if 0
	save('/tmp/A.mat', 'A')
return
end


%% test data
if 1 || ~isvar('x'), printm 'setup data'
	if length(N) == 2
		x = zeros(N);
		x(5:25,10:25) = 1;
		x(15:20,15:20) = 2;
		x(15,5) = 2;
		if N(2) == 1
			x = x(:,5,1);
		end
	else
		rng(0)
		x = rand(N);
	end
	im(2, x, '\x')

	yd = Ad * x;

	if length(N) == 2
		n1 = ([0:N(1)-1]/N(1) - 0.5)*2*pi;
		if N(2) == 1 % 1D case
			yd_g = interp1(omega(:,1), yd, n1);
		else
			n2 = ([0:N(2)-1]/N(2) - 0.5)*2*pi;
			[nn1 nn2] = ndgrid(n1, n2);
			yd_g = griddata(omega(:,1), omega(:,2), yd,  nn1, nn2);
			yd_g(isnan(yd_g)) = 0;
		end
		im(3, abs(yd_g), '$|\y_d|$'), cbar
	end
end


%% build gram
if 1, printm 'Gnufft gram'
	wi = [1:size(omega,1)]';
	T = build_gram(A, wi);
	y2 = T * x(:);
	y1 = A' * (wi .* (A * x(:)));
	max_percent_diff y1 y2
%	equivs(y1, y2)
prompt
end


%% compare forward
if 1, printm 'forward'
	yn = A * [x(:) x(:)]; % test with two
	yn = yn(:,1);

	if length(N) == 2
		if N(2) == 1 % 1D case
			yn_g = interp1(omega(:,1), yn, n1);
		else
			yn_g = griddata(omega(:,1), omega(:,2), yn,  nn1, nn2);
			yn_g(isnan(yn_g)) = 0;
		end

		im(4, abs(yn_g), '$|\y_g|$'), cbar
		im(5, abs(yd_g - yn_g), '$|\y_g - \y_d|$'), cbar
	end
	max_percent_diff yd yn
end


%% compare adjoint
if 1, printm 'adjoint'
	yb = ones(size(omega,1), 1);
	xd = Ad' * yb;
	xd = iembed(Ad, xd);
	xn = A' * [yb yb]; % test with two
	xn = xn(:,1);
	xn = iembed(A, xn);

	im(7, fftshift(abs(xd)), '|back dtft|'), cbar
	im(8, fftshift(abs(xn)), '|back nufft|'), cbar
	im(9, fftshift(abs(xn-xd)), '|back err|'), cbar
	max_percent_diff xd xn
end
