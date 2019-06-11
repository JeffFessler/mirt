% gg_Gnufft_test.m
% Test the Gnufft object (vs exact Gdsft)

% create Gnufft class object
% lets try a 1-d comparison instead
npow = [5:10];
%npow = 8*ones(10);
dispimgs = 0;
doexact = 1;
%npow = 10;
%Jval = [3:12];

for ii = 1:length(npow)
%for ii = 1:length(Jval)

	disp 'setup Gnufft_test'
	N = [2^npow(ii) 2^npow(ii)];
	%N = [2^npow 2^npow];
	J = [6 6];
	%J = [Jval(ii) Jval(ii)];
	% N = [32 1]; J = [5 1];
	K = 2*N;
	%omega = linspace(0, 10*2*pi, 4^npow(ii))'; % crude spiral:
	omega = linspace(0, 10*2*pi, 201)';
	omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);
	pl = 330;
	%if im, clf, subplot(pl+1), plot(omega(:,1), omega(:,2), '.'), end
	if N(2) == 1
		args = {omega(:,1), N(1), J(1), K(1), N(1)/2}; omega(:,2) = 0;
	else
		%args = {omega, N, J, K, N/2, 'table', 2^10, 'minmax:kb'};
		args = {omega, N, J, K, N/2};
	end
	% mask = true(N); mask(1,1) = false;
	% todo: mask!
	G = Gnufft(args);

	% set up greengard nufft
	Msp = 6;R = 2;
	%Msp = Jval(ii);R = 2;n_shift = [];
	gg_G = gg_Gnufft(omega,N(1),Msp,R);
	if doexact
		Gd = Gdsft(omega, N, 'n_shift', N/2);
	end
	%
	% test data
	%
	disp 'setup data'
	if N(1) == 32
		x = zeros(N);
		x(5:25,10:25) = 1;
		x(15:20,15:20) = 2;
		x(15,5) = 2;
	else
		x = zeros(N);
		x(5*floor(N(1)/32):25*floor(N(1)/32),10*floor(N(1)/32):25*floor(N(1)/32)) = 1;
		x(15*floor(N(1)/32):20*floor(N(1)/32),15*floor(N(1)/32):20*floor(N(1)/32)) = 2;
		x(15*floor(N(1)/32),5*floor(N(1)/32)) = 2;
	end

	if N(2) == 1
		x = x(:,5);
	end
	if dispimgs,im(pl+1, x, 'x'),end

	%%%%% forward exact dft %%%%%%%%%
	if doexact
		if length(npow) > 1 || ii == 1 % only need to do it once if only
																	% one data length
			tic
				yd = Gd * x(:);
			dtfttime(ii) = toc;
		end
	else
		yd = zeros(length(omega),1);
	end


	if dispimgs
	n1 = ([0:N(1)-1]/N(1) - 0.5)*2*pi;
	if N(2) == 1 % 1D case
		yd_g = interp1(omega(:,1), yd, n1);
	else
		n2 = ([0:N(2)-1]/N(2) - 0.5)*2*pi;
		[nn1 nn2] = ndgrid(n1, n2);
		yd_g = griddata(omega(:,1), omega(:,2), yd,	nn1, nn2);
		yd_g(isnan(yd_g)) = 0;
	end
	im(pl+2, abs(yd_g), '|y_d|'), colorbar
	end

	% compare forward
	tic
	yn = G * x(:); % test with two
	nforwtime(ii) = toc;
	tic
	gg_yn = gg_G * x(:);
	gg_nforwtime(ii) = toc;

	if dispimgs
	if N(2) == 1 % 1D case
		yn_g = interp1(omega(:,1), yn, n1);
		gg_yn_g = interp1(omega(:,1), gg_yn, n1);
	else
		yn_g = griddata(omega(:,1), omega(:,2), yn,	nn1, nn2);
		yn_g(isnan(yn_g)) = 0;
		gg_yn_g = griddata(omega(:,1), omega(:,2), gg_yn,	nn1, nn2);
		gg_yn_g(isnan(gg_yn_g)) = 0;
	end

		im(pl+5, abs(yn_g), '|y_g|'), colorbar
		im(pl+6, abs(yd_g - yn_g), '|y_g - y_d|'), colorbar
		im(pl+8, abs(gg_yn_g), '|ggy_g|'), colorbar
		im(pl+9, abs(yd_g - gg_yn_g), '|ggy_g - y_d|'), colorbar
	end

	max_percent_diff(yd, yn, mfilename)
	nforwerr(ii) = nrmse(yd,yn,ones(size(yd)));
	max_percent_diff(yd, gg_yn, mfilename)
	gg_nforwerr(ii) = nrmse(yd,gg_yn,ones(size(yd)));

	% compare adjoint
	if dispimgs,figure,end
	yb = ones(size(omega,1), 1);
	%%%%%% back exact dft
	if doexact
		if length(npow) > 1 || ii == 1 % only need to do it once if only
																	% one data length
			tic
			xd = Gd' * yb;
			idtfttime(ii) = toc;

			xd = embed(xd, G.arg.mask);
		end
	else
		xd = zeros(size(G.arg.mask));
	end


	tic
	xn = G' * yb;
	nbacktime(ii) = toc;
	xn = embed(xn, G.arg.mask);
	tic
	gg_xn = gg_G' * yb;
	gg_nbacktime(ii) = toc;

	if dispimgs
		im(pl+1, abs(xd), '|back dtft|'), colorbar
		im(pl+4, abs(xn), '|back nufft|'), colorbar
		im(pl+5, abs(xn-xd), '|back err|'), colorbar
		im(pl+7, abs(gg_xn), '|back nufft| - greengard'), colorbar
		im(pl+8, abs(gg_xn-xd), '|back err|'), colorbar
	end

	max_percent_diff(xd, xn, mfilename)
	nbackerr(ii) = nrmse(xd,xn,ones(size(xd)));
	max_percent_diff(xd, gg_xn, mfilename)
	gg_nbackerr(ii) = nrmse(xd,gg_xn,ones(size(xd)));

end

if length(npow) > 1 % then we evaluated over data size
	figure,hold on
	plot(npow,nforwtime,'r')
	plot(npow,gg_nforwtime,'g');
	plot(npow,nbacktime,'r*')
	plot(npow,gg_nbacktime,'g*')

	xlabel 'log_2 N'
	ylabel 'Time, sec'
	legend('KB table forward','GG forward','KB table back','GG back')
	grid on

	figure,hold on
	plot(npow,nforwerr,'r')
	plot(npow,gg_nforwerr,'g')
	plot(npow,nbackerr,'r*')
	plot(npow,gg_nbackerr,'g*')

	xlabel 'log_2 N'
	ylabel 'NRMSE'
	legend('KB table forward','GG forward','KB table back','GG back')
	grid on
end
if isvar('Jval')
	if length(Jval) > 1 % then we evaluated over spreading width
		figure,hold on
		plot(Jval,nforwtime,'r')
		plot(Jval,gg_nforwtime,'g');
		plot(Jval,nbacktime,'r*')
		plot(Jval,gg_nbacktime,'g*')

		xlabel 'J'
		ylabel 'Time, sec'
		legend('KB table forward','GG forward','KB table back','GG back')
		grid on

		figure,hold on
		plot(Jval,nforwerr,'r')
		plot(Jval,gg_nforwerr,'g')
		plot(Jval,nbackerr,'r*')
		plot(Jval,gg_nbackerr,'g*')

		xlabel 'J'
		ylabel 'NRMSE'
		legend('KB table forward','GG forward','KB table back','GG back')
		grid on
	end
end
