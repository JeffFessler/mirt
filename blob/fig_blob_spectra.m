% fig_blob_spectra.m
% 1D blob-like basis functions and their spectra

% interestingly, best J=2 blob approximates raised cosine in 1d (exactly?).

t = linspace(-0,3,201)';
br = rect(t);
bs = nufft_sinc(t);
bc = (1+cos(pi*t))/2 .* (abs(t) < 1);

dim = 1;
kb_J = 3; kb_m = 2; kb_alf = 2.472;
kb_J = 4; kb_m = 2; kb_alf = 2.461;
kb_J = 4; kb_m = 2; kb_alf = 8.687;
kb_J = 4; kb_m = 2; kb_alf = 11.17;
kb_J = 2; kb_m = 2; kb_alf = 2.502;
bkf = @(t) kaiser_bessel(t, kb_J, kb_alf, kb_m);
bk = bkf(t);
nn = [-kb_J/2+1:kb_J/2-1]
bkn = bkf(nn)

% cubic spline
%cs = (x) (1-2.5*abs(x).^2+1.5*abs(x).^3) .* (abs(x) <= 1) + ...
% (2-4*abs(x)+2.5*abs(x).^2-0.5*abs(x).^3) .* (1 <= abs(x) & abs(x) <= 2);

% cubic B-spline
b3f = @(t) 1/6*(2 - abs(t)).^3 .* (1 < abs(t) & abs(t) < 2) ...
	+ (2/3 - t.^2 + 1/2 * abs(t).^3) .* (abs(t) <= 1);

b1 = (1-abs(t)) .* rect(t/2);

% canonical interpolators
[b3i H3i] = deconv_z([1/6 4/6 1/6], b3f, t); % cubic B-spline
[bki Hki] = deconv_z(bkn, bkf, t); % KB

clf, subplot(211)
plot(t, bs, ':', t, bc, '-', t, bk, '--', t, br, '-.', ...
	t, b3i, '-', ...
	t, bki, '-')
%	t, b1, '-', ...
legend('sinc', 'cos', 'kb', 'rect', 'b3', 'b1')
axisy([-0.25 1])

subplot(212)
u = linspace(-0,3,1001)';
Bs = rect(u);
Bc = nufft_sinc(2*u) + 0.5*nufft_sinc(2*(u-1/2)) + 0.5*nufft_sinc(2*(u+1/2));
Br = nufft_sinc(u);
Bk = kaiser_bessel_ft(u, kb_J, kb_alf, kb_m, dim);
%B3 = nufft_sinc(u).^4 ./ (2/3 + 1/3 * cos(2*pi*u)); % cubic b-spline blob
B3i = nufft_sinc(u).^4 ./ H3i(2*pi*u); % cubic b-spline blob
Bki = Bk ./ Hki(2*pi*u);
B1 = nufft_sinc(u).^2; % linear

db = @(x) 20*log10(abs(x) + eps);
Bs = db(Bs);
Bc = db(Bc);
Br = db(Br);
Bk = db(Bk);
Bki = db(Bki);
B3i = db(B3i);
B1 = db(B1);
plot(u, Bs, ':', u, Bc, '-', u, Bki, '--', u, Br, '-.', u, B3i, '-', u, B1, '-')
xlabel u
%axisy(1e-4, 2)
axisy(-70, 2)
if 0
	hold on, plot(kk/Np/dt, Bd, 'y+'), hold off
	axisx(minmax(u))
end


if 0 % verify analytical via dft
	N = length(t);
	Np = 4*N;
	bd = padsym2(bc, Np, 1);
	dt = t(2) - t(1);
	Bd = dt * fftshift(reale(fft(ifftshift(bd))));
	kk = [-(N-1)/2:(N-1)/2]';
	kk = [-Np/2:Np/2-1]';
	%Bd = Bd(kk+1+Np/2);
end

% outer product (separable) of raised cosine is nearly symmetric!
if 0
	%bc = 1-abs(t);
	c = [0:10]/10;
	c = (1+cos(pi*c))/2;
	bo = bc * bc';
	clf, contour(t, t, bo, c), axis square
	[xx yy] = ndgrid(t,t);
	rr = sqrt(xx.^2+yy.^2);
	br = (1+cos(pi*rr))/2 .* rect(rr/2);
	hold on
	contour(t, t, br, c, 'g-')
	hold off
	im([br; bo])
end
