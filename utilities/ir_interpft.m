  function yy = ir_interpft(xx, K, dim)
%|function yy = ir_interpft(xx, K, dim)
%|
%| Interpolate dimension "dim" (default 1) of input vector or array xx
%| using an FFT-based approach so that it has K > size(xx, dim) samples.
%| Unlike matlab or octave interpft, if the input is real, this version 
%| always produces a real output.
%|
%| 2013-02-17, Jeff Fessler, University of Michigan

if nargin == 1 && streq(xx, 'test'), ir_interpft_test, return, end 
if nargin < 2, ir_usage, end
if nargin < 3
	dim = 1;
end

if isrow(xx)
	yy = ir_interpft(xx.', K).';
return
end

isize = size(xx);
if dim ~= 1
	order = 1:numel(isize);
	order(1) = dim;
	order(dim) = 1;
	xx = permute(xx, order);
end

N = size(xx,1);
xx = reshapee(xx, N, []); % [N *]

Xk = fft(xx);

osize = [K isize(2:end)];

if rem(N,2) && ~rem(K,2) % odd N, even K
	Yk = zeros(osize, class(Xk));
	kh = [0:floor(N/2)]; % right half
	Yk(1+kh,:) = Xk(1+kh,:);
	kh = [ceil(N/2):N-1]; % left half
	Yk(1+kh+K-N,:) = Xk(1+kh,:);
	yy = K/N * ifft(Yk); % [K *]

elseif rem(K,2)
	fail 'only even K done'

else % N and K even
	Yk = zeros(osize, class(Xk));
	kh = [0:N/2-1]; % right half
	Yk(1+kh,:) = Xk(1+kh,:);
	kh = [N/2+1:N-1]; % left half
	Yk(1+kh+K-N,:) = Xk(1+kh,:);
	kh = N/2; % middle point (at +/- pi)
	Xmid = Xk(1+kh,:); % these three lines are missing from interpft!
	Yk(1+kh,:) = Xmid/2; % trick
	Yk(1+kh+K-N,:) = Xmid/2; % "" trick: no conjugate!
	yy = K/N * ifft(Yk); % [K *]
end

yy = reshape(yy, osize);
if dim ~= 1
	yy = permute(yy, order);
end

if isreal(xx)
	yy = reale(yy);
end



% ir_interpft_test
function ir_interpft_test

N = 4; K = 10; % stress test using non-integer magnification
%N = 5; K = 12;
%N = 14; K = 22;
ff = @(x) (1i) * ...
	(5 - 6 * cos(2*pi * (1/N) * x - pi/7) + 7 * cos(2*pi * (2/N) * x) ...
	+ 0i * cos(2*pi * (2/N) * x));
mag = K/N;
n0 = [0:N-1]';
n1 = [0:K-1]';
x0 = n0;
x1 = n1 / mag;
f0 = ff(x0); % coarse
f1 = ff(x1); % fine

gn = interpft(f0, K);
hn = ir_interpft(f0, K);
if isreal(f0)
%	gn = reale(gn); % fails due to suboptimal interpft()
	hn = reale(hn);
end

equivs(hn, f1, 'fail', 0)

xx = linspace(0, N, 10*N+1);
fx = ff(xx);

Fk = fft(f0);
Gk = fft(gn) * N/K;
Hk = fft(hn) * N/K;
%keyboard

if 0
	Fk = reale(Fk);
	Gk = reale(Gk);
	Hk = reale(Hk);
end

if im
	im plc 2 1
	im subplot 1
	plot(xx, real(fx), ':', x0, real(f0), 'o', ...
        	x1, real(gn), 'x', x1, real(hn), '+')
	xtick([0:N])
	legend('Re f(x)', 'Re f(n)', 'Re g[n]', 'Re h[n]', ...
		'location', 'northeast')

	im subplot 2
	plot(xx, imag(fx), ':', x0, imag(f0), 'o', ...
        	x1, imag(gn), 'x', x1, imag(hn), '+')
	xtick([0:N])
	xlabelf('x')
	legend('Im f(x)', 'Im f(n)', 'Im g[n]', 'Im h[n]', ...
		'location', 'southeast')
	title 'The imaginary part of interpft() output is undesirable'
end
