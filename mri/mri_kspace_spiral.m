 function [kspace, omega, gradxy] = mri_kspace_spiral(varargin)
%function [kspace, omega, gradxy] = mri_kspace_spiral(varargin)
%|
%| k-space spiral trajectory based on GE 3T scanner constraints
%| options (name / value pairs)
%|	fov	field of view in cm; default 22
%|	gamp	max gradient amplitude; default 2.2 G/cm
%|	gslew	max gradient slew rate; default 180 mT/m/ms
%|	dt	time sampling interval; default 5e-6
%|	N	size of reconstructed image; default 64
%|	Nl	# of interleaves; default 1
%|	Nt	# of time points; default is 3770 for fov=22
%|
%| out
%|	kspace [Nt 2]	kspace trajectory [kx ky] in cycles/fov
%|		or [Nt 2 Nl] if Nl > 1
%|	omega [Nt 2]	"" in radians
%|	gradxy [Nt 2]	gradient waveforms in (units?)
%|
%| based on m-files from valur that he got from brad who got them from doug...

if nargin == 1 && streq(varargin{1}, 'test')
	mri_kspace_spiral_test, return
end

% defaults for spiral trajectory
arg.dt = 5e-6;	% sampling interval
arg.fov = 22;	% in cm
gamp = 2.2;	% gradient amplitude G/cm
gslew = 180;	% gradient slew mT/m/ms
arg.N = 64;	% size of the reconed image [square]
arg.Nt = [];	% # of time sampling points
arg.Nl = 1;	% # of interleves
arg.rotamount = 0; % ?
arg.rev_flag = 0; % reverse spiral?

arg = vararg_pair(arg, varargin);

% provide sensible default Nt
if isempty(arg.Nt)
	if arg.fov == 20
		arg.Nt = 4026;
	elseif arg.fov == 22
		arg.Nt = 3770;
	else
		arg.Nt = 0; % let algorithm choose
		warning 'unknown FOV: specify Nt?'
	end
end

if arg.fov > 100
	warning 'fov > 100; use cm not mm!'
end

% generate spiral k-space trajectory
[kx, ky, gx, gy] = ...
	genkspace(arg.fov, arg.N, arg.Nt, arg.Nl, gamp, gslew, ...
		arg.dt, arg.rotamount, arg.rev_flag);

if arg.Nl == 1
	kspace = [kx ky];
	gradxy = [gx gy];
else
	kspace = permute(cat(3, kx, ky), [1 3 2]);
	gradxy = permute(cat(3, gx, gy), [1 3 2]);
end
omega = 2*pi*kspace / arg.N;

if max(omega(:)) > pi, error 'bad spiral', end

end % mri_kspace_spiral()


% genkspace
%
% Generate the proper length of k-space trajectory.
% Linearly interpolate the output of genspi() to the correct length
% and take care of the rotations for the interleaves.
%
% ld is the length of the data (total, for all interleaves)
% nint is the number of interleaves
%
% Brad Sutton, University of Michigan

function [kx, ky, gx, gy] = ...
	genkspace(FOV, N, ld, nint, gamp, gslew, tsamp, rotamount, rev_flag)

nk = ld/nint;
if round(nk) ~= nk
	fail('Input should have num data pts/number of interleaves must be int')
end

flag = (nk == 0); % auto determine number of k-space points

[Gx, Gy, kxi, kyi, sx, sy, dts] = genspi(FOV, N, nint, gamp, gslew);

% resample genspi waveforms to "tsamp" over same total duration

ik = 0:(length(kxi)-1);
tk = (0:(dts/tsamp*length(kxi)))*tsamp;
kxt = interp1(ik*dts, kxi, tk)';
kyt = interp1(ik*dts, kyi, tk)';

ig = 0:(length(Gx)-1);
tg = (0:(dts/tsamp*length(Gx)))*tsamp;
gxt = interp1(ig*dts, Gx, tg)';
gyt = interp1(ig*dts, Gy, tg)';

if flag
	nk = length(kxt) - 2;
elseif nk > length(kxt)
	fail('nk = %d, length(kxt) = %d', nk, length(kxt))
end

kx = zeros(nk,nint);
ky = zeros(nk,nint);
kxo = kxt(1:nk);
kyo = kyt(1:nk);

gx = zeros(nk,nint);
gy = zeros(nk,nint);
gxo = gxt(1:nk);
gyo = gyt(1:nk);

%{
	if length(kxi)==nk;
		kxo = kxi.';
		kyo = kyi.';
	else

		kxo(1) = kxi(1);
		kyo(1) = kyi(1);

		nki = length(kxi);
		sprintf('Interpolating %d kspace pts to %d',nki,nk)

		sf = (nki-1)/(nk-1); % scaling factor
		for ii = 1:(nk-2);
			ind = ii*sf;
			kxo(ii+1) = kxi(floor(ind)+1)*(1+floor(ind)-ind)+kxi(ceil(ind)+1)*(ind-floor(ind));
			kyo(ii+1) = kyi(floor(ind)+1)*(1+floor(ind)-ind)+kyi(ceil(ind)+1)*(ind-floor(ind));
		end
		kxo(end) = kxi(end);
		kyo(end) = kyi(end);
	end
%}

% rotate matrix for proper orientation
phir = -rotamount*pi/2;
kxop = kxo*cos(phir) - kyo*sin(phir);
kyop = kyo*cos(phir) + kxo*sin(phir);
gxop = gxo*cos(phir) - gyo*sin(phir);
gyop = gyo*cos(phir) + gxo*sin(phir);

if rev_flag
	kxop = -flipud(kxop);
	kyop = -flipud(kyop);
	warn 'may need to flip gx too'
end

kx(:,1) = kxop;
ky(:,1) = kyop;
gx(:,1) = gxop;
gy(:,1) = gyop;
if nint > 1
%	printm('Performing %d rotations', nint)
	phi = 2*pi/nint;
	for ii = 1:(nint-1)
		kx(:,ii+1) = kxop*cos(ii*phi) - kyop*sin(ii*phi);
		ky(:,ii+1) = kyop*cos(ii*phi) + kxop*sin(ii*phi);
		gx(:,ii+1) = gxop*cos(ii*phi) - gyop*sin(ii*phi);
		gy(:,ii+1) = gyop*cos(ii*phi) + gxop*sin(ii*phi);
	end
end

end % genkspace()


% genspi()
% this is translation of C code from scanner, exactly what is played
% out to gradients at 4us.
%
% multi-shot spiral design
% uses Duyn's approximate slewrate limited design
% augmented with archimedian gmax limit
%
% inputs (args)
%	D = FOV, cm
%	N = matrix size
%	Tmax = longest acquisition allowed, s
%	dts = output sample spacing, s
%	gtype = trajectory type
% inputs (CVs)
%	Nl = number of interleaves
%	gamp = design grad max, G/cm
%	gslew = design slew rate, mT/m/ms
%	nramp = number of rampdown points
% out
%	Gx, Gy
%	grev
% time is in sec
%
%	rev 0 12/26/98	original
%	rev 1 4/15/99	little better calc of ts
%
% borrowed from Doug Noll, Univ. of Michigan
% modified to take more input cv's
function [Gx, Gy, kx, ky, sx, sy, gts] = genspi(D, N, Nl, gamp, gslew)

%% Predefined variables

GRESMAX = 21000;
if ~exist('gamp','var')
	gamp=2.2; %3.50; % 2.2 for both 1.5 T and 3 T data
end
if ~exist('gslew','var')
	gslew=180 % 200 % 180 for 3T data and 120 (150) for 1.5 T data
end
nramp = 0;
% nramp = 100;
MAX_PG_WAMP = 32766;

gts = 4e-6; % [sec]

Tmax = GRESMAX*gts;

dts = gts;
opfov = D;

%%

gamma = 2*pi*4.257e3;
gambar = gamma/(2*pi);

gx = zeros(1,2*GRESMAX);
gy = zeros(1,2*GRESMAX);

q = 5;
S0 = gslew*100;
dt = dts*.5;

% slew-rate limited approximation

Ts = .666667 / Nl * sqrt(((pi*N)^3) / (gamma*D*S0));
if (Ts > Tmax) fail('slew limited readout too long'), end

a2 = N*pi / (Nl*(Ts^(.666667)));
a1 = 1.5*S0/a2;
beta = S0*gamma*D/Nl;
Gmax = a1*(Ts^.333333);
gmax = 0;

% printm('Ts=%g', Ts)

t = 0:dt:Ts;
x = t .^ 1.333333;
theta = (t.^2) .* (.5 * beta ./ (q + .5*beta./a2 .* x));
y = q + .5 * beta / a2 * x;
dthdt = t .* (beta .* (q + .166667 * beta / a2 * x) ./ (y.*y));
c = cos(theta);
s = sin(theta);
gx = (Nl/(D*gamma)) .* dthdt .* (c - theta .* s);
gy = (Nl/(D*gamma)) .* dthdt .* (s + theta .* c);
gabs = abs(gx + 1i * gy);
% cut short if over peak
gmax = abs(gamp ./ (theta+eps) + 1i * gamp);
l1 = length(t) - sum(gabs > gmax);
ts = t(l1);
thetas = theta(l1);

% gmax limited approximation

l3 = 0;
T = ts;
if Gmax > gamp
	T = ((pi*N/Nl)*(pi*N/Nl) - thetas*thetas) / (2*gamma*gamp*D/Nl) + ts;
	if T > Tmax
		fail('gmax limited readout too long')
	end
	t = (ts+dt):dt:T;
	theta = sqrt(thetas*thetas + (2*gamma*gamp*D) .* (t-ts) / Nl);
	c = cos(theta);
	s = sin(theta);
	ind2 = l1 + (1:length(t));
	gx(ind2) = gamp .* (c ./ theta - s);
	gy(ind2) = gamp .* (s ./ theta + c);
	l3 = length(t);
end

l2 = l1 + l3;

Gx = gx(1:2:l2); % or gx(1:2:l2)*MAX_PG_WAMP/gamp
Gy = gy(1:2:l2); % or gy(1:2:l2)*MAX_PG_WAMP/gamp
g = Gx + 1i * Gy; % slew rate vector
s = diff(g) / (gts*1000); % grad vector
Kx = cumsum([0 Gx]) * gts * opfov * gambar;
Ky = cumsum([0 Gy]) * gts * opfov * gambar;
k = Kx + 1i * Ky; % kspace vector
t = 0:gts:T; % time vector
% matrix = max(abs(k)) * 2;
maxg = max(abs(g));
maxs = max(abs(s));
maxt = max(t) * 1000;

kx = real(k);
ky = imag(k);
sx = real(s);
sy = imag(s);

end % genspi()


function mri_kspace_spiral_test
	N = 64;
	[k0, o0, g0] = mri_kspace_spiral;
	k20 = mri_kspace_spiral('fov', 20, 'Nt', 0);
	k22 = mri_kspace_spiral('fov', 22, 'Nt', 0);
	if im
		im plc 2 2

		im subplot 1
		plot(k0(:,1), k0(:,2), 'b.')
		hold on
		plot(k20(:,1), k20(:,2), 'r.', k22(:,1), k22(:,2), 'g.')
		hold off
		axis([-1 1 -1 1]*N/2), axis square
		xtick([-1 0 1] * N/2); ytick([-1 0 1] * N/2)

		im subplot 2
		plot(g0, '.-')
		xlim([1, size(g0,1)])
		xtick([1, size(g0,1)])
		ylim([-1,1]*2.2)

		im subplot 3
		[kl,~,gl] = mri_kspace_spiral('Nl', 5); % interleaves
		plot(kl(:,1,1), kl(:,2,1), 'b.-')
		hold on
		for ii=2:5
			plot(kl(:,1,ii), kl(:,2,ii), '.-')
		end
		hold off
		axis([-1 1 -1 1]*N/2), axis square
		xtick([-1 0 1] * N/2); ytick([-1 0 1] * N/2)

		im subplot 4
		plot(squeeze(gl(:,1,:)), '.-')
	end
	
end % mri_kspace_spiral_test()
