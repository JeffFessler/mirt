  function fm = de_ftab_fm(sll, mac, Ide)
%|function fm = de_ftab_fm(sll, mac, Ide)
%| calculate f_m(s1, s2, ...) (nonlinear BH function) using 
%| in
%|	sll	cell{LL} or [(Nd),LL] material density integrals
%|				{s1, s2, ..., sL} (ndgrid so same size [(Nd)])
%|	mac	[ne,LL]		mass atten coefficients
%|	Ide	[ne,MM]		spectra; same energy samples as mac
%| out
%|	fm	[(Nd),MM]	nonlinear BH function values
%|
%| Copyright 2006-05-18, Jeff Fessler, University of Michigan

if nargin == 1 && streq(sll, 'test'), de_ftab_fm_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

[ne LL] = size(mac);
MM = size(Ide,2);
if ne ~= size(Ide,1), error 'bad size Ide', end

if iscell(sll)
	sll = stackup(sll{:}); % [(Nd),LL]
	if ndims(sll) ~= LL+1, error 'bad size sll', end
end
Ns = size(sll);
if Ns(end) ~= LL, error 'bad size sll', end
Nd = Ns(1:(end-1));

Isum = sum(Ide); % [M,1]
sll = reshapee(sll, [], LL); % [*Nd, L]

fm = zeros([prod(Nd) MM]); % [*Nd, M]
tmp = 1;
for ll=1:LL
	tmp = tmp .* exp(-sll(:,ll) * mac(:,ll)'); % [*Nd, ne]
end
fm = tmp * Ide; % [*Nd, ne] * [ne,MM] -> [*Nd,MM]
fm = -log(fm ./ repmat(Isum, [prod(Nd) 1]));
if min(sll(:)) >= 0
	fm = zero_tiny_negative(fm);
end
fm = reshape(fm, [Nd MM]);


%
% de_ftab_fm_test()
%
function de_ftab_fm_test
sl{1} = linspace(0, 50, 26);
sl{2} = linspace(0, 30, 31);
xray = xray_read_spectra('ps1');
mtype = {'water', 'bone'};
mac = xray_read_atten(mtype, xray.en);
sll = ndgrid_jf('cell', sl{:});

%cpu etic
fm = de_ftab_fm(sll, mac, xray.Ide);
im clf, im(sl{1}, sl{2}, fm, 'fm'), cbar
%cpu etoc 'new'

% check old way:
if 0
cpu etic
fo = zeros([Nd MM]);
for i1=1:size(fm,1)
	for i2=1:size(fm,2)
		for mm=1:MM
			tmp = 0;
			for ll=1:LL
				tmp = tmp + mac(:,ll) * sll{ll}(i1,i2);
			end
			fo(i1,i2,mm) = -log( (Ide(:,mm)' * exp(-tmp)) / I(mm) );
		end
	end
end
fo = zero_tiny_negative(fo);
cpu etoc 'old'

max_percent_diff(fm,fo)
end
