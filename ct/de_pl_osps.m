 function xs = de_pl_osps(x, Gb, ymi, Im, rmi, R, ftab, Gtab, masseff, niter, pixmax, gi, denom)
%function xs = de_pl_osps(x, Gb, ymi, Im, rmi, R, ftab, Gtab, masseff, niter, pixmax, gi, denom)
% The Dual-Energy PL-OS-SPS algorithm for transmission Poisson problem
% (ordered subsets, separable paraboloidal surrogates)
% uses fast precomputed curvatures described in paper fessler::
% in:
%	x		[np,2] initial guess
%	Gb		Gblock object
%	ymi		[nb,na,2] raw polyenergetic measurements
%	Im		[2] source intensities
%	rmi		[nb,na,2] (or scalar) background/scatter
%	R		penalty object (see Robject.m)
%	ftab		F table
%	Gtab		G table
%	masseff		[2,2] effective mass atten coef's
%			2 energies by 2 materials
%	niter		# iterations (including initial guess)
%			fix: allow relaxation parameters here!
%	pixmax		upper constraint for pixel values
%		can be scalar (e.g. 'inf') or an array the size of x
%	gi		sum(Gt)'
%		denom		[np,2] (if precomputed denominator)
% out:
%	x [np,2,niter]	updated image vectors each iteration
%
% Copyright 2001-04-28, Jeff Fessler, The University of Michigan

if nargin < 6, help(mfilename), error(mfilename), end

if ~isvar('Im') | isempty(Im)
	Im = [1 1];
end
if ~isvar('rmi') | isempty(rmi)
	rmi = zeros(size(ymi));
end

trl_check(ymi(:,:,1), [], rmi(:,:,1));
trl_check(ymi(:,:,2), [], rmi(:,:,2));

if ~isvar('R') | isempty(R), error 'R required', end

nblock = block_ob(Gb, 'n');
if ~isvar('niter')	| isempty(niter),	niter = 1;	end
if ~isvar('pixmax')	| isempty(pixmax),	pixmax = inf;	end
if ~isvar('chat')	| isempty(chat),	chat = 0;	end

starts = subset_start(nblock);
[nb, na, nm] = size(ymi);

%
% precompute denominator
%
if ~isvar('denom'), error 'use de_pl_denom', end

%
% loop over iterations
%
xs = zeros([size(x) niter]);	% [np,2,niter]
x = max(x,0);
xs(:,:,1) = x;
for it=2:niter
	printf('DE-PL-OSPS iteration %d', it)

	%
	% loop over subsets
	%
	for iset=1:nblock
		iblock = starts(iset);
		ia = iblock:nblock:na;

		% s=G*x "mass integrals"
		s1 = Gb{iblock} * x(:,1);
		s2 = Gb{iblock} * x(:,2);

		% predicted meas. means
		s1 = reshape(s1, nb, length(ia));
		s2 = reshape(s2, nb, length(ia));
		fh = ftab.feval(ftab, s1, s2);
		yb1 = Im(1) * exp(-fh(:,:,1)) + rmi(:,ia,1);
		yb2 = Im(2) * exp(-fh(:,:,2)) + rmi(:,ia,2);

		dh1 = ymi(:,ia,1) ./ yb1 - 1;	% m=1
		dh2 = ymi(:,ia,2) ./ yb2 - 1;	% m=2

	if 0
		G11 = interp2(Gtab.s1, Gtab.s2, Gtab.Gm1(:,:,1)', s1, s2);
		G21 = interp2(Gtab.s1, Gtab.s2, Gtab.Gm1(:,:,2)', s1, s2);
		G12 = interp2(Gtab.s1, Gtab.s2, Gtab.Gm2(:,:,1)', s1, s2);
		G22 = interp2(Gtab.s1, Gtab.s2, Gtab.Gm2(:,:,2)', s1, s2);

		% undo nonlinearity associated with table of Gml
		G11 = -exp(-G11);
		G12 = -exp(-G12);
		G21 = -exp(-G21);
		G22 = -exp(-G22);

	else
		% Gml using polynomial model!
		% note: 'coef' is indexed by m, but d_l
if it*iset == 2, warning 'not tested', end
		d11 = reshape(ftab.basis_d1(s1(:), s2(:)) * ftab.coef(:,1), ...
				size(s1));
		d21 = reshape(ftab.basis_d1(s1(:), s2(:)) * ftab.coef(:,2), ...
				size(s1));
		d12 = reshape(ftab.basis_d2(s1(:), s2(:)) * ftab.coef(:,1), ...
				size(s2));
		d22 = reshape(ftab.basis_d2(s1(:), s2(:)) * ftab.coef(:,2), ...
				size(s2));
		G11 = -(yb1-rmi(:,ia,1)) .* d11;
		G21 = -(yb2-rmi(:,ia,2)) .* d21;
		G12 = -(yb1-rmi(:,ia,1)) .* d12;
		G22 = -(yb2-rmi(:,ia,2)) .* d22;
	end

		% g_l is sum over m
		g1 = G11 .* dh1 + G21 .* dh2;
		g2 = G12 .* dh1 + G22 .* dh2;

		grad = Gb{iblock}' * [g1(:) g2(:)];	% [np,2]

%%		for isub=1:1
			for ll=1:2
				pgrad(:,ll) = R{ll}.cgrad(R{ll}, x(:,ll));
				Rdenom(:,ll) = R{ll}.denom(R{ll}, x(:,ll));
			end

			num = nblock * grad - pgrad;
			den = denom + Rdenom;	% semi-constant denom

			x = x + num ./ den;		% the update!
%%			x = x + ((x-xold) .* Lden + num) ./ den; % the update!
			x = max(x,0);		% enforce nonnegativity
			x = min(x,pixmax);	% enforce upper bound constraint
%%		end
	end

	if chat, printf('Range %g %g', min(x), max(x)), end
	xs(:,:,it) = x;
end
