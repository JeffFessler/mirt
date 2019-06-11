 function R = Rbuild(ctype, kappa, nbrs, beta, potential, delta, ...
		match_aspire, dist_power)
%function R = Rbuild(ctype, kappa, nbrs, beta, potential, delta, ...
%		match_aspire, dist_power)
%
% build a regularization "object" for penalized-likelihood reconstruction.
%
% general form of nonquadratic penalty function:
%	R(x) = \sumk \potk([Cx]_k)
% for quadratic case, \potk(t) = t^2/2 in which case
%	R(x) = x' C' C x / 2
%
% penalty gradient is C' D C x, where D = diag{\wpot_k([Cx]_k)}
%
% in
%	ctype			e.g. 'leak', see C2sparse()
%	kappa	[nx,ny]		(or mask), see C2sparse()
%	nbrs			e.g. 8, see C2sparse()
%	beta			global regularization parameter
%	potential		e.g. 'huber', see potential_func()
%	delta			potential parameter, see potential_func()
%	match_aspire		modify denominator to match aspire?
%				trick: or "-1" for no denominator
%				since R.E and R.denom only needed for SPS
%	dist[ance]_power	see C2sparse()
% out
%	R structure has the following inline "methods" 
%	R.penal(R, x)	evaluates R(x)
%	R.cgrad(R, x)	evaluates \cgrad R(x) (column gradient)
%	R.denom(R, x)	evaluates denominator for separable surrogate
%	R.dercurv(R, Cx) derivatives and curvatures for non-separable parabola
%			surrogates, returning [nk, 2] array.
%	R.diag(R)	diagonal of Hessian of R (at x=0), for preconditioners.
%
% Copyright 2002-3-14, Jeff Fessler, The University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

warning 'Rbuild_new is obsolete: use Robject'

if ~isvar('nbrs') | isempty(nbrs), nbrs = 8; end
if ~isvar('beta') | isempty(beta), beta = 1; end
if ~isvar('potential') | isempty(potential), potential = 'quad'; end
if ~isvar('delta') | isempty(delta), delta = inf; end
if ~isvar('match_aspire') | isempty(match_aspire)
	match_aspire = false;
end
if ~isvar('dist_power')
	dist_power = [];
end

if streq('ctype', 'none')
	R.C = 0;
	error 'todo: unpenalized'
end

chat = false;

R.ctype = ctype;
R.nbrs = nbrs;
R.potential = potential;
R.delta = delta;
R.isquad = streq(R.potential, 'quad');
if match_aspire == -1
	R.has_sps = false;
else
	R.has_sps = true;
end

%
% build sparse differencing matrix "C"
%
[R.C wt] = C2sparse(ctype, kappa, nbrs, chat, dist_power);
R.C = R.C(:,kappa(:)~=0);	% compact form

%
% trick: for quadratic penalty, absorb 'wt' into C
%
if R.isquad
	R.C = spdiag(sqrt(beta * wt)) * R.C;
else
	R.wt = beta * wt;	% scale wt's by beta
end
clear beta wt

%
% desired potential function
%
R.pot = potential_func(potential, delta);

%
% jth penalty separable surrogate curvature is \sumk |\ckj| \ck \wpotk
% where \ck = \sumj |\ckj|.  Define E by E_{jk} = |\ckj| \ck
%
if R.has_sps % this is only needed for separable quad. surr. algorithms!
	if R.isquad
		R.denom_vec = abs(R.C)' * sum(abs(R.C'))';	% E * 1_K
	else
		R.E = ( spdiag(sum(abs(R.C'))) * abs(R.C) )';
%		R.E = ( spdiag(sum(R.C' ~= 0)) * (R.C .^2) )';	% probably poor
	end
end


%
% make penalty curvatures consistent with ASPIRE using depierro=2 factor
% (only needed for debugging) (only slightly inefficient otherwise)
%
if R.has_sps & match_aspire
	n_per_k = sum(R.C' ~= 0);
	if max(n_per_k) == 2
		if any(n_per_k & (n_per_k ~= 2))
			warning 'using depierro = 2 version'
			n_per_k(:) = 2;
		end
	end
	R.E = ( spdiag(n_per_k) * (R.C .^2) )';
	clear n_per_k
	if R.isquad
		R.denom_vec = sum(R.E')';
	end
end


%
% inline functions
%
R.dercurv = @Rbuild_dercurv;
if R.isquad
	R.penal = inline('sum(abs(R.C * x).^2)/2', 'R', 'x');
	R.cgrad = inline('R.C'' * (R.C * x)', 'R', 'x');
	if R.has_sps
		R.denom = inline('R.denom_vec', 'R', 'x');
	end
	R.diag = inline('sum(transpose(abs(R.C.^2)))', 'R');
else
	R.penal = inline('sum(R.wt .* R.pot.potk(R.pot, R.C * x))', 'R', 'x');
	R.cgrad = 'R.C'' * (R.wt .* R.pot.dpot(R.pot, R.C * x))';
	R.cgrad = inline(R.cgrad, 'R', 'x');
	if R.has_sps
		R.denom = 'R.E * (R.wt .* R.pot.wpot(R.pot, R.C * x))';
		R.denom = inline(R.denom, 'R', 'x');
	end
	R.diag = inline('transpose(abs(R.C.^2)) * R.pot.wpot(R.pot, 0)', 'R');
end


%
% Compute both cgrad and denom (of separable surrogate) efficiently.
% Unused for now, but could be used if inlines are too inefficient
% since both R.cgrad and R.denom use C*x so there is redundancy.
% What we really need is an inline that has two output arguments.
% No, the feval with a function_handle will suffice!
%
function [cgrad, denom] = Rbuild_cgrad_denom(R, x)
if R.isquad
	cgrad = R.C' * (R.C * x);
	denom = R.denom;
else
	Cx = R.C * x;
	wx = R.wt .* R.pot.wpot(R.pot, Cx);
	cgrad = R.C' * (wx .* Cx);
	denom = R.E * wx;
end


%
% evaluate \dpoti and \wpoti
%
function [deriv, curv] = Rbuild_dercurv(R, Cx)
if R.isquad
	deriv = Cx;
	curv = ones(size(Cx));
else
	deriv = R.wt .* R.pot.wpot(R.pot, Cx) .* Cx; 
	curv = R.wt .* R.pot.wpot(R.pot, Cx);
end
