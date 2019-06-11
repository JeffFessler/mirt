 function [C, wjk] = C2sparse(type, kappa, nbrs, chat, dist_power)
%function [C, wjk] = C2sparse(type, kappa, nbrs, chat, dist_power)
%
% create 2D penalty matrix C with 2nd-order pixel neighborhood
% that performs 1st-order finite differences.
% R = beta * C'*D(wjk)*C is the penalty Hessian (in quadratic case)
% in:
%	type	string		choices:
%				'leak' - allow penalty across mask edge
%					(consistent with aspire)
%				'tight' - no penalty across mask edges
%	kappa	[n1 n2]		if logical, then just support (aka mask)
%					if double, then as in fessler:96:srp
%	nbrs	scalar		4 or 8 neighbors
%	dist_power scalar	scale wjk by 1/distance^dist_power
%		default is 1, so diagonal neighbors wjk = 1 / sqrt(2)^1
%		but i think dist_power=2 may be better for wjk since
%		i think xj-xk should be scale by the distance
% out:
%	C [n1*n2*nbrs/2,n1*n2]	most rows have a 1 and a -1
%	wjk [n1*n2*nbrs/2,1]	penalty weights:
%				kappa_j * kappa_k * [1 or 1/sqrt(2)^?]
%
% Caller may want to do:	C = C(:,mask(:));
%
% Copyright 2002-1-30, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(type, 'test'), C2sparse_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end


if ~isvar('chat'), chat = 0; end
if ~isvar('dist_power') || isempty(dist_power), dist_power = 1; end

if streq(type, 'leak')
	leak = true;
elseif streq(type, 'tight')
	leak = false;
else
	error 'bad type'
end

if nbrs == 8
	offset1 = [-1 +0 -1 +1];
	offset2 = [+0 -1 -1 -1];
elseif nbrs == 4
	offset1 = [-1 +0];
	offset2 = [+0 -1];
else
	error 'bad nbrs'
end

[n1 n2 n3] = size(kappa);
if n3 > 1, error '3d not done', end


%
% build matrix, stacking up a block for each neighbor
%
j1 = 1:n1;
j2 = 1:n2;
[j1 j2] = ndgrid(j1, j2);
jj = sub2ind([n1 n2], j1, j2);

C = [];
wjk = zeros(n1*n2, length(offset1));

%
% loop over all neighboring pixels
%
for io=1:length(offset1)
	o1 = offset1(io);
	o2 = offset2(io);

	k1 = j1 + o1;	% k: index of neighbor
	k2 = j2 + o2;

	%
	% only k's within image rectangle can have '-1' entries
	%
	krect = (k1 >= 1) & (k1 <= n1) & (k2 >= 1) & (k2 <= n2);

	if any(krect(:))
		if 0 && chat
			printf('doing offset %d of %d', io, length(offset1))
		end
		kk = ones(size(k1));	% dummy '1' value for non-rect cases
		kk(krect) = sub2ind([n1 n2], k1(krect), k2(krect));

		mask = kappa(:) > 0;
		kmask = zeros(size(mask));
		kmask(krect(:)) = mask(kk(krect));

		if leak		% put 1 or -1 if either j or k is within mask
			jvalue = mask;
			kvalue = kmask;
		else		% both j and k must be within mask
			jvalue = mask & kmask;
			kvalue = jvalue;
		end
		Ctmp	= sparse(1:n1*n2, jj(:), jvalue, n1*n2, n1*n2) ...
			- sparse(1:n1*n2, kk(:), kvalue, n1*n2, n1*n2);

		wjk(:,io) = (jvalue | kvalue) / sqrt(o1^2 + o2^2)^dist_power;

		%
		%	include kappa effect, but handle "leak" kappas carefully
		%	cf: wjk(:,io) *= kappa(jj) .* kappa(kk)
		%
		kappa2 = zeros(size(mask));
		jk = mask & kmask;
		kappa2(jk) = kappa(jk) .* kappa(kk(jk));
		jk = mask & ~kmask;
		kappa2(jk) = kappa(jk).^2;
		jk = ~mask & kmask;
		kappa2(jk) = kappa(kk(jk)).^2;
		wjk(:,io) = wjk(:,io) .* kappa2(:);

		C = [C; Ctmp];
	else
		if chat
			printf('skipping offset %d of %d', io, length(offset1))
		end
		C = [C; sparse(n1*n2,n1*n2)];
	end
end

wjk = wjk(:);

%
% optional display
%
if chat
	im clf
	display_it(type, C, wjk, kappa, n1, n2, nbrs, [131 132 133])
end

%
% internal display
%
function display_it(type, C, wjk, kappa, n1, n2, nbrs, where)
if ~im, return, end

subplot(where(1)), spy(C), title(['C ' type])
if ncol(C) < 500
	im(where(2), C'*C, 'C''C')
else
	im(where(2), embed(diag(C'*C), kappa>0), 'Diag(C''C)'), cbar
end
tmp = reshape(wjk, [n1 n2 nbrs/2]);
im(where(3), tmp, 'wjk', [0 max(wjk(:))]), cbar


%
% self test
%
function C2sparse_test
n1 = 5; n2 = 4;
kappa = ones(n1,n2);
nbrs = 8;
[C wjk] = C2sparse('leak', kappa, nbrs, 0);
im clf
display_it('leak', C, wjk, kappa, n1, n2, nbrs, [141 242 246])
[C wjk] = C2sparse('tight', kappa, nbrs, 0);
display_it('tight', C, wjk, kappa, n1, n2, nbrs, [143 244 248])
