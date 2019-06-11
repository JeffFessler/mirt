 function [mat, nexp] = filtmat(arg0, arg1, arg2, arg3, arg4)
%function mat = filtmat('1d', kernx, nx)
%			1d nonperiodic
%function mat = filtmat({'1d' 'per'}, kernx, nx)
%			1d periodic
%function mat = filtmat('2d,mul,sep', kernx, kerny, [nx ny])
%			2d multiplicatively separable nonperiodic
%function mat = filtmat('3d,mul,sep', kernx, kerny, kernz, [nx ny nz])
%			3d multiplicatively separable nonperiodic
%			kernel(x,y,z) = kernx(x) * kerny(y) * kernz(z)
%function [mat, nexp] = filtmat({'3d,mul,sep' 'expand'}, kernx, kerny, kernz, [nx ny nz])
%			3d multiplicatively separable nonperiodic, expanded
%function mat = filtmat('1d,causal', kernx, nx)
%			1d causal filter, ala convmtx
%
%	Form matrix to do 3d linear filtering - odd kernel length.
%	expanding is like what happens in 'conv' - output is longer

if ~nargin, help filtmat, error arg, end
flag_expand = any(strcmp(arg0, 'expand'));

if nargin == 3 && strcmp(arg0, '1d,causal')
	mat = filtmat1causal(arg1, arg2);
	return

elseif nargin == 3 && any(strcmp(arg0, '1d'))
	if any(strcmp(arg0, 'per'))
		mat = filtmat1per(arg1, arg2);
	else
		mat = filtmat1(arg1, arg2, flag_expand);
	end
	return

elseif nargin == 4 && any(strcmp(arg0, '2d,mul,sep'))
	mat = filtmat('3d,mul,sep', arg1, arg2, 1, [arg3 1]);
	return

elseif nargin == 5 && any(strcmp(arg0, '3d,mul,sep'))
	k.x = arg1;
	k.y = arg2;
	k.z = arg3;
	n.x = arg4(1);
	n.y = arg4(2);
	n.z = arg4(3);
else
	nargin, help filtmat, error arg
end

	if flag_expand

		[p.x, o.x] = filtmat1(k.x, n.x, flag_expand);
		[p.y, o.y] = filtmat1(k.y, n.y, flag_expand);
		[p.z, o.z] = filtmat1(k.z, n.z, 0);	% no z expansion

		p.x = kron(speye(n.y*n.z), p.x);
		p.y = kron(speye(n.z), kron(p.y, speye(o.x)));
		p.z = kron(p.z, speye(o.x*o.y));

		if 0
			n.xyz = n.x * n.y * n.z;
			t = (o.y-n.y)/2 * o.x;
			p.x = [sparse(t,n.xyz); p.x; sparse(t,n.xyz)];
			t = zeros(o.x,o.y);
			t([1:n.x]+(o.x-n.y)/2-1,:) = 1;
			tt = sparse(o.x*o.y, n.x*n.y);
			tt(find(t(:)),:) = p.y;
			p.y = tt;
			p.z = sparse(1,1);
			if n.z ~= 1, error notdone, end
		end
	else
		p.x = filtmat1(k.x, n.x, flag_expand);
		p.y = filtmat1(k.y, n.y, flag_expand);
		p.z = filtmat1(k.z, n.z, flag_expand);

		p.x = kron(speye(n.y*n.z), p.x);
		p.y = kron(speye(n.z), kron(p.y, speye(n.x)));
		p.z = kron(p.z, speye(n.x*n.y));
		o = n;
	end

%	mat = p.x + p.y + p.z;
	if flag_expand, warning('fix: expand may not work for multiplicative?'), end
	mat = p.z * p.y * p.x;
	nexp = [o.x o.y o.z];

function mat = filtmat1causal(kern, n)
	nk = length(kern);
	if nk < n, error n, end
	if n ~= nk, error 'not done', end
	t = kern(:,ones(1,n));
	t = [t; zeros(nk,nk)];
	t = t(:);
	t((end-nk+1):end) = [];
	mat = reshape(t, 2*nk-1,nk);
	mat = mat(1:nk,:);

function [mat, nexp] = filtmat1(kern, n, flag_expand)
%	1d matrix that does linear filtering - odd kernel length.
%	flag_expand means output has more entries than input (n+nk-1)

	nk = length(kern);
	ic = (nk+1)/2;		% center index
	if round(ic) ~= ic,	error('odd kernel only'), end

	if flag_expand
		n = n + nk-1;
	end
	nexp = n;

	mat = diag(ones(n,1)*kern(ic));
	for iv=1:(ic-1)
		mat = mat + diag(ones(n-iv,1) * kern(ic-iv), -iv) ...
			  + diag(ones(n-iv,1) * kern(ic+iv), iv);
	end

	if flag_expand
		mat = mat(:,ic:(end-ic+1));
	end
	mat = sparse(mat);

function mat = filtmat1per(kern, nn)
%	1d matrix that does periodic linear filtering.  odd or even kernel ok.

	nk = length(kern);
	ic = (nk+1)/2;

	mat = spdiag(ones(nn,1)*kern(ic));
	for kk=1:(ic-1)
		mat = mat ...
			+ diag(ones(nn-kk,1) * kern(ic-kk),	-kk) ...
			+ diag(ones(kk   ,1) * kern(ic-kk),	nn-kk);
	end

	for kk=1:floor(nk/2)
		mat = mat ...
			+ diag(ones(nn-kk,1) * kern(ic+kk),	kk) ...
			+ diag(ones(kk   ,1) * kern(ic+kk),	kk-nn);
	end
	mat = sparse(mat);
