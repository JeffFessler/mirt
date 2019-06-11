 function y = mtimes(ob, x)
%function y = mtimes(ob, x)	for Gnufft object
% y = G * x	or x = G' * y

if ob.apower ~= 1, error 'power not done', end


%
% partial projection or backprojection (for ordered subsets)
%
if ob.is_subref
	error 'subref not done'

%
% full "projection"
%
elseif ~ob.is_transpose
	if ob.is_masked
		x = embed(x, ob.mask);
	end

	%
	% full dsft, which doesn't have shift implemented anyway
	%
%	if ob.is.dsft
%		x = reshape(x, ob.Nd);
%		y = dtft_mex('forward', ob.omega', x, int32(ob.nthread));

	%
	% nufft
	%
%	else
		x = reshape(x, [ob.Nd numel(x)/prod(ob.Nd)]);
		y = nufft(x, ob.st);
%	end


%
% full "back-projection"
%
else

%	if ob.is.dsft
%		if isreal(y)
%			y = complexify(y);
%			warning 'faking complex'
%		end
%		y = dtft_mex('adjoint', ob.omega', y, ...
%			int32([ob.nx ob.ny]'), int32(ob.nthread));
%		y = ob.pixel_size * y;

%	else
		y = nufft_adj(x, ob.st);
%	end

	if ob.is_masked
		y = y(ob.mask);
	end

	y = reshape(y, [ob.dims(1) numel(y) / ob.dims(1)]); % [*Nd,L]
end
