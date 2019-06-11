 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
%	y = G * x	or x = G' * y

if ob.is.empty
	error empty
end

if ob.power ~= 1, error notdone, end


%
%	partial projection or backprojection (for ordered subsets)
%
if ob.is.subref
	error 'subref not done'

%	if ~ob.is.transpose
%		error 'G(:,ii) not done, only Gt(:,ii) is done!'
%	end
%
%	%
%	%	Gt(:,ii)' * x		partial forward projection
%	%
%	if ob.is.transpose_after_sub
%		if ob.is.masked
%			x = embed(x, ob.mask);
%		end
%		y = wtfmex('dsc,proj', ob.arg', single(x), ob.mask, ...
%			int32(ob.ia_start), int32(ob.ia_inc), ob.chat);
%		y = y(:,(ob.ia_start+1):ob.ia_inc:ob.na);
%
%
%	%
%	%	Gt(:,ii) * y		partial back projection
%	%
%	else
%		y = zeros(ob.nb,ob.na);
%		ia = (ob.ia_start+1):ob.ia_inc:ob.na;
%		y(:,ia) = reshape(x, ob.nb, length(ia));
%		y = wtfmex('dsc,back', ob.arg', single(y), ob.mask, ...
%			int32(ob.ia_start), int32(ob.ia_inc), ...
%			int32(ob.nthread), ob.chat);
%		if ob.is.masked
%			y = y(ob.mask);
%		end
%	end


%
%	full "projection"
%
elseif ~ob.is.transpose
	if ob.is.masked
		x = embed(x, ob.mask);
	end
	y = conv2(x, ob.psf, 'same');

%
%	full "back-projection"
%
else
	y = conv2(x, ob.psf, 'same');
	if ob.is.masked
		y = y(ob.mask);
	end
end
	y = y(:);
