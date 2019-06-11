 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
%	y = G * x	or x = G' * y
%	Copyright 2002-2-20	Jeff Fessler	The University of Michigan

%
%	scalar * G
%
if isa(ob, 'double') & length(x) == 1 & isa(x, 'Gtomo2_dsc')
	y = x;
	y.scale = ob;
	return
end

if ob.apower ~= 1, error notdone, end


%
%	partial projection or backprojection (for ordered subsets)
%
if ob.is_subset
	%
	%	Gt(:,ii)' * x		partial forward projection
	%
	if ~ob.is_transpose
		if ob.is_masked
			x = embed(x, ob.mask);
		end
		y = wtfmex('dsc,proj', ob.arg', single(x), ob.mask, ...
			int32(ob.ia_start), int32(ob.ia_inc), ...
			int32(ob.nthread), ob.chat);
		y = y(:,(ob.ia_start+1):ob.ia_inc:ob.na);


	%
	%	Gt(:,ii) * y		partial back projection
	%
	else
		y = zeros(ob.nb,ob.na);
		ia = (ob.ia_start+1):ob.ia_inc:ob.na;
		y(:,ia) = reshape(x, ob.nb, length(ia));
		y = wtfmex('dsc,back', ob.arg', single(y), ob.mask, ...
			int32(ob.ia_start), int32(ob.ia_inc), ...
			int32(ob.nthread), ob.chat);
		if ob.is_masked
			y = y(ob.mask);
		end
	end


%
%	full projection
%
elseif ~ob.is_transpose
	if ob.is_masked
		x = embed(x, ob.mask);
	end
	y = wtfmex('dsc,proj', ob.arg', single(x), ob.mask, ...
		int32(0), int32(1), int32(ob.nthread), ob.chat);


%
%	full back-projection
%
else
	y = wtfmex('dsc,back', ob.arg', single(x), ob.mask, ...
		int32(0), int32(1), int32(ob.nthread), ob.chat);
	if ob.is_masked
		y = y(ob.mask);
	end
end

y = ob.scale * double(y(:));
