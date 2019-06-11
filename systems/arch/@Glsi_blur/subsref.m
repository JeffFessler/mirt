 function out = subsref(ob, arg)
%function out = subsref(ob, arg)
%	evaluate call of the form ob.ref or ob(ref,:)
%	fix: needs a lot of work to handle ordered subsets!

if arg.type == '.'
%	eval(sprintf('out = ob.%s;', arg.subs))
%	out = getfield(ob, arg.subs);
	out = getfield(struct(ob), arg.subs);

elseif arg.type == '()'
	sub = arg.subs;
	ob.is.subref = true;

	if ob.power ~= 1, error 'power: impossible', end

	%	G(:,:)
	if length(sub) == 2 & sub{1} == ':' & sub{2} == ':'
		error 'G(:,:) not done'

	%
	%	Gt(:,ii) - for ordered subsets
	%
	elseif length(sub) == 2 & sub{1} == ':' ...
		& ndims(sub{2}) == 2 & ob.is.transpose
		nb = ob.nb;
		ii = sub{2};
		ia = 1 + round((ii(1,:) - 1) / nb);
		if any(ii(:) ~= col(outer_sum(1:nb,(ia-1)*nb)))
			error 'bad subref'
		end
		ob.ia_start = ia(1) - 1;
		ob.ia_inc = ia(2) - ia(1);
		if any(diff(ia) ~= ob.ia_inc)
			error 'nonuniform angles'
		end
		ob.is.subref = 1;
		out = ob;

	%	G(:,j)
	elseif length(sub) == 2 & sub{1} == ':'
		jj = sub{2};


		%
		%	handle G(:,mask(:)) as a special case
		%
		if numel(jj(:)) == ob.nx*ob.ny & all(jj(:) == ob.mask)
			ob.is.masked = 1;
			ob.is.subref = false;
			out = ob;

		else

			if any(jj < 1 | jj > ob.nx*ob.ny)
				error(sprintf('bad j'))
			end
			out = zeros(ob.nb*ob.na,length(jj));
			for nn=1:length(jj)
				x = zeros(ob.nx, ob.ny);
				x(jj(nn)) = 1;
				out(:,nn) = ob * x(:);
			end
			out = sparse(out);
		end

	%	G(i,:)
	elseif length(sub) == 2 & sub{2} == ':'
		ii = sub{1};
		if any(ii < 1 | ii > ob.nb*ob.na)
			error(sprintf('bad i'))
		end
		out = zeros(ob.nx*ob.ny,length(ii));
		for nn=1:length(ii)
			y = zeros(ob.nb, ob.na);
			y(ii(nn)) = 1;
			out(:,nn) = ob' * y(:);
		end
		out = sparse(out);

	else
		arg.subs
		error todo
	end

else
	error(sprintf('type %s notdone', arg.type))
end
