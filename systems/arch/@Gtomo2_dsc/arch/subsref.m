 function out = subsref(ob, args)
%function out = subsref(ob, args)
%	evaluate call of the form ob.ref or ob(ref,:)
%	fix: needs a lot of work to handle ordered subsets!

%
%	enable handling of multiple subscripts, e.g., ob.field()
%
if length(args) > 1
	arg = args(1);
else
	arg = args;
end

%
%	ob.?
%
if arg.type == '.'
	out = getfield(struct(ob), arg.subs);

%
%	ob(?)
%
elseif arg.type == '()'
	subs = arg.subs;

	if ~isempty(ob.index1) | ~isempty(ob.index2)
		error 'subscripts after subscripts perhaps not done'
	end

	%
	%	G(:,:)
	%
	if length(subs) == 2 & subs{1} == ':' & subs{2} == ':'
		out = ob(:,[1:ob.dims(2)]);

	%
	%	Gt(:,ii) - for ordered subsets
	%
	elseif length(subs) == 2 & subs{1} == ':' ...
		& ndims(subs{2}) == 2 & ob.is_transpose
		nb = ob.nb;
		ii = subs{2};
		ia = 1 + round((ii(1,:) - 1) / nb);
		if any(ii(:) ~= col(outer_sum(1:nb,(ia-1)*nb)))
			error 'bad subref'
		end
		ob.ia_start = ia(1) - 1;
		ob.ia_inc = ia(2) - ia(1);
		if any(diff(ia) ~= ob.ia_inc)
			error 'nonuniform angles'
		end
		ob.is_subset = true;
		out = ob;

	%
	%	G(:,j)
	%
	elseif length(subs) == 2 & subs{1} == ':'
		jj = col(subs{2});
		if islogical(jj)
			if length(jj) ~= ob.dims(2)
				error 'bad column logical length'
			end
		elseif any(jj < 1 | jj > ob.dims(2))
			error 'bad column index'
		end

		%
		%	handle G(:,mask(:)) as a special case
		%
		if length(jj) == ob.nx*ob.ny & all(jj == ob.mask(:)) ...
				& ~ob.is_masked & ~ob.is_transpose
			ob.is_masked = true;
			ob.dims(2) = sum(ob.mask(:));
			out = ob;

		else
			out = zeros(ob.dims(1),length(jj));
			for nn=1:length(jj)
				x = zeros(ob.dims(2),1);
				x(jj(nn)) = 1;
				out(:,nn) = ob * x;
			end
			out = sparse(out);
		end

	%
	%	G(i,:)
	%
	elseif length(subs) == 2 & subs{2} == ':'
		ii = col(subs{1});
		if islogical(ii)
			if length(ii) ~= ob.dims(1)
				error 'bad row logical length'
			end
		elseif any(ii < 1 | ii > ob.dims(1))
			error(sprintf('bad row index %s', int2str(ii)))
		end

		out = zeros(ob.dims(2), length(ii));
		for nn=1:length(ii)
			y = zeros(ob.dims(1), 1);
			y(ii(nn)) = 1;
			out(:,nn) = ob' * y;
		end
		out = sparse(out);

	else
		arg.subs
		error todo
	end

else
	error(sprintf('type %s notdone', arg.type))
end

if length(args) > 1
	out = subsref(out, args(2:end));
end
