 function out = subsref(ob, arg)
%function out = subsref(ob, arg)
%	evaluate call of the form ob.ref or ob(ref,:)

%
%	ob.?
%
if arg.type == '.'
	out = getfield(struct(ob), arg.subs);


%
%	ob(?)
%
elseif arg.type == '()'

	if ~isempty(ob.index1) | ~isempty(ob.index2)
		error 'subscripts after subscripts perhaps not done'
	end

	sub = arg.subs;

	%
	%	G(:,:)
	%
	if length(sub) == 2 & sub{1} == ':' & sub{2} == ':'
		out = ob.G;

	%
	%	G(:,j)
	%
	elseif length(sub) == 2 & sub{1} == ':'
		jj = col(sub{2});
		if islogical(jj)
			if length(jj) ~= ob.dims(2)
				error 'bad column logical length'
			end
		else
			bad = jj < 1 | jj > ob.dims(2);
			if any(bad)
				printf('bad column indeces:')
				jj(bad)'
				error 'subsref problem'
			end
		end

		%
		%	handle G(:,mask(:)) as a special case
		%
		if length(jj) == ob.nx*ob.ny & all(jj == ob.mask(:)) ...
				& ~ob.is_masked & ~ob.is_transpose
			ob.is_masked = true;
			ob.dims(2) = sum(ob.mask(:));
			ob.G = ob.G(:,ob.mask(:));
			out = ob;

		else
%			G = ob.G;
			out = subsref(ob.G, arg);
		end

	%
	%	G(i,:)
	%
	elseif length(sub) == 2 & sub{2} == ':'
		out = subsref(ob.G, arg);

	else
		arg.subs
		error todo
	end

else
	error(sprintf('type %s notdone', arg.type))
end
