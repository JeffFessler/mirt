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
	sub = arg.subs;

	if ~isempty(ob.index1) | ~isempty(ob.index2)
		error 'subscripts after subscripts perhaps not done'
	end

	%	fix: needs to handle ordered subsets!

	%
	%	G(:,:)
	%
	if length(sub) == 2 & sub{1} == ':' & sub{2} == ':'
		out = ob(:,[1:ob.dims(2)]);

	%
	%	G(:,j)
	%
	elseif length(sub) == 2 & sub{1} == ':'
		jj = col(sub{2});
		if islogical(jj)
			if length(jj) ~= ob.dims(2)
				error 'bad column logical length'
			end
		elseif any(jj < 1 | jj > ob.dims(2))
			error(sprintf('bad column index %s', int2str(jj)))
		end

		%
		%	handle G(:,mask(:)) as a special case
		%
		if length(jj) == ob.nx*ob.ny*ob.nz & all(jj == ob.mask(:)) ...
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
	elseif length(sub) == 2 & sub{2} == ':'
		ii = col(sub{1});
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
