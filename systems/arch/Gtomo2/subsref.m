 function out = subsref(ob, args)
%function out = subsref(ob, args)
% handle subscript references like ob.ref or ob(ref,:)
% Copyright 2002-2-20, Jeff Fessler, The University of Michigan

%
% handle multiple subscripts, e.g., ob.field()
%
if length(args) > 1
	arg = args(1);
else
	arg = args;
end

%
% ob.?
%
if arg.type == '.'
	out = getfield(struct(ob), arg.subs);

%
% ob(?)
%
elseif arg.type == '()'
	subs = arg.subs;

	if ~isempty(ob.index1) | ~isempty(ob.index2)
		error 'subscripts after subscripts perhaps not done'
	end

	%
	% G(:,:) or G(:,j)
	%
	if length(subs) == 2 & ischar(subs{1}) & streq(subs{1}, ':')
		if ischar(subs{2}) & streq(subs{2}, ':')
			jj = [1:ob.dims(2)]';
		elseif isnumeric(subs{2}) | islogical(subs{2})
			jj = col(subs{2});
		else
			error 'bad G(:,?)'
		end

		if islogical(jj)
			if length(jj) ~= ob.dims(2)
				error 'bad column logical length'
			end
		else
			bad = jj < 1 | jj > ob.dims(2);
			if any(bad)
				printf('bad column indeces:')
				minmax(jj(bad))
				'subsref problem'
				keyboard
			end
		end

		%
		% handle G(:,mask(:)) as a special case
		%
		if ~ob.is_masked & ~ob.is_transpose & islogical(jj) ...
			& all(jj == ob.mask(:))
			% trick: if jj is logical, then it better have
			% have the same length as ob.dims(2)
			% length(jj) == ob.dims(2)
			ob.is_masked = true;
			ob.dims(2) = sum(ob.mask(:));
			out = ob;

		else
			try
				out = zeros(ob.dims(1),length(jj));
			catch
				printm 'masking problem?'
				keyboard
			end
			for nn=1:length(jj)
				x = zeros(ob.dims(2),1);
				x(jj(nn)) = 1;
				out(:,nn) = ob * x;
			end
			out = sparse(out);
		end

	%
	% G(i,:)
	%
	elseif length(subs) == 2 & subs{2} == ':'
		ii = col(subs{1});
		if islogical(ii)
			if length(ii) ~= ob.dims(1)
				error 'bad row logical length'
			end
		elseif any(ii < 1 | ii > ob.dims(1))
			error(sprintf('bad row index %s', num2str(ii)))
		end

		out = zeros(ob.dims(2), length(ii));
		for nn=1:length(ii)
			y = zeros(ob.dims(1), 1);
			y(ii(nn)) = 1;
			out(:,nn) = ob' * y;
		end
		out = sparse(out);

	else
		disp(arg.subs)
		error todo
	end

else
	error(sprintf('type %s notdone', arg.type))
end

if length(args) > 1
	out = subsref(out, args(2:end));
end
