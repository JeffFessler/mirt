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

switch(arg.type)

%
% ob.?
%
case '.'
	out = struct(ob);
	try
		out = out.(arg.subs);
	catch
		out
		error(['No field ' arg.subs])
	end

%
% ob(?)
%
case '()'
	subs = arg.subs;

	%
	% G(:,:) or G(:,j)
	%
	if length(subs) == 2 && ischar(subs{1}) && streq(subs{1}, ':')
		if ischar(subs{2}) && streq(subs{2}, ':')
			jj = [1:ob.dim(2)]';
		elseif isnumeric(subs{2}) || islogical(subs{2})
			jj = col(subs{2});
		else
			error 'bad G(:,?)'
		end

		if islogical(jj)
			if length(jj) ~= ob.dim(2)
				error 'bad column logical length'
			end
		else
			bad = jj < 1 | jj > ob.dim(2);
			if any(bad)
				printm('bad column indeces:')
				minmax(jj(bad))
				printm('subsref problem')
				keyboard
			end
		end

		%
		% do G(:,?) by matrix multiplication
		%
		for nn=1:length(jj)
			x = zeros(ob.dim(2),1);
			x(jj(nn)) = 1;
			tmp = ob * x;
			if (nn == 1)
				out = zeros(ob.dim(1),length(jj), class(tmp));
			end
			out(:,nn) = tmp;
		end

	%
	% G(i,:)
	%
	elseif length(subs) == 2 && streq(subs{2}, ':')
		ii = col(subs{1});
		if islogical(ii)
			if length(ii) ~= ob.dim(1)
				error 'bad row logical length'
			end
		elseif any(ii < 1 | ii > ob.dim(1))
			fail('bad row index %s', num2str(ii))
		end

		for nn=1:length(ii)
			y = zeros(ob.dim(1), 1);
			y(ii(nn)) = 1;
			tmp = (ob' * y)';
			if nn == 1
				out = zeros(length(ii), ob.dim(2), class(tmp));
			end
			out(nn,:) = tmp;
		end

	else
		printm('Fatrix subsref of type () called with these args:')
		disp(arg.subs)
		fail('That Fatrix subsref type is not done.  Do you mean {1}?')
	end


%
% ob{?} (for ordered-subsets / block methods)
%
case '{}'
	if length(arg.subs) ~= 1, fail('"{?}" usage'), end

	out = ob;
	if isempty(ob.nblock), error 'not a block object', end

	% note: user can use subset_starts to select blocks in other orders
	% in which case this "iblock" is really "istart"
	out.iblock = arg.subs{1};
	if out.iblock < 1 || out.iblock > ob.nblock
		error 'bad block index'
	end

	odim = ob.arg.odim;
	na = odim(end);
	ia = out.iblock:ob.nblock:na;
	out.dim = [prod(odim(1:end-1))*length(ia) ob.dim(2)]; % 2009-6-5

otherwise
	fail('type %s notdone', arg.type)
end

if length(args) > 1
	out = subsref(out, args(2:end));
end
