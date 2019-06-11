 function out = subsref(ob, args)
%function out = subsref(ob, args)
% handle subscript references like ob.ref or ob(ref,:)
% Copyright 2002-2-20, Jeff Fessler, University of Michigan

out = fatrix2_subsref1(ob, args(1)); % handle first subscript

if numel(args) > 1 % handle multiple subscripts, e.g., ob.field()
	out = subsref(out, args(2:end));
end


% fatrix2_subsref1()
function out = fatrix2_subsref1(ob, arg)

switch(arg.type)

case '.' % ob.?
	out = fatrix2_subsref_dot(struct(ob), arg);
	return


case '()' % ob(?)
	subs = arg.subs;

	% A(:,:) or A(:,j)
	if length(subs) == 2 && streq(subs{1}, ':')
		out = fatrix2_subsref_colon(ob, subs{2});

	% A(i,:)
	elseif length(subs) == 2 && streq(subs{2}, ':')
		out = fatrix2_subsref_colon(ob', subs{1})'; % trick: (A'(:,i))'

	% A(i,j)
	elseif length(subs) == 2
		out = fatrix2_subsref_colon(ob, subs{2}); % A(:,j) first
		if numel(subs{2}) > 1
			out = full(out); % fatrix2 to array
		end
		out = out(subs{1},:); % A(i,j) using normal array indexing

	% A(:)
	elseif length(subs) == 1 && streq(subs{1}, ':')
		out = col(fatrix2_subsref_colon(ob, ':')); % col(A(:,:))

	% A(k)
	elseif length(subs) == 1 && isnumeric(subs{1})
		[ii, jj] = ind2sub(size(ob), subs{1});
		out = zeros(numel(jj),1);
		for ll=1:numel(jj)
			tmp = fatrix2_subsref_colon(ob, jj(ll)); % A(:,j)
			out(ll) = tmp(ii(ll)); % A(i,j)
		end

	else
		printm('A subsref of type () called with these args:')
		pr subs
		error('That subsref type is unknown.  Do you mean {1}?')
	end


case '{}' % ob{?} (for ordered-subsets / block methods)
	if length(arg.subs) ~= 1, fail('"{?}" usage'), end

	out = fatrix2_subsref_brace(ob, arg.subs{1});

otherwise
	fail('subscript type %s unknown', arg.type)
end
