 function out = fatrix2_subsref_colon(ob, sub2, varargin)
%function out = fatrix2_subsref_colon(ob, sub2, varargin)
%
% handle subscript references like ob(:,sub2)
% This will called from ../subsref with (ob, subs{2})
%
% Copyright 2010-12-02, Jeff Fessler, University of Michigan

if streq(sub2, ':') % ob(:,:)
	% do by col if thin or if user-defined forw() with default back()
	user_forw =	~isequal(ob.handle_forw, @fatrix2_def_forw) && ...
			~isequal(ob.handle_forw, @fatrix2_def_back);
	do_by_col = ((size(ob,1) >= size(ob,2)) && user_forw) || ...
		(isequal(ob.handle_back, @fatrix2_def_back));

	if nargin > 2 % for full(ob, 'col') in test_adjoint.m
		if numel(varargin) == 1 && ischar(varargin{1})
			switch varargin{1}
			case 'col'
				do_by_col = true;
			case 'row'
				do_by_col = false;
			otherwise
				fail('unknown option "%s"', varargin{1})
			end
		else
			fail 'unknown arguments'
		end
	end
	if do_by_col % thin (or 'col'): do each column
		jj = 1:size(ob,2);
		out = fatrix2_subsref_colon2_vectors(ob, jj);
	else % wide (or 'row'): do each row
		ii = 1:size(ob,1);
		out = fatrix2_subsref_colon2_vectors(ob', ii)';
	end
	return

elseif islogical(sub2)
	tmp = col(sub2);
	if length(tmp) ~= size(ob,2);
		fail('bad column logical length')
	end
	jj = find(tmp);

elseif isnumeric(sub2)
	if isempty(sub2)
		out = [];
		return
	end

	jj = col(sub2);

	if 0 % eliminated this check to save time; could be an option
		bad = jj < 1 | jj > size(ob,2) | jj ~= round(jj);
		if any(bad)
			printm('bad column indeces:')
			minmax(jj(bad))
			printm('subsref problem')
			keyboard
		end
	end

else
	pr sub2
	fail('bad A(:,?)')
end

if numel(jj) == 1 % A(:,scalar)
	out = fatrix2_subsref_colon2_vectors(ob, jj); % return vector
else % A(:,array)
	out = fatrix2_subsref_list(ob, jj); % make a new fatrix2
end


% fatrix2_subsref_colon2_vectors()
% column vectors [A(:,?)] by matrix multiplication
function out = fatrix2_subsref_colon2_vectors(ob, jj)
z = zeros(size(ob,2), 1, 'single'); % trick: saves memory
x = z; x(jj(1)) = 1;
tmp = ob * x;
out = zeros(size(ob,1), length(jj), class(tmp));
out(:,1) = tmp;

for nn=2:length(jj)
	x = z;
	x(jj(nn)) = 1;
	out(:,nn) = ob * x;
end


% fatrix2_subsref_list()
% create modified fatrix2 based on index list for A(:,list)
function out = fatrix2_subsref_list(ob, list)

out = ob;
if isempty(ob.imask)
	out.imask = false([ob.idim 1]);
	out.imask(list) = true;
else
	tmp = false(ob.size(2), 1);
	tmp(list) = true;
	out.imask = iembed(ob, tmp);
end
out.size = [ob.size(1) numel(list)];
