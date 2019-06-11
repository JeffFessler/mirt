  function out = ir_nargout(fun)
%|function out = ir_nargout(fun)
%|
%| version of nargout() that tries to work with freemat and octave,
%| basically by defaulting to "1"
%|
%| 2012-09-23, Jeff Fessler

if isfreemat % freemat: does not like nargout()
	out = 1;

elseif ir_is_octave

%	http://hg.savannah.gnu.org/hgweb/octave/rev/20cb178716ee
	if ~ischar(fun) % octave 3.6.3 does not like function handles
		fun_in = fun;
		fun = func2str(fun);
	end

	if fun(1) == '@'
		out = 1; % assume @(x) f(x) returns 1 output
	else
		try
			out = nargout(str);
		catch % octave does not like private functions
%			warn('speculating nargout("%s")=1', str)
			out = 1;
		end
	end

else % matlab
	out = nargout(fun);
end
