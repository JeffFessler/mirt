 function label(varargin)
% usage:
%	label x string
%	label y string
%	label tex ...	todo

narg = length(varargin);
if narg < 2
	error 'need x or y'
end

props = {'interpreter', 'none'};

if streq(varargin{1}, 'x')
	arg2 = varargin{2};
	xlabel(arg2, props{:})

elseif streq(varargin{1}, 'y')
	arg2 = varargin{2};
	ylabel(arg2, props{:})

else
	error 'not done'
end

function xlabeltex(arg)
%function xlabeltex(arg)
% xlabel embedded in \tex{} for psfrag
xlabel(['\tex[t][t]{' arg '}'], 'interpreter', 'none')

 function ylabeltex(arg)
%function ylabeltex(arg)
% ylabel embedded in \tex{} for psfrag
ylabel(['\tex[B][B]{' arg '}'], 'interpreter', 'none')
