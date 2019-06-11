 function ir_display_struct(ob, varargin)
%function ir_display_struct(ob, varargin)
%|
%| display a structure concisely
%|
%| option
%|	'pre'		prefix
%|
%| Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if streq(ob, 'test'), ir_display_struct_test, return, end

arg.pre = ' ';
arg = vararg_pair(arg, varargin);

%pn = jf_protected_names;
% display(struct(pn))
% tmp = pn.struct_recurse(ob);

if ~isstruct(ob)
	fail('%s not a struct', class(ob))
end

printf([inputname(1) ':'])
%printf('pre = "%s"', arg.pre)
ir_display_struct_do(ob, arg.pre)


% ir_display_struct_do()
function ir_display_struct_do(ob, pre)

fnames = fieldnames(ob);

nmax = 1; % find max length field name
for ii=1:length(fnames)
	nmax = max(nmax, length(fnames{ii}));
end
%nmax = min(nmax, 40);

format = sprintf('%d', nmax);
format = [pre '%' format 's: %s'];

for ii=1:length(fnames)
	fname = fnames{ii};
	tmp = ob.(fname);

	if isempty(tmp)
		str = sprintf('[] (%s)', class(tmp));
	else
		switch class(tmp)
		case 'function_handle'
			str = func2str(tmp);
			if (str(1) ~= '@')
				str = ['@' str];
			end
		case 'char'
			str = ['"' tmp '"'];
		case 'struct'
			str = '(struct)';
		case 'cell'
			str = ['{' num2str(size(tmp)) '}'];
		case {'double', 'single', 'logical'}
			if numel(tmp) == 1
				str = num2str(tmp);
			else
				str = ['[' num2str(size(tmp)) ']'];
				if numel(tmp) < 5
					str = [str ' : ' num2str(tmp(:)')];
				end
			end
			str = [str sprintf(' (%s)', class(tmp))];
		otherwise
			str = sprintf('? (%s)', class(tmp));
		end
	end

	printf(format, fname, str)
	if isstruct(tmp)
		white = sprintf('%d', nmax+2); % 1 more because of ':'
		white = [pre '%' white 's'];
		white = sprintf(white, ' ');
		ir_display_struct_do(tmp, white);
	end
end

function ir_display_struct_test
s.a = 'a';
s.b = 2;
s.c = {1,2,3};
s.d = ones(3);
s.e = 1+2i;
s.f1 = @(x) x+1;
s.f2 = @ir_display_struct;
s.g = complex(1);
t.s = 'sub';
u.v = 'subsub';
t.u = u;
s.t = t;
ir_display_struct(s)
