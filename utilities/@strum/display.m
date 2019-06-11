 function display(ob)
%function display(ob)
% "display" method for this class

name = inputname(1);
printm('"%s" is an object of class "%s":', name, class(ob))
ob = struct(ob);
ir_display_struct(ob);

if ~isempty(ob.meth) && ~isempty(ob.docs)
	printf('\n methods with comments:')
	display_method(ob.meth, ob.docs)
end
return

fnames = fieldnames(ob);
for ii=1:length(fnames)
	fname = fnames{ii};
	if isstruct(ob.(fname))
		printf('%s.%s :', inputname(1), fname)
		if streq(fname, 'meth') && ~isempty(ob.docs)
			display_method(ob.meth, ob.docs)
		else
%			ir_display(ob.(fname))
		end
	end
end


function display_method(meth, docs)
mnames = fieldnames(meth);

nmax = 1; % find max length method name
for im=1:length(mnames)
	nmax = max(nmax, length(mnames{im}));
end
%nmax = min(nmax, 40);

hmax = 1; % find max length handle name
for im=1:length(mnames)
	hmax = max(hmax, length( func2str(meth.(mnames{im})) ));
end
%hmax = min(hmax, 40);
%pr '[nmax hmax]'

format = [' %' sprintf('%d', nmax) 's: @%' sprintf('%d', hmax) 's %s'];

for im=1:length(mnames)
	tmp = func2str(meth.(mnames{im}));
	if tmp(1) == '@', tmp = tmp(2:end); end % implicit functions
	printf(format, mnames{im}, tmp, docs{im})
end


function ir_display(ob)
display(class(ob))

switch class(ob)
case 'struct'
	ir_display_struct(ob)
return

otherwise
	disp(ob)
end
