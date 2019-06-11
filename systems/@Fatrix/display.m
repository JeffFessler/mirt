 function display(ob)
%function display(ob)
% "display" method for this class

name = inputname(1);
printf('"%s" is an object of class "%s":', name, class(ob))
ob = struct(ob);
disp(ob)

fnames = fieldnames(ob);
for ii=1:length(fnames)
	fname = fnames{ii};
	if isstruct(ob.(fname))
		printf('%s.%s :', inputname(1), fname)
		disp(ob.(fname))
	end
end
