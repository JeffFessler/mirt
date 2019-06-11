 function char_array_write(file, args)
%function char_array_write(file, args)
%write char array to file

fid = fopen(file, 'w');
if (fid == -1), error 'cannot open file', end
for ii=1:nrow(args)
	arg = args(ii,:);
%	arg(arg == 0) = []
	arg = sprintf('%s\n', arg);
%(arg == 0) = []
	fwrite(fid, arg, 'char');
end
if fclose(fid), error 'fclose?', end
