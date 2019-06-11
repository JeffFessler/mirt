 function x = ir_read_mat(file, pick)
%function x = ir_read_mat(file, pick)
%|
%| Read a matlab matrix from file, e.g., 'data.mat'
%| If pick is an integer, it indexes which field in the file is desired.
%| Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if ~isvar('pick') || isempty(pick), pick = 0; end

s = load(file);
names = fieldnames(s);
name = names{1};
if length(names) > 1
	if ~isa(pick, 'numeric')
		error 'only numeric picking done now'
	end
	if pick < 1
		disp(names)
		warn('using default first array: %s', name)
	elseif pick > length(names)
		error 'pick too large'
	else
		name = names{pick};
		printm('Using array %s', name)
	end
end
x = s.(name);
