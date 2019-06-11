function jf_whos_nan
% script to find variables in work space that have nan values

st = evalin('caller', 'whos');

for ii=1:length(st)
	cl = st(ii).class;
	name = st(ii).name;

	switch cl
	case {'double', 'single'}
		try
			tmp = sprintf('sum(isnan(%s(:)))', name);
			tmp = evalin('caller', tmp);
			if tmp
				printm('%d nan in %s', tmp, name);
			end
		catch
			printm(['unsure: ' name])
		end

	case 'struct'
		tmp = sprintf('%s', name);
		tmp = evalin('caller', tmp);

		jf_whos_nan_struct_do(tmp, name)
%		out = jf_struct_recurse([], tmp);
%		pr out
	end
end


function jf_whos_nan_struct_do(x, in_prefix)

for ix = 1:numel(x) % in case x is a struct array
	y = x(ix);
	names = fieldnames(y);
	for ii = 1:numel(names)
		name = names{ii};
		z = y.(name);

		if ii == 1
			prefix = in_prefix;
		else
			prefix = sprintf([in_prefix '(%d)'], ii);
		end
		
		if isempty(z)
			% ignore

		elseif isnumeric(z) || islogical(z)
			if any(isnan(z(:)))
				printm('%s.%s : %d NaN', prefix, name, ...
					sum(isnan(z(:))))

			end

		elseif isstruct(z)
			tmp = [prefix '.' names{ii}];
			jf_whos_nan_struct_do(z, tmp)

		end
	end
end
