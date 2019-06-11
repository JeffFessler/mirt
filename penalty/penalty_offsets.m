  function offsets = penalty_offsets(offsets, isize)
%|function offsets = penalty_offsets(offsets, isize)
%|
%| determine the 'offsets' argument needed by regularizers
%|
%| in
%|	offsets		empty or string or array
%|			'1d' [1] (1d default)
%|			'2d:hvd' [1 nx nx+1 nx-1] (2d default)
%|			'3d:hvu' [1 nx nx*ny] (3d default)
%|			'3d:26'	all 13 neighbors in 3d
%|	isize [N]	image size
%| out
%|	offsets [L]	penalty offsets, e.g., [1 nx nx+1 nx-1]
%|
%| Copyright 2006-12-6, Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

if isempty(offsets)
	switch length(isize)
	case 1
		if isize(1) <= 1
			offsets = [];
			warn 'single pixel so no offset'
		else
			offsets = [1];
		end

	case 2
		if isize(2) == 1 % 1d
			if isize(1) <= 1
				offsets = [];
				warn 'single pixel so no offset'
			else
				offsets = [1];
			end
		else % 2d
			nx = isize(1);
			offsets = [1 nx nx+1 nx-1];
		end

	case 3 % bare-bones 3D
		nx = isize(1);
		ny = isize(2);
		offsets = [1 nx nx*ny];

	otherwise
		error 'only 2D and 3D done'
	end

elseif ischar(offsets)
	switch offsets
	case {'ident', 'identity', 'I'}
		offsets = [0]; % trick

	case '1d'
		if max(isize) <= 1
			offsets = [];
			warn 'single pixel so no offset'
		else
			offsets = [1];
		end

	case {'2d,hv', '2d:hv'}
		nx = isize(1);
		ny = isize(2);
		if nx == 1 && ny == 1
			offsets = [];
			warn 'single pixel so no offset'
		elseif nx == 1 || ny == 1 % 1d
			offsets = [1];
		else
			offsets = [1 nx];
		end

	case {'2d,hvd', '2d:hvd'}
		nx = isize(1);
		ny = isize(2);
		if nx == 1 && ny == 1
			offsets = [];
			warn 'single pixel so no offset'
		elseif nx == 1 || ny == 1
			offsets = [1];
		else
			offsets = [1 nx nx+1 nx-1];
		end

	case '3d:hvu'
		nx = isize(1);
		ny = isize(2);
		offsets = [1 nx nx*ny];
		if nx == 1 || ny == 1
			warn 'not done'
		end

	case '3d:26' % all 26 neighbors (13 pairs)
		if length(isize) ~= 3, error '3d:26 expects 3d image', end
		nx = isize(1);
		ny = isize(2);
		nz = isize(3);
		offsets = [1 nx+[0 1 -1] ...
			nx * ny + col(outer_sum([-1:1]', [-1:1]*nx))'];
		if nx == 1 || ny == 1 || nz == 1
			warn 'not done'
		end

	otherwise
		fail('bad offsets string "%s"', offsets)
	end
end

% if offsets has a zero, it must be just a single zero (for identity)
if any(offsets == 0) && any(offsets)
	error 'identity is offsets=[0] only'
end
