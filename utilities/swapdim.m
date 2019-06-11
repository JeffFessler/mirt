 function y = swapdim(x, dim1, dim2)
%function y = swapdim(x, dim1, dim2)
%	swap the given dimensions of x
%	eg. swapdim(x, 1, 2) transposes each slice

error 'obsolete: use permute()'

if nargin < 3
	help swapdim
	error args
end

	t = size(x);
	t([dim1 dim2]) = t([dim2 dim1]);
	y = zeros(t);

	if dim1 == 1 && dim2 == 2 && ndims(x) == 3
		for iz=1:size(x,3)
			y(:,:,iz) = x(:,:,iz)';
		end

	elseif dim1 == 1 && dim2 == 2 && ndims(x) == 4
		for iz=1:size(x,3)
			for i4=1:size(x,4)
				y(:,:,iz,i4) = x(:,:,iz,i4)';
			end
		end

	elseif dim1 == 2 && dim2 == 3 && ndims(x) == 3
		for iz=1:size(x,3)
			y(:,iz,:) = x(:,:,iz);
		end

	else
		error('not done')
	end
