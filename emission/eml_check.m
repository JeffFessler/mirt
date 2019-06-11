 function eml_check(yi, ci, ri, varargin)
%function eml_check(yi, ci, ri [, 'os', nblock])
% check the validity of emission data from Poisson transmission model
% if 'os' option, then check that data is a 2d matrix, since only the
% second index is subsetized.  so for SPECT or 3D PET this means one must
% reshape each projection view into a column, e.g., before calling OS recon do:
% reshape(yi, [size(yi,1)*size(yi,2) size(yi,3)])
% or just:
% reshaper(yi, '2d') 

if nargin < 3, error(mfilename), end

if iscell(yi)
%	if length(yi) ~= length(ci)
	for ii=1:length(yi)
		eml_check_data(yi{ii}, ci{ii}, ri{ii});
	end
return
end


eml_check_data(yi, ci, ri)

is_fbp = 0;
if nargin > 3
	arg = varargin{1};
	if streq(arg, 'fbp')
		is_fbp = 1;
	elseif streq(arg, 'os')
		if ndims(yi) ~= 2
			error(['data must have 2 dimensions for OS.  see ' mfilename])
		end
%		if size(yi,2) ~= varargin{2}
%		end
	else
		error 'bug'
	end
end

if ndims(yi) > 2 && ~is_fbp
	error 'most routines need column yi(:), except OS routines need "2D" yi, using "reshaper()"'
end


%
% eml_check_data(yi, ci, ri)
%
function eml_check_data(yi, ci, ri)

if any(size(ci) > 1)
	if ~isequal(size(yi), size(ci)), error 'yi vs ci size', end
end
if any(size(ri) > 1)
	if ~isequal(size(yi), size(ri)), error 'yi vs ri size', end
end

if any(ci(:) < 0), error 'ci''s must be nonnegative', end
if any(yi(:) < 0), error 'yi''s must be nonnegative', end
if any(ri(:) < 0), error 'ri''s must be nonnegative', end

if any(yi(:) & ~(ci(:) + ri(:))), error 'model mismatch', end
