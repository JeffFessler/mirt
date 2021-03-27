function  [A_mc] = Asense(A, sens)
%function [A_mc] = Asense(A, sens)
%|
%| Add sensitivity maps to a fatrix
%|
%| Inputs
%|   A      [nd np]    system matrix fatrix object
%|   sens   [nc 1]     cell containing sensitivty maps for each coil
%|     OR   [nx ny nz nc]    matrix containing 3d sensitivity maps 
%|     OR   [nx ny nc]       matrix containing 2d sensitivity maps 
%|   
%| Outputs
%|   A_mc   [nd*nc np]
%|
%|
%|
%|  Sizes:
%|    nd - number of data points for single coil
%|    nc - number of coils
%|    np - number of image pixels
%|
%| Melissa Haskell, April 2020

% check inputs
if nargin < 2, help(mfilename), error('Not enough input arguments.'), end

% if sens is a matrix, convert to cell
if ~iscell(sens)
    size_v = size(sens);
    if numel(size_v) == 2
        nc = 1;
        sens = squeeze(mat2cell(sens, size_v(1), size_v(2), ones(nc,1)));
    else
        nc = size_v(end);
    end
    if numel(size_v) == 3
        sens = squeeze(mat2cell(sens, size_v(1), size_v(2), ones(nc,1)));
    elseif numel(size_v) == 4
        sens = squeeze(mat2cell(sens, size_v(1), size_v(2), size_v(3), ones(nc,1)));
    end
end

% check that size of sense maps and A match up
np = numel(sens{1});
if size(A,2) ~= np
	fail('# of columns in A must be the same size as each sens map')
end

% create cell with each elements of A scaled by the sens map, then
% concatenate
nc = numel(sens);
A_mc_cell = cell(nc,1);
for ii = 1:nc
    A_mc_cell{ii} = A*Gdiag(sens{ii});
end
A_mc = cat(1,A_mc_cell{:});


end




