 function [W Wgx Wgy Wgz] = makeW(Bcell, Acell)
%function [W Wgx Wgy Wgz] = makeW(Bcell, Acell)
%
% Create a Fatrix W operation of deformation based on rect basis 
%
% in:
%	'Bcell'		array of cells of fatrix operation B for warping {x}{y}{z} 
%	'Acell'		array of cells of warping coefficients Alpha {x}{y}{z} 
%
% out:
%	'W'		fatrix operation for deformation
%	'Wgx'		fatrix operation for deformation
%	'Wgy'		fatrix operation for deformation
%	'Wgz'		fatrix operation for deformation
%
% Copyright August 2006, Se Young Chun and Jeff Fessler, University of Michigan 

for i = 1 : length(Bcell)
	argW.B{i} = Bcell{i};
	argW.BA{i} = Bcell{i}*single(Acell{i});
end

argW.deg = 3;

W = Fatrix([Bcell{1}.arg.ig.np Bcell{1}.arg.ig.np], argW, 'forw', ...
	@makeW_forward, 'back', @makeW_transpose, 'ufun', @makeW_warpgrid);
Wgx = Fatrix([Bcell{1}.arg.ig.np Bcell{1}.arg.ig.np], argW, 'forw', ...
	@makeWgx_forward, 'back', @makeWgx_transpose);
Wgy = Fatrix([Bcell{1}.arg.ig.np Bcell{1}.arg.ig.np], argW, 'forw', ...
	@makeWgy_forward, 'back', @makeWgy_transpose);
if (length(Bcell) == 3)
	Wgz = Fatrix([Bcell{1}.arg.ig.np Bcell{1}.arg.ig.np], argW, 'forw', ...
		@makeWgz_forward, 'back', @makeWgz_transpose);
end


%%%
function IM = makeW_forward(argW, u)
u = single(reshape(u, argW.B{1}.arg.ig.dim));
if (argW.B{1}.arg.ig.is3 == 0) % 2D case
        IM = BsplCo2ValMirr(u, argW.B{1}.arg.ig.dim, [0 0], [1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:)}, argW.deg);
else
        IM = BsplCo2ValMirr(u, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
end
IM = IM(:);



%%%
function CO = makeW_transpose(argW, IM)
IM = single(reshape(IM, argW.B{1}.arg.ig.dim));
if (argW.B{1}.arg.ig.is3 == 0) % 2D case
        CO = BsplCo2ValTranMirr(IM, argW.B{1}.arg.ig.dim, [0 0], [1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:)}, argW.deg);
else
        CO = BsplCo2ValTranMirr(IM, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
end
CO = CO(:);



%%%
function GDx = makeWgx_forward(argW, u)
u = single(reshape(u, argW.B{1}.arg.ig.dim));
if (argW.B{1}.arg.ig.is3 == 0) % 2D case
        GDx = BsplCo2GdXMirr(u, argW.B{1}.arg.ig.dim, [0 0], [1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:)}, argW.deg);
else
        GDx = BsplCo2GdXMirr(u, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
end
GDx = GDx(:);



%%%
function CO = makeWgx_transpose(argW, IM)
IM = single(reshape(IM, argW.B{1}.arg.ig.dim));
if (argW.B{1}.arg.ig.is3 == 0) % 2D case
        CO = BsplCo2GdXTranMirr(IM, argW.B{1}.arg.ig.dim, [0 0], [1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:)}, argW.deg);
else
        CO = BsplCo2GdXTranMirr(IM, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
end
CO = CO(:);



%%%
function GDy = makeWgy_forward(argW, u)
u = single(reshape(u, argW.B{1}.arg.ig.dim));
if (argW.B{1}.arg.ig.is3 == 0) % 2D case
        GDy = BsplCo2GdYMirr(u, argW.B{1}.arg.ig.dim, [0 0], [1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:)}, argW.deg);
else
        GDy = BsplCo2GdYMirr(u, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
end
GDy = GDy(:);



%%%
function CO = makeWgy_transpose(argW, IM)
IM = single(reshape(IM, argW.B{1}.arg.ig.dim));
if (argW.B{1}.arg.ig.is3 == 0) % 2D case
        CO = BsplCo2GdYTranMirr(IM, argW.B{1}.arg.ig.dim, [0 0], [1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:)}, argW.deg);
else
        CO = BsplCo2GdYTranMirr(IM, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
		, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
end
CO = CO(:);



%%%
function GDz = makeWgz_forward(argW, u)
u = single(reshape(u, argW.B{1}.arg.ig.dim));
GDz = BsplCo2GdZMirr(u, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
	, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
GDz = GDz(:);



%%%
function CO = makeWgz_transpose(argW, IM)
IM = single(reshape(IM, argW.B{1}.arg.ig.dim));
CO = BsplCo2GdZTranMirr(IM, argW.B{1}.arg.ig.dim, [0 0 0], [1 1 1]...
	, {argW.BA{1}(:), argW.BA{2}(:), argW.BA{3}(:)}, argW.deg);
CO = CO(:);



%%%
function warpgrid = makeW_warpgrid(W, space)
if nargin < 2
	space = 8; % jf
end

grid = zeros(W.arg.B{1}.arg.ig.data.dim, 'single');

if W.arg.B{1}.arg.ig.data.is3 == 0 % B:2D case
        grid(1:space:end, :, :) = 1;
        grid(:, 1:space:end, :) = 1;

else % B:3D case
        grid(1:space:end, :, :) = 1;
        grid(:, 1:space:end, :) = 1;
        grid(:, :, 1:space:end) = 1;
end
Cgrid = BsplVal2CoMirr(grid, 3);
IM = BsplCo2ValMirr(Cgrid, W.arg.B{1}.arg.ig.dim, [0 0], [1 1]...
		, {W.arg.BA{1}(:), W.arg.BA{2}(:)}, 3);
warpgrid = reshape(IM, W.arg.B{1}.arg.ig.data.dim);
