  function [B Bgx Bgy Bgz] = makeB(ig, kg, varargin)
%|function [B Bgx Bgy Bgz] = makeB(ig, kg, [?])
%|
%| Create a Fatrix B operation of warping based on cubic B-spline
%|
%| in:
%|	'ig'		image geometry 
%|	'kg'		knot geometry
%|
%| out:
%|	B		fatrix operation for warping
%|	Bgx		fatrix operation for warping
%|	Bgy		fatrix operation for warping
%|	Bgz		fatrix operation for warping
%|
%| Copyright August 2006, Se Young Chun and Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

argB.power = 0;
argB.ig = ig;
argB.kg = kg;
[argB.kernelx argB.kernelgx] = makeKernel(3, kg.mx);
[argB.kernely argB.kernelgy] = makeKernel(3, kg.my);
if (kg.is3 == 1) 
	[argB.kernelz argB.kernelgz] = makeKernel(3, kg.mz);
end

if (nargin == 2)
	B = Fatrix([ig.np kg.np], argB, 'forw', @makeB_forward, ...
		'back', @makeB_transpose, 'power', @makeB_power);
	Bgx = Fatrix([ig.np kg.np], argB, 'forw', @makeBgx_forward, ...
		'back', @makeBgx_transpose);
	Bgy = Fatrix([ig.np kg.np], argB, 'forw', @makeBgy_forward, ...
		'back', @makeBgy_transpose);
	if (kg.is3 == 1) 
		Bgz = Fatrix([ig.np kg.np], argB, 'forw', @makeBgz_forward, ...
			'back', @makeBgz_transpose);
	end

else
	if (kg.is3 == 1) 
		argB.WarpCell = varargin(1:3);
		Bgz = Fatrix([ig.np kg.np],argB,'forw',@makeBgz_forwardwarp, ...
			'back', @makeBgz_transposewarp);
	else
		argB.WarpCell = varargin(1:2);
	end
	B = Fatrix([ig.np kg.np], argB, 'forw', @makeB_forwardwarp, ...
		'back', @makeB_transposewarp);
	Bgx = Fatrix([ig.np kg.np], argB, 'forw', @makeBgx_forwardwarp, ...
		'back', @makeBgx_transposewarp);
	Bgy = Fatrix([ig.np kg.np], argB, 'forw', @makeBgy_forwardwarp, ...
		'back', @makeBgy_transposewarp);
end


	
%%%
function DFM = makeB_forward(argB, alpha)
alpha = double(reshape(alpha, argB.kg.dim));
if (argB.kg.is3 == 0) % 2D case
	DFM = BsplCo2ValZeroFilt(alpha, [argB.ig.nx, argB.ig.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], {argB.kernelx, argB.kernely});
else
	DFM = BsplCo2ValZeroFilt(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], ...
		{argB.kernelx, argB.kernely, argB.kernelz});
end
DFM = single(DFM(:));



%%%
function ALP = makeB_transpose(argB, imgco)
imgco = double(reshape(imgco, argB.ig.dim));
if (argB.kg.is3 == 0) % 2D case
	ALP = BsplCo2ValTranZeroFilt(imgco, [argB.kg.nx argB.kg.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], {argB.kernelx, argB.kernely});
else
	ALP = BsplCo2ValTranZeroFilt(imgco,[argB.kg.nx argB.kg.ny argB.kg.nz],...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], ...
		{argB.kernelx, argB.kernely, argB.kernelz});
end
ALP = single(ALP(:));



%%%
function DFM = makeB_forwardwarp(argB, alpha)
alpha = single(reshape(alpha, argB.kg.dim));
if (argB.kg.is3 == 0) % 2D case
        DFM = BsplCo2ValZero(alpha, [argB.ig.nx, argB.ig.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], argB.WarpCell);
else
        DFM = BsplCo2ValZero(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
end
DFM = DFM(:);



%%%
function ALP = makeB_transposewarp(argB, imgco)
imgco = single(reshape(imgco, argB.ig.dim));
if (argB.kg.is3 == 0) % 2D case
	ALP = BsplCo2ValTranZero(imgco, [argB.kg.nx, argB.kg.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], argB.WarpCell);
else
	ALP = BsplCo2ValTranZero(imgco, [argB.kg.nx, argB.kg.ny argB.kg.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
end
ALP = ALP(:);



%%%
function DFM = makeBgx_forward(argB, alpha)
alpha = double(reshape(alpha, argB.kg.dim));
if (argB.kg.is3 == 0) % 2D case
	DFM = BsplCo2ValZeroFilt(alpha, [argB.ig.nx, argB.ig.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], {argB.kernelgx, argB.kernely});
else
	DFM = BsplCo2ValZeroFilt(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], ...
		{argB.kernelgx, argB.kernely, argB.kernelz});
end
DFM = single(DFM(:));



%%%
function ALP = makeBgx_transpose(argB, imgco)
imgco = double(reshape(imgco, argB.ig.dim));
if (argB.kg.is3 == 0) % 2D case
	ALP = BsplCo2ValTranZeroFilt(imgco, [argB.kg.nx argB.kg.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], ...
		{fliplr(argB.kernelgx), fliplr(argB.kernely)});
else
	ALP = BsplCo2ValTranZeroFilt(imgco,[argB.kg.nx argB.kg.ny argB.kg.nz],...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], ...
		{fliplr(argB.kernelgx), fliplr(argB.kernely), ...
		fliplr(argB.kernelz)});
end
ALP = single(ALP(:));



%%%
function DFM = makeBgx_forwardwarp(argB, alpha)
alpha = single(reshape(alpha, argB.kg.dim));
if (argB.kg.is3 == 0) % 2D case
        DFM = BsplCo2GdXZero(alpha, [argB.ig.nx, argB.ig.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], argB.WarpCell);
else
        DFM = BsplCo2GdXZero(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
end
DFM = DFM(:);



%%%
function ALP = makeBgx_transposewarp(argB, imgco)
imgco = single(reshape(imgco, argB.ig.dim));
if (argB.kg.is3 == 0) % 2D case
	ALP = BsplCo2GdXTranZero(imgco, [argB.kg.nx, argB.kg.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], argB.WarpCell);
else
	ALP = BsplCo2GdXTranZero(imgco, [argB.kg.nx, argB.kg.ny argB.kg.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
end
ALP = ALP(:);



%%%
function DFM = makeBgy_forward(argB, alpha)
alpha = double(reshape(alpha, argB.kg.dim));
if (argB.kg.is3 == 0) % 2D case
	DFM = BsplCo2ValZeroFilt(alpha, [argB.ig.nx, argB.ig.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], {argB.kernelx, argB.kernelgy});
else
	DFM = BsplCo2ValZeroFilt(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], ...
		{argB.kernelx, argB.kernelgy, argB.kernelz});
end
DFM = single(DFM(:));



%%%
function ALP = makeBgy_transpose(argB, imgco)
imgco = double(reshape(imgco, argB.ig.dim));
if (argB.kg.is3 == 0) % 2D case
	ALP = BsplCo2ValTranZeroFilt(imgco, [argB.kg.nx argB.kg.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], ...
		{fliplr(argB.kernelx), fliplr(argB.kernelgy)});
else
	ALP = BsplCo2ValTranZeroFilt(imgco,[argB.kg.nx argB.kg.ny argB.kg.nz],...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], ...
		{fliplr(argB.kernelx), fliplr(argB.kernelgy), ...
		fliplr(argB.kernelz)});
end
ALP = single(ALP(:));



%%%
function DFM = makeBgy_forwardwarp(argB, alpha)
alpha = single(reshape(alpha, argB.kg.dim));
if (argB.kg.is3 == 0) % 2D case
        DFM = BsplCo2GdYZero(alpha, [argB.ig.nx, argB.ig.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], argB.WarpCell);
else
        DFM = BsplCo2GdYZero(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
end
DFM = DFM(:);



%%%
function ALP = makeBgy_transposewarp(argB, imgco)
imgco = single(reshape(imgco, argB.ig.dim));
if (argB.kg.is3 == 0) % 2D case
	ALP = BsplCo2GdYTranZero(imgco, [argB.kg.nx, argB.kg.ny], ...
		[argB.kg.offset_x argB.kg.offset_y], ...
		[argB.kg.mx argB.kg.my], argB.WarpCell);
else
	ALP = BsplCo2GdYTranZero(imgco, [argB.kg.nx, argB.kg.ny argB.kg.nz], ...
		[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
		[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
end
ALP = ALP(:);



%%%
function DFM = makeBgz_forward(argB, alpha)
alpha = double(reshape(alpha, argB.kg.dim));
DFM = BsplCo2ValZeroFilt(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
	[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
	[argB.kg.mx argB.kg.my argB.kg.mz], ...
	{argB.kernelx, argB.kernely, argB.kernelgz});
DFM = single(DFM(:));



%%%
function ALP = makeBgz_transpose(argB, imgco)
imgco = double(reshape(imgco, argB.ig.dim));
ALP = BsplCo2ValTranZeroFilt(imgco,[argB.kg.nx argB.kg.ny argB.kg.nz],...
	[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
	[argB.kg.mx argB.kg.my argB.kg.mz], ...
	{fliplr(argB.kernelx), fliplr(argB.kernely), ...
	fliplr(argB.kernelgz)});
ALP = single(ALP(:));



%%%
function DFM = makeBgz_forwardwarp(argB, alpha)
alpha = single(reshape(alpha, argB.kg.dim));
DFM = BsplCo2GdZZero(alpha, [argB.ig.nx, argB.ig.ny argB.ig.nz], ...
	[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
	[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
DFM = DFM(:);



%%%
function ALP = makeBgz_transposewarp(argB, imgco)
imgco = single(reshape(imgco, argB.ig.dim));
ALP = BsplCo2GdZTranZero(imgco, [argB.kg.nx, argB.kg.ny argB.kg.nz], ...
	[argB.kg.offset_x argB.kg.offset_y argB.kg.offset_z], ...
	[argB.kg.mx argB.kg.my argB.kg.mz], argB.WarpCell);
ALP = ALP(:);



%%%
function ob = makeB_power(ob, sup)

if (sup == 1)

elseif ((sup == 2) && (ob.arg.power == 0))
	ob.arg.power = 1;
else
	error('only squares of each element are supported');
end
ob = ob;
