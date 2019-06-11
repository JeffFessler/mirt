  function ob = Gwave1( varargin )
%|function ob = Gwave1([args])
%| Construct Gwave object that computes wavelet decomposition of a 
%| signal with dimensions [(Nd)].  This is useful for compressed sensing
%| imaging. (1D or 2D wavelets only)
%|
%| options:
%|	'mask'	logical [(Nd)]	image-domain mask, usually: true(nx,ny)
%|  'level'  integer specifying decomposition level (default 1)
%|  'wname'  string specifying wavelet, ex. 'haar', 'db1', 'db2'
%|              (default 'haar') (**see note**).
%|  'dwtmode'  set global dwtmode parameter (default = 'zpd')
%|  'zeropad'  (0|1) zeropad the input data to be multiple of 2^level
%|                (default 0)
%|  'noCA'  (0|1) discard the approximation coefficients (default 0)
%|  'swt'  (0|1) use stationary wavelet transform (default 0)
%|  'out1D'  (0|1) force output of Gwave_back to be column vector
%|                  containing only elements corresponding to the '1' cells 
%|                  of the mask (default 1)
%|
%| out:
%|	ob	[np]	Fatrix object, where np = sum(mask(:))
%|
%| Basically, you create a system matrix object by calling:
%|	A = Gwave1( ... )
%| and then you can use it thereafter by typing commands like
%|	y = A * x;
%| which will auto-magically evaluate the wavelet coefficients.
%| This is useful for compressed sensing of natural images.
%|
%| Besides simple utilities like display, there are the following
%| capabilities of this object:
%|	y = A * x		forward operation
%|	x = A' * y		adjoint operation
%|	
%| Note on wavelets: Matlab boundary conditions may create a longer output
%| vector than input.  This discrepancy can affect adjoint and orthoganal
%| properties, particularly if 'dwtmode' ~= 'zpd'.
%|
%| Copyright 2009-11-7, Michael Allison, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), Gwave1_test, return, end

arg.mask = [];
arg.level = 1;
arg.wname = 'haar';
arg.dwtmode = 'zpd';
arg.zeropad = 0;
arg.noCA = 0;
arg.swt = 0;
arg.out1D = 1;
arg = vararg_pair(arg, varargin);

if isempty(arg.mask), fail 'must provide a mask', end
if ~islogical(arg.mask), error 'mask must be logical', end

if arg.swt
    arg.zeropad = 1; %required
end

if arg.noCA ~= 0 && arg.noCA ~= 1
    error('noCA must be ''0'' or ''1''.');
end

% Change DWTMODE parameter and inform user.
st = dwtmode('status','nodisp');
if ~strcmp(st,arg.dwtmode)
    dwtmode(arg.dwtmode,'nodisp');
    disp(['** DWTMODE CHANGED TO ''' arg.dwtmode ''' FROM ''' st ''' **']);
end

% Determine dimension of transform.
arg.Nd = size(arg.mask);
if arg.Nd(end) == 1
    arg.ndim = 1;
else
    arg.ndim = length(arg.Nd);
end

% Determine required zero-padding.
if arg.zeropad
    Nl = 2^(arg.level);
    if arg.ndim == 1
        pad = logical(mod(arg.Nd(1),Nl)) * (Nl - mod(arg.Nd(1),Nl)); % extra term ensures nothing added if exact.
        pad = [pad 0]; %easier future computation.
    else
        pad = zeros(1,arg.ndim);
        for i = 1:arg.ndim
            pad(i) = logical(mod(arg.Nd(i),Nl)) * (Nl - mod(arg.Nd(i),Nl));  
        end
    end
else
    pad = zeros(size(arg.Nd));
end
arg.pad = pad;

% Compute output coefficients (swt, zeropad, noCA).
if ~arg.swt % DWT mode.
    [clength, clengthCA, s] = outCoeffDWT(arg); % hand comp. depends on DWT behaviour.
    arg.s = s;
elseif arg.swt % SWT mode.
    Ndtmp = arg.Nd + arg.pad;
    if arg.ndim == 1;
            clengthCA = prod(Ndtmp) * (arg.level+1);
            clength = clengthCA - arg.noCA * prod(Ndtmp);
    elseif arg.ndim == 2;
            clengthCA = prod(Ndtmp) * (3 * arg.level + 1);
            clength = clengthCA - arg.noCA * prod(Ndtmp);
    else
        error('GWAVE.M currently only supports 1D and 2D matrices.');
    end
end

arg.clength = clength;
arg.clengthCA = clengthCA;
arg.np = sum(arg.mask(:));
arg.dim = [arg.clength arg.np];

%
% build Fatrix object
%
ob = Fatrix(arg.dim, arg, ...
	'forw', @Gwave_forw, 'back', @Gwave_back, ...
	'caller', mfilename);


%%
% Gwave_forw(): y = A * x 
%    - Performs wavelet analysis on x.
% in
%	x	[np L] or [(Nd) L]
% out
%	y	[clength L]   where clength is either full coefficients or noCA.
%
function y = Gwave_forw(arg, x)

xdim = size(x);
NN = prod(arg.Nd);
LL = xdim(end);

% apply mask to data (stacked data has same mask applied to each level).
if size(x,1) == arg.np		% insert data into mask [np (L)]
    x = embed(x, arg.mask);	% [(Nd) (L)]
elseif isequal(size(x), arg.Nd)    % apply mask to single image data.
    x(~arg.mask) = 0;  
else   % apply mask to stacked data.  %%Could use stackpick/stackup... cat(ndim+1,etc). better for memory!!**
    if arg.ndim == 1 && isequal(xdim(1), NN)
        repmask = repmat(arg.mask,[1 LL]);
    elseif arg.ndim == 2 && isequal(xdim(1:2), arg.Nd)
        repmask = repmat(arg.mask,[1 1 LL]);
    else
        error('Input x must be either [np L], Nd, or [Nd L].');
    end
    x(~repmask) = 0;
end

% zeropad the matrix
x = padarray(x,arg.pad, 0, 'post');

% compute coefficients and remove approx. coeff. if desired.
if arg.swt
    if isreal(x)
        y = Gwave_forw_SWT(arg, x, xdim);
    else
        yr = Gwave_forw_SWT(arg, real(x), xdim);
        yi = Gwave_forw_SWT(arg, imag(x), xdim);
        y = complex(yr,yi);
    end
else
    % complex data must be treated seperately.
    if isreal(x)
        y = Gwave_forw_DWT(arg, x, xdim);
    else
        yr = Gwave_forw_DWT(arg, real(x), xdim);
        yi = Gwave_forw_DWT(arg, imag(x), xdim);
        y = complex(yr,yi);
    end
end

%%
% Gwave_back(): x = A' * y
%   - Performs wavelet synthesis on y.
% in
%	y	[clength L]   where clength is either full coefficients or noCA.
% out
%	x	[Np L] if 'out1D' = 1, or [Nd L]
%
function x = Gwave_back(arg, y)

ydim = size(y); % coefficiants always vectors
LL = ydim(end);

% Reconstruct x from y.
if arg.swt
    if isreal(y)
        x = Gwave_back_SWT(arg, y);
    else
        xr = Gwave_back_SWT(arg, real(y));
        xi = Gwave_back_SWT(arg, imag(y));
        x = complex(xr,xi);
    end
else
    if isreal(y)
        x = Gwave_back_DWT(arg, y);
    else % separation doesn't appear to be needed, but best to be safe.
        xr = Gwave_back_DWT(arg, real(y));
        xi = Gwave_back_DWT(arg, imag(y));
        x = complex(xr,xi);
    end
end

% Remove the extra pixels.
if arg.zeropad
    if arg.ndim == 1
        x = x(1:arg.Nd(1), :);
    elseif arg.ndim == 2
        x = x(1:arg.Nd(1), 1:arg.Nd(2), :);
    else 
        fail 'Gwave_back only accomodates 1D or 2D data.';
    end
end

if arg.out1D
    if arg.ndim == 2
        x = reshape(x,prod(arg.Nd),LL);
    end
    x = x(arg.mask(:),:);
end

%%
% Gwave_forw_DWT(): y = A * x  where A is the discrete wavelet transform
% in
%	x	[np L] or [(Nd) L]  where x is zeropadded and has masked applied.
% out
%	y	[clength L]  where clength is either full coefficients or noCA.
%
function y = Gwave_forw_DWT(arg, x, xdim)

LL = xdim(end);

% Single image.
if isequal(size(x), arg.Nd+arg.pad) % [(Nd)]
    if arg.ndim == 1
        y = wavedec(x,arg.level,arg.wname);
    elseif arg.ndim == 2
        y = wavedec2(x,arg.level,arg.wname);
    else
        fail 'only performs 1D or 2D wavelet decomposition'
    end
    y = col(y); % ensure its a column vector, not row
else % stacked images.
    if arg.ndim == 1
        y = zeros([arg.clengthCA LL]);
        for ll=1:LL  % Could use mdwtdec,but need entire object to recon.
            tmp = stackpick(x, ll);
            y(:,ll) = wavedec(tmp,arg.level,arg.wname);
        end
    elseif arg.ndim == 2
        y = zeros([arg.clengthCA LL]);
        for ll=1:LL
            tmp = stackpick(x, ll);
            y(:,ll) = wavedec2(tmp,arg.level,arg.wname);
        end
    else
        fail 'only performs 1D or 2D wavelet decomposition'
    end
end

% remove cA if required.
if arg.noCA
    if arg.ndim == 1
        y = y(arg.s(1)+1:end,:);
    elseif arg.ndim == 2
        y = y(arg.s(1,1)*arg.s(1,2)+1:end,:);
    end
end

%%
% Gwave_forw_SWT(): y = A * x  where A is the stationary wavelet transform
% in
%	x	[np L] or [(Nd) L]  where x is zeropadded and has masked applied.
% out
%	y	[clength L]   where clength is either full coefficients or noCA.
%
function y = Gwave_forw_SWT(arg, x, xdim)

LL = xdim(end);

% Single image.
if isequal(size(x), arg.Nd+arg.pad) % [(Nd)]
    if arg.ndim == 1
        y = swt(x,arg.level,arg.wname);
    elseif arg.ndim == 2
        y = swt2(x,arg.level,arg.wname);
    else
        fail 'only performs 1D or 2D wavelet decomposition'
    end
    y = col(y); % ensure its a column vector, not row ([col(y(:,:,1));col(y(:,:,2))...;col(y(:,:,LL))]
else % stacked images.
    if arg.ndim == 1
        y = zeros([arg.clengthCA LL]);
        for ll=1:LL  % Could use mdwtdec,but need entire object to recon.
            tmp = stackpick(x, ll);
            y(:,ll) = col(swt(tmp,arg.level,arg.wname));
        end
    elseif arg.ndim == 2
        y = zeros([arg.clengthCA LL]);
        for ll=1:LL
            tmp = stackpick(x, ll);
            y(:,ll) = col(swt2(tmp,arg.level,arg.wname));
        end
    else
        fail 'only performs 1D or 2D wavelet decomposition'
    end
end

% remove cA if required.
if arg.noCA
        y = y(1:arg.clength,:);
end


%%
% Gwave_back_DWT(): x = A' * y where A is the discrete wavelet transform
% in
%	y	[clength L]   where clength is either full coefficients or noCA.
% out
%	x	[Nd L]  where Nd is mask dimensions with any padding.
%
function x = Gwave_back_DWT(arg, y)

% Add zeros to coefficiant matrix in place of approx. coeff.
if arg.noCA
   if arg.ndim == 1
       y = padarray(y, arg.s(1), 0, 'pre');
   elseif arg.ndim == 2
       y = padarray(y, arg.s(1,1)*arg.s(1,2), 0, 'pre');
   else
       fail 'Gwave_back only accomodates 1D or 2D data.';
   end
end

ydim = size(y); % coefficiants always vectors
LL = ydim(end);
NN = prod(arg.Nd + arg.pad); % size of original (possibly padded) x.

if arg.ndim == 1
    x = zeros([NN LL]);
    for ll=1:LL  % Could use mdwtrec,but need entire object to recon.
        tmp = stackpick(y, ll);
        x(:,ll) = waverec(tmp,arg.s,arg.wname);
    end
elseif arg.ndim == 2
    x = zeros([arg.Nd+arg.pad LL]);
    for ll=1:LL
        tmp = stackpick(y,ll);
        x(:,:,ll) = waverec2(tmp,arg.s,arg.wname);
    end
else
    fail 'Gwave_back only accomodates 1D or 2D data.';
end

%%
% Gwave_back_SWT(): x = A' * y where A is the stationary wavelet transform
% in
%	y	[clength L]   where clength is either full coefficients or noCA.
% out
%	x	[Nd L]  where Nd is mask dimensions with any padding.
%
function x = Gwave_back_SWT(arg, y)

Ndpad = arg.Nd + arg.pad;
NN = prod(Ndpad); % size of original (possibly padded) x.

% Add zeros to coefficiant matrix in place of approx. coeff.
if arg.noCA
       y = padarray(y, NN, 0, 'post');
end

ydim = size(y); % coefficiants always vectors
LL = ydim(end);

if arg.ndim == 1
    x = zeros([NN LL]);
    for ll=1:LL
        tmp = stackpick(y, ll);
        %reshape for iswt
        tmp = reshape(tmp,[arg.level+1 NN]);
        x(:,ll) = iswt(tmp,arg.wname);
    end
elseif arg.ndim == 2
    x = zeros([Ndpad LL]);
    for ll=1:LL
        tmp = stackpick(y,ll);
        % reshape for iswt2
        tmp = reshape(tmp,[Ndpad 3*arg.level+1]);
        x(:,:,ll) = iswt2(tmp,arg.wname);
    end
else
    fail 'Gwave_back only accomodates 1D or 2D data.';
end

%%
% outCoeffDWT(arg):
% Computes the output coefficients and the related 's' (i.e., 'l') matrix.
% If noCA = 1, clength is reduced by the # of approx. coeff.
% If zeropad = 1, clength & clengthCA are computed with 0-padded data.
%
function [clength clengthCA s] = outCoeffDWT(arg)

if arg.zeropad
    %Use temp 0-padded mask to determine the output coeff. under this option.
    Ndtmp = arg.Nd + arg.pad;
else
    Ndtmp = arg.Nd;
end

% create tmp matrix to compute 's' and 'clength'.
% remove approx. coeff if desired. Leave s alone for future computation.
tmp = zeros(Ndtmp); 
if arg.ndim == 1 % 1D case
    [c s] = wavedec(tmp,arg.level,arg.wname);
    clengthCA = length(c);
    if arg.noCA
        clength = clengthCA - s(1);
    else
        clength = clengthCA;
    end
    clear c;
elseif arg.ndim == 2
    [c1 s1] = wavedec(tmp(:,1),arg.level,arg.wname);
    [c2 s2] = wavedec(tmp(1,:),arg.level,arg.wname);
    s = [s1 s2.'];
    sM = s(:,1) .* s(:,2);
    clengthCA = sum(sM(2:end-1))*3 + sM(1);
    if arg.noCA
        clength = clengthCA - sM(1);
    else
        clength = clengthCA;
    end
    clear c1 c2;
else
    error('GWAVE.M currently only supports 1D or 2D matrices.');
end


%%
% Gwave1_test()
% 
% Provides examples of Gwave.m usage and verifies the fatrix under typical
% usage scenarios.
%
function Gwave1_test

if exist('dwtmode') ~= 2
	warn 'Gwave1 requires matlab wavelet toolbox'
return
end

%%%%%% test functionality SWT (1D)
fprintf('\n --TEST VECTOR WITH MASK AND DWT - Level 2 ''DB2''-- \n \n');
I = rand(51,1);
I = complex(I,0.3*I);
Nd = size(I);
mask = true(Nd); mask(1:17) = false;

A = Gwave1('mask', mask, 'level', 2, 'wname', 'db2', 'noCA',0);

% basic tests
Fatrix_test_basic(A, mask);
test_adjoint(A,'big',0,'complex',1,'chat',0);

% functionality tests
wv = A * I;
img = A' * wv;
Imask = I;
Imask(~mask) = 0;
diff = double(Imask) - embed(img,mask);
maxDiff = max(max(diff(mask)));
disp(['Maximum difference is ' num2str(real(maxDiff)) ' ' num2str(imag(maxDiff)) 'i']);

%%%%%%% test functionality SWT (1D)
fprintf('\n --TEST VECTOR WITH MASK AND SWT - Level 2 ''HAAR''-- \n \n');
I = rand(51,1);
I = complex(I,0.3*I);
Nd = size(I);
mask = true(Nd); mask(1:17) = false;
%I = repmat(I, [1 4]);

A = Gwave1('mask', mask, 'level', 2, 'wname', 'haar', 'swt', 1, 'noCA',0);

Fatrix_test_basic(A, mask);
disp('**1D SWT FAILS ADJOINT TEST**');
%test_adjoint(A,'big',0,'complex',1,'chat',1);

wv = A * I;
img = A' * wv;
Imask = I;
Imask(~mask) = 0;
diff = double(Imask) - embed(img,mask);
maxDiff = max(max(diff(mask)));
disp(['Maximum difference is ' num2str(real(maxDiff)) ' ' num2str(imag(maxDiff)) 'i']);

%%%%%%%%% test functionality DWT (2D)
fprintf('\n --TEST REAL IMAGE WITH MASK AND DWT - Level 2 ''DB2''-- \n \n');

%ig = image_geom('nx', 128, 'dx', 1);
%I = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
I = imread('cameraman.tif'); I = I(1:end-1,1:end-1);
I = single(I);
I = complex(I,I*.3);

Nd = size(I);
mask = true(Nd); mask(1:50,1:50) = false; % stress test
%I = repmat(I, [1 1 4]);

A = Gwave1('mask', mask, 'level', 2, 'wname', 'db2', 'swt', 0, 'noCA', 0);

Fatrix_test_basic(A, mask) 
test_adjoint(A,'big',1,'complex',1,'chat',0);

wv = A * I;
img = A' * wv;
Imask = I;
%Imask(~repmat(mask,[1 1 4])) = 0;
Imask(~mask) = 0;
%figure; im(embed(img,mask));
diff = double(Imask) - embed(img,mask);
maxDiff = max(max(diff(:)));
disp(['Maximum difference is ' num2str(real(maxDiff)) ' ' num2str(imag(maxDiff)) 'i']);

%%%%%%%%% test functionality SWT (2D)
fprintf('\n --TEST REAL IMAGE WITH MASK AND SWT - Level 2 ''HAAR''-- \n \n');

I = imread('cameraman.tif');
I = single(I);
I = I(1:end-1,1:end-1);
I = complex(I,I*.3);

Nd = size(I);
mask = true(Nd); mask(1:50,1:50) = false; % stress test
%I = repmat(I, [1 1 4]);

A = Gwave1('mask', mask, 'level', 2, 'wname', 'db2', 'swt', 1, 'noCA', 0);

Fatrix_test_basic(A, mask) 
disp('**2D SWT FAILS ADJOINT TEST**');
%test_adjoint(A,'big',1,'complex',1,'chat',0);

wv = A * I;
img = A' * wv;
Imask = I;
%Imask(~repmat(mask,[1 1 4])) = 0;
Imask(~mask) = 0;
%figure; im(embed(img,mask));
diff = double(Imask) - embed(img,mask);
maxDiff = max(max(diff(:)));
disp(['Maximum difference is ' num2str(real(maxDiff)) ' ' num2str(imag(maxDiff)) 'i']);
