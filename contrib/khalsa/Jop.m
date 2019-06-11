 function ob = Jop(varargin)
%function ob = Jop(varargin)
%J = Jop({'J', 4, 'L', 4, 'kspace', kspace, 'fov', fov, ...
%     'kn.ktype', 'kaiser', 'kn.kb_alf', kb_alf, 'kn.kb_m', kb_m});
%
% Creates fatrix J for use in operations of the form y = J * x.  Designed
% for use in iterative algorithms for estimating density compensation
% weights for conjugate phase MRI reconstruction.
%
% J is an MxM symmetric banded matrix, with entries corresponding to a
% gridding kernel, C, convolved with itself, evaluated at values
% corresponding to the difference between kspace sample locations, i.e.
%       J(i,l) = C(*)C(v_i - v_l)
% where C is the gridding kernel (usually a Kaiser-Bessel), (*) represents 
% the convolution operator, and kspace sample locations are {v_m}, 
% m = 1,..., M.
%
% in: 
%   kspace          kspace sample locations
%   fov             Field of View
%
%  optional inputs:
%   J               neighborhood for interpolation kernel
%   L               samples per integer
%   kn.*            various information to define kernel
%                   e.g. kn.ktype = 'kaiser', kn.kb_alf = 10, kn.kb_m = 0
% out:
%   J   [nd nd]     Fatrix operator
%
% K. Khalsa, Mar. 2006


% default assignments:

arg.J = 4; 
arg.L = 4;
arg.kspace = [];
arg.fov = [];           
arg.del = [];           

arg.kn.ktype = 'kaiser';
arg.kn.kernel = {};     %*** find out what to do with this, see Gn.arg.st
arg.kn.kb_alf = [];     %
arg.kn.kb_m = [];       %

% pair input arguments
if iscell(varargin)
    varargs = varargin{1};
    arg = vararg_pair(arg, varargs);
else
    printf('input needs to be in a cell array.');
    help Jop;
    return;
end



if isempty(arg.kspace)
    error('kspace sample locations are a required input argument');
else
    arg.M = size(arg.kspace, 1);
    arg.dim = [arg.M arg.M];
end

arg.del = 1 ./ (arg.fov * arg.L);

kappa = linspace(-arg.J/2, arg.J/2, arg.J * arg.L + 1)';
%kappa = [-arg.J/2 : 1/arg.L : arg.J/2]';   %equivalently


if (arg.kn.ktype == 'kaiser')

    if xor(isempty(arg.kn.kb_alf), isempty(arg.kn.kb_m))
        printf('For Kaiser-Bessel kernel, both m and alpha are required, or neither');
        return;
    elseif (isempty(arg.kn.kb_alf) && isempty(arg.kn.kb_m))
        % use parameters from Jackson's 1991 gridding paper
        arg.kn.kb_m = 0;
        ww = [1.5 2.0 2.5 3.0 3.5 4.0];
        jack_alf = [6.6875 9.1375 11.5250 13.9086 16.2734 18.5547];
        arg.kn.kb_alf = interp1(ww, jack_alf, arg.J, 'linear', 'extrap');
        
    end
    % keyboard
    arg.kn.C = kaiser_bessel(kappa, arg.J, arg.kn.kb_alf, ...
        arg.kn.kb_m);
    if size(arg.kspace, 2) == 2
        arg.kn.C = arg.kn.C * arg.kn.C';    % 2D kernel
    end
    arg.kn.c0 = kaiser_bessel_ft(0, arg.J, arg.kn.kb_alf, arg.kn.kb_m);

elseif (arg.kn.ktype == 'function_handle')
    % figure out what to do here
    % see Gn.arg.st.kernel = {[1x1] function_handle}... ??
else
    error('kn.ktype must be kaiser or function_handle');
end


% Build Fatrix object
ob = Fatrix(arg.dim, arg, 'forw', @J_forw, 'back', @J_forw, ...
    'caller', mfilename);
% transpose multiplication = same as fwd multiplication b/c we
% stipulate that all kernels must be real and symmetric

%-------------------------- 
% multiplication
function y = J_forw(arg, x)

if (size(x,1) ~= arg.dim(2))
    error('dimension mismatch in matrix vector multiplication');
end

% 1D case
if size(arg.kspace, 2) == 1
    k0 = abs(min(arg.kspace / arg.del) - arg.J);
    k = round(arg.kspace / arg.del + k0);
    N = max(k) + arg.J;
    xtmp = full(sum(sparse(1:arg.M, k, x, arg.M, N), 1))';

    CC = conv(arg.kn.C, arg.kn.C);
    Jxtmp = conv(xtmp, CC);
    Jxtmp = arg.del * Jxtmp;
    Jxtmp = Jxtmp / max(CC);    % adjust scaling?
    y = Jxtmp(k + arg.J * arg.L);


    % 2D case
elseif size(arg.kspace, 2) == 2

    k01 = abs(min(arg.kspace(:,1) / arg.del(1)) - arg.J);
    k02 = abs(min(arg.kspace(:,2) / arg.del(2)) - arg.J);
    k0 = [k01, k02];

    k1 = round(arg.kspace(:,1) / arg.del(1) + k0(1));
    k2 = round(arg.kspace(:,2) / arg.del(2) + k0(2));
    N1 = max(k1); N2 = max(k2);

    xtmp = full(sparse(k1, k2, x, N1, N2));

    CC = conv2(arg.kn.C, arg.kn.C);
    Jxtmp = conv2(xtmp, CC, 'same');
    ind = sub2ind(size(Jxtmp), k1,k2);

    y = Jxtmp(ind);
    %  y = y / max(CC(:));  % adjust scaling?
else
    error('only 1D and 2D currently supported')
end
y = y / ((arg.kn.c0)^2); % adjust scaling??
% keyboard


% %plot to see if it's working right 2D
%     figure(11), clf
%     subplot(211)
%     plot3(arg.kspace(:,1), arg.kspace(:,2), x, 'x');
%     title('ks locs before regridding'), xlabel('kx'), ylabel('ky')
%     subplot(212)
%     plot3(k1, k2, x, 'x'), title('ks locs after regridding'), xlabel('kapx')
% 
% 
%   figure(12), clf
%   subplot(211), imagesc(xtmp), axis square, colorbar, title('xtmp')
%   subplot(212), imagesc(Jx), axis square, colorbar, title('Jx')
% keyboard

%%plot to see if it's working right 1D
%     figure(11), clf
%     subplot(211)
%     stem(arg.kspace, x), title('wi before regridding'), xlabel('kx')
%     subplot(212)
%     stem(1:N, xtmp), title('wi after regridding'), xlabel('kap')

%   figure(12), clf
%     subplot(211), stem(1:length(Jxkap), Jxkap), title('Jxkap')
%     subplot(212), stem(arg.kspace, y), title('y = J * x'), xlabel('kx')
% keyboard



