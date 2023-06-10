function [wmap, wconv] = call_fmap(yik, etime, func_name, varargin)
% This function is a wrapper for fmap_est_qm.m and fmap_est_pcg_ls.m to
% make them callable in a similar fashion as mri_field_map_reg.m

%| inputs:
%|	yik	[(N) nset ncoil]	complex images at different echo times
%|				(these should be undistorted)
%|				(N) can be (nx,ny) or (nx,ny,nz) etc.
%|	etime	[nset]		vector of echo times [usually in seconds]
%|  func_name   string specifying which function to call
%|
%| fmap_est_pcg_ls options
%|	stepper {'qs',# its}	monotonic line search parameters (default: {})
%|	niter			# of iterations (def: 30)
%|	maskR	[(np)]	logical reconstruction mask (required!)
%|	order			order of the finite diff matrix (def: 2)
%|	l2b             regularization parameter (2^) (def: -6)
%|	gammaType		CG direction: PR = Polak-Ribiere (default) or FR = Fletcher-Reeves
%|	precon			Preconditioner: 'diag', 'chol', 'ichol' (def: 'ichol')
%|	reset			# of iterations before resetting direction (def: inf)
%|  df              delta f value in water-fat imaging (def: 0)
%|  relamp          relative amplitude in multipeak water-fat  (def: 1)
%|  tol             tolerance for ichol (def: 1e-3)
%|
%| fmap_est_qm options
%|	niter			# of iterations (default: 30)
%|	maskR	[(np)]	logical support mask (required)
%|	order			order of the finite diff reg. (def: 2)
%|	l2b             regularization parameter (2^) (def: -6)
%|  hess            'diag' (default) or 'chol'
%|  df              delta f value in water-fat imaging (def: 0)
%|  relamp          relative amplitude in multipeak water-fat  (def: 1)
%|
%| out
%|	wmap	[(N)]	regularized estimate of field map [rad/sec]
%|	wconv	""	conventional estimate from first two scans
%|
%| Written 2021-06-03, Melissa Haskell, University of Michigan

% default optional arguments
arg.quiet = true;
arg.winit = [];
arg.sens = [];
arg.maskR = [];
arg = vararg_pair(arg, varargin, 'allow_new', 1);


% determine data size & reshape as needed
if numel(size(yik)) == 3
    % 2d single coil
    [nlin, ncol, necho] = size(yik); size_v = [nlin, ncol];
    y = reshape(yik,[nlin*ncol,1,necho]);
    arg.sens = ones(nlin*ncol,1);
elseif numel(size(yik)) == 4
    if isempty(arg.sens)
        % 3d single coil
        [nlin, ncol, nsli, necho] = size(yik); size_v = [nlin, ncol, nsli];
        y = reshape(yik,[nlin*ncol*nsli,1,necho]);
        arg.sens = ones(nlin*ncol*nsli,1);
    else
        % 2d multicoil
        [nlin, ncol, necho, ncoil] = size(yik); size_v = [nlin, ncol];
        y = reshape(permute(yik,[1 2 4 3]),[nlin*ncol,ncoil,necho]);
    end
else
    % 3d multicoil
    [nlin, ncol, nsli, necho, ncoil] = size(yik); size_v = [nlin, ncol, nsli];
    y = reshape(permute(yik,[1 2 3 5 4]),[nlin*ncol*nsli,ncoil,necho]);
end

% initialize mask
if isempty(arg.maskR)
    if ~exist('nsli','var')
        arg.maskR = true(nlin,ncol);
    else
        arg.maskR = true(nlin,ncol,sli);
    end
end

% conventional field map estimate from first two scans
wconv = angle(stackpick(yik,2) .* conj(stackpick(yik,1))) ...
	/ (etime(2) - etime(1));

if isempty(arg.winit), arg.winit = wconv; end
w = arg.winit(:);


% prep optional arguments to be passed into fmap function
arg_fmap = arg;
arg_fmap = rmfield(arg_fmap,'sens');   % remove sensitivity maps, since it's now an input
arg_fmap = rmfield(arg_fmap,'winit');  % remove image initialization, since it's now an input
arg_fmap = rmfield(arg_fmap,'quiet');  % remove option only used in this function
varargin_fmap = arg2varargin(arg_fmap);

% call function
command = strcat("wmap_structure = ",func_name,"(w, y, etime, arg.sens, varargin_fmap{:});");
if arg.quiet
    evalc(command);
else
    eval(command)
end
wmap = reshape(wmap_structure.ws(:,end),size_v);
end



function [varargin] = arg2varargin(arg)
% converts a variable input arguments structure, arg,  into a cell of
% fieldnames and inputs to be passed along to another function as a
% varargin cell, i.e. [out] = myfunc(in1, in2, varargin{:})
%
% IMPORTANT NOTE! when passing the output of this function, varargin, to
% the next fucntion, the "{:}" at the end is necessary
%
% Inputs:
%   arg     struct of variable inputs
%
% Outputs:
%   vargin  cell of fieldname strings and inputs
%
% Written 2021-06-03, Melissa Haskell, University of Michigan

varargin = [fieldnames(arg), struct2cell(arg)]';


end
