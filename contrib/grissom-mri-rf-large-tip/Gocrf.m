 function [ob,m,mz] = Gocrf(k,fov,N,mdomswitch,dt,varargin)
%function [ob,m,mz] = Gocrf(k,fov,N,m0,mdomswitch,dt,varargin)
%
% Construct Gocrf object for fast optimal control large-tip-angle 
% multidimensional (and parallel) MR RF pulse design. Used in 
% conjunction with an optimization routine, this object can be 
% used to design small pulse updates that, when summed with an 
% underlying pulse, result in an accurate excitation.
% The vectors m and mz are the results of initial Bloch sim at 
% object creation.
%
% The object is created by calling:
%    Gocrf = Gocrf( ... );
% and it is then used by commands like:
%    y = Gocrf * x;
% which will evaluate perturbations to the magnetization
% produced by pulse perturbations x. It also passes x with 
% unit weighting, so that if G is the sub-object relating 
% x to its perturbations, then y = [G;I]*x
%
% in: 
%    k            [Nt,d]   excitation k-space trajectory
%                          (= -gam int_t^T g(t') dt')
%                          (inverse fov units)
%    fov          [d,1]    field of view in each dimension
%    N            [d,1]    size of design grid in each
%                          dimension
%    mdomswitch   [1]      switch to evaluate magnetization
%                          or spinors (spinors not tested yet)
%    dt           [1]      RF sampling period 
%
% options: 
%    m0           [3,1]    initial magnetization state
%                          (e.g., if equilibrium, [0 0 1] (default))
%    'sens'       [prod(N),Ncoils]   transmit (B1+) sensitivities
%    'd'          [[N]]    desired excitation pattern, in 
%                          arbitrary flip angle units. this is used
%                          to determine sub-space over which to
%                          simulate
%    'indmask'    [[N]]    mask to choose sub-space points from
%    'baseB1'     [Ncoils*Nt]    baseline RF pulse(s), used for
%                                initial Bloch simulation during 
%                                object creation
%    'g'          [Nt,d]   gradient waveforms, for Bloch simulation
%                          (gauss/fov units)
%    'Nsvd'       [1]      # of baseline excitation expansion terms
%                          (default is 4)
%    'a0B','b0B'  [Nt,Nsvd]      baseline temporal basis function
%                                matrices
%    'a0Ct','b0Ct' [prod(N),Nsvd]   baseline spatial coefficient
%                                   matrices
%    'a0','b0'    [prod(N),1]    final baseline spinors
%    'nufft'      cell     nufft arguments (passed to Gmri_SENSE)
%    'RFbasis'    [Nt,Nbasis]    RF pulse basis functions for
%                                parameterized pulse design
% 
% options for field-corrected pulse design (untested in this object!):
%    'ti'         [Nt,1]   sample times (0:T in seconds)
%    'we'         [[N]]    field_map (Hz)
%    'L'          [1]      # of approximation terms (default is 4)
% 
% After building this object, the user can Bloch-simulate
% a new pulse and update the object using:
%     [ob,m,mz] = feval(ob.arg.update, ob, baseB1);
%
% or a Bloch simulation (without object update) can be performed
% using:
%     [m,mz,a0,b0] = feval(ob.arg.blochsim,baseB1,arg,pos,sens,we);
%
% OR an update can be made (without Bloch sim) using:
%     ob = feval(ob.arg.new_B_Ct,ob,a0B,a0Ct,b0B,b0Ct,a0,b0)
% 
% Copyright 2008-10-8, Will Grissom, Stanford University

if nargin < 6, help(mfilename), error(mfilename), end

%
% default object
%
% Gmri objects 
arg.a0G = [];arg.b0G = [];
% svd matrices
arg.a0B = [];arg.a0Ct = [];arg.b0B = [];arg.b0Ct = [];
% underlying final alpha0, beta0
arg.a0 = [];arg.b0 = [];
% initial mag (global)
arg.bsq = 0; % 
arg.exc = 0; % just return mxy
arg.rfb = 0; % no rf basis
arg.m0 = [0 0 1];
arg.mdomswitch = mdomswitch;
arg.we = [];
arg.nufft_args = {N,6*ones(size(N)),2*N,N/2,'table',2^10,'minmax:kb'};    
arg.new_B_Ct = @Gocrf_new_B_Ct; 
arg.blochsim = @Gocrf_blochsim;
arg.subindcalc = @Gocrf_subindcalc;
arg.update = @Gocrf_update;
arg.subind = []; % set of points to Bloch sim and do SVD on
arg.nsubbins = 10; % default # of histogram bins to choose sub-points from
arg.simpos = []; % spatial locations for Bloch sim
arg.sens = []; % need sens for Bloch sims
arg.d = []; % desired pattern (for subind calc)
arg.indmask = []; % mask (for subind calc)
arg.simmask = []; % mask for Bloch simulation
arg.dt = dt; % sampling period
arg.ti = []; % time points (for field-map corrected design)
arg.fov = fov; 
arg.k = k; % excitation k-space traj 
arg.N = N; % spatial dimensions
arg.L = 4;
arg.Nsvd = 4;
arg.subindlen = 15*arg.Nsvd; % default # of sub-points for svd sim
% field map bases and coeffs
arg.B = [];arg.Ct = [];
arg.baseB1 = []; % base B1 pulses (for Bloch sim)
arg.g = []; % gradient waveforms (for Bloch sim)
arg.fmod = []; % modulation freq (for eg dual band designs)
arg.aux1 = 1; % vector to multiply a or mxy by
arg.aux2 = 1; % vector to multiply b or b^2 or mz by
arg.RFbasis = [];
ob.Nbasis = 0;
% Bloch-simulated patterns (will be empty if basis and coeff
% matrices are provided)
m = [];mz = [];

gambar = 4257;             % gamma/2pi in Hz/g
gam = gambar*2*pi;         % gamma in radians/g

arg = vararg_pair(arg, varargin, 'subs', ...
	{'basis', 'basis_args'; 'nufft', 'nufft_args'});

% if not provided, get the subindices for Bloch sim
if isempty(arg.subind) & ~isempty(arg.d) & ~isempty(arg.indmask)
  arg = Gocrf_subindcalc(arg,arg.d,arg.indmask);
elseif isempty(arg.subind)
  disp ['No sub indices or means to calculate them provided. update() ' ...
        'will not work'];
end

ktmp = -[arg.k - repmat(arg.k(1,:)/2,[size(arg.k,1) 1])];
if isempty(arg.we) || sum(abs(arg.we(:))) == 0  % no field map
  
	if ~isempty(arg.sens)
		arg.a0G = Gmri_SENSE(ktmp,true(N),'fov',fov, ...
			'basis',{'dirac'},'nufft',arg.nufft_args, ...
			'exact',0,'sens',arg.sens*(-1i*gam*arg.dt/2))';
		arg.b0G = Gmri_SENSE(ktmp,true(N),'fov',fov, ...
			'basis',{'dirac'},'nufft',arg.nufft_args, ...
			'exact',0,'sens',arg.sens*(1i*gam*arg.dt/2))';
	else
%{
		arg.a0G = Gmri_SENSE(ktmp,true(N),'fov',fov,'basis',{'dirac'},'nufft',arg.nufft_args, ...
                         'exact',0,'sens',ones(prod(N),1)*(-1i*gam*arg.dt/2))';
		arg.b0G = Gmri_SENSE(ktmp,true(N),'fov',fov,'basis',{'dirac'},'nufft',arg.nufft_args, ...
                         'exact',0,'sens',ones(prod(N),1)*(1i*gam*arg.dt/2))';
%}
		% jf version:
		arg.a0G = (-1i*gam*arg.dt/2)' * ...
			Gmri(ktmp, true(N), 'fov', fov, 'basis', {'dirac'}, ...
				'nufft',arg.nufft_args)';
		arg.b0G = (1i*gam*arg.dt/2)' * ...
			Gmri(ktmp, true(N), 'fov', fov, 'basis', {'dirac'}, ...
				'nufft', arg.nufft_args)';

	end

else % non-zero field map
  
	if ~isempty(arg.sens)
		arg.a0G = Gmri_SENSE(ktmp,true(N),'fov',fov,...
			'basis',{'dirac'},'nufft',arg.nufft_args, ...
			'exact',0,'sens',arg.sens*(-1i*gam*arg.dt/2), ...
			'ti',-arg.ti+arg.ti(end)/2,...
			'zmap',1i*2*pi/2*arg.we,'L',arg.L)';
		arg.b0G = Gmri_SENSE(ktmp,true(N),'fov',fov,...
			'basis',{'dirac'},'nufft',arg.nufft_args, ...
			'exact',0,'sens',arg.sens*(1i*gam*arg.dt/2), ...
			'ti',-arg.ti+arg.ti(end)/2, ...
			'zmap',1i*2*pi/2*arg.we,'L',arg.L)';
	else
%{
		arg.a0G = Gmri_SENSE(ktmp,true(N),'fov',fov,...
			'basis',{'dirac'},'nufft',arg.nufft_args, ...
			'exact',0,'sens',ones(prod(N),1)*(-1i*gam*arg.dt/2), ...
			'ti',-arg.ti+arg.ti(end)/2, ...
			'zmap',1i*2*pi/2*arg.we,'L',arg.L)';
		arg.b0G = Gmri_SENSE(ktmp,true(N),'fov',fov, ...
			'basis',{'dirac'},'nufft',arg.nufft_args, ...
			'exact',0,'sens',ones(prod(N),1)*(1i*gam*arg.dt/2), ...
 			'ti',-arg.ti+arg.ti(end)/2, ...
			'zmap',1i*2*pi/2*arg.we,'L',arg.L)';
%}
		% jf version:
		arg.a0G =(-1i*gam*arg.dt/2)' * ...
			Gmri(ktmp,true(N),'fov',fov,...
				'basis',{'dirac'},'nufft',arg.nufft_args, ...
				'ti',-arg.ti+arg.ti(end)/2, ...
				'zmap',1i*2*pi/2*arg.we,'L',arg.L)';
		arg.b0G = (1i*gam*arg.dt/2)' * ...
			Gmri(ktmp,true(N),'fov',fov, ...
				'basis',{'dirac'},'nufft',arg.nufft_args, ...
 				'ti',-arg.ti+arg.ti(end)/2, ...
				'zmap',1i*2*pi/2*arg.we,'L',arg.L)';
	end

	arg.B = arg.a0G.arg.B;
	arg.Ct = arg.a0G.arg.Ct;
end

tmp = arg.a0G;
[dim1,dim2] = size(tmp);
if arg.exc
  arg.dim = [dim1 + 2*dim2, 2*dim2];
else
  arg.dim = [2*dim1+2*dim2, 2*dim2];
end

if ~isempty(arg.RFbasis)

  arg.rfb = 1;  
  arg.Nbasis = size(arg.RFbasis,2);
  
  % assuming that there is only one coil when we parameterize the
  % pulse, adjust object dims
  if arg.exc
    arg.dim = [dim1 + 2*arg.Nbasis, 2*arg.Nbasis];
  else
    arg.dim = [2*dim1 + 2*arg.Nbasis, 2*arg.Nbasis];
  end
  
  % allocate space for basis perturbations
  arg.Gat = zeros(dim1,arg.Nbasis);
  arg.Gbt = zeros(dim1,arg.Nbasis);
  
end
               
arg.is.empty = false;

% now build object
%ob = Fatrix(arg.dim, arg, 'caller', mfilename, ...
%	'forw', @Gocrf_forw, 'back', @Gocrf_back);
ob = fatrix2('odim', arg.dim(1), 'idim', arg.dim(2), ...
	'arg', arg, ...
	'forw', @Gocrf_forw, 'back', @Gocrf_back); % jf

% if SVD matrices not supplied, do the necessary Bloch sims
if isempty(ob.arg.a0B)
  % get Bloch simulation grid
  ndgridarg = ['-fov(1)/2:fov(1)/N(1):fov(1)/2-fov(1)/N(1)'];
  posarg = ['x1'];
  for ii = 2:length(N)
    ndgridarg = [ndgridarg sprintf(',-fov(%d)/2:fov(%d)/N(%d):fov(%d)/2-fov(%d)/N(%d)',ii,ii,ii,ii,ii,ii)];
    posarg = [posarg sprintf(',x%d',ii)];
  end
  eval(sprintf('[%s] = ndgrid(%s);',posarg,ndgridarg));
  ob.arg.simpos = [x1(:)];
  for ii = 2:length(N)
    ob.arg.simpos = [ob.arg.simpos eval(sprintf('x%d(:)',ii))];
  end
  % do the Bloch sims and expansion calc
  [ob,m,mz] = Gocrf_update(ob,ob.arg.baseB1);
  % we don't need base B1 anymore
  ob.arg.baseB1 = [];
end
  

%
% Gocrf_new_B_Ct()
% new values of (svd) bases and coefficients
%
function G = Gocrf_new_B_Ct(G, a0B, a0Ct, b0B, b0Ct, a0, b0)

% apply frequency modulation to basis matrices
if ~isempty(G.arg.fmod)
  t = [0:G.arg.dt:(size(G.arg.k,1)-1)*G.arg.dt] - size(G.arg.k,1)*G.arg.dt/2;
  Amod = repmat(exp(1i*2*pi*G.arg.fmod*t(:)),[1 G.arg.Nsvd]);
  a0B = a0B .* Amod;
  b0B = b0B .* Amod;
end

if ~isempty(G.arg.we) & sum(abs(G.arg.we(:))) ~= 0
  % multiply out new expansion with field map expansion
  newa0B = [];newa0Ct = [];
  newb0B = [];newb0Ct = [];
  for ll = 1:G.arg.L
    newa0B = [newa0B spdiag(G.arg.B(:,ll))*a0B];
    newa0Ct = [newa0Ct spdiag(G.arg.Ct(:,ll))*a0Ct];
    newb0B = [newb0B spdiag(G.arg.B(:,ll))*b0B];
    newb0Ct = [newb0Ct spdiag(G.arg.Ct(:,ll))*b0Ct];
  end
  % update G with provided B and C.' matrices
  G.arg.a0G = feval(G.arg.a0G.arg.new_B_Ct,G.arg.a0G,newa0B,newa0Ct);
  G.arg.b0G = feval(G.arg.b0G.arg.new_B_Ct,G.arg.b0G,newb0B,newb0Ct);
else
  % set the expansion directly
  G.arg.a0G = feval(G.arg.a0G.arg.new_B_Ct,G.arg.a0G,a0B,a0Ct);
  G.arg.b0G = feval(G.arg.b0G.arg.new_B_Ct,G.arg.b0G,b0B,b0Ct);
end

G.arg.a0 = a0;
G.arg.b0 = b0;


% 
% Gocrf_forw(): y = G * x
%
function v = Gocrf_forw(arg, u)

% assume the RF vector is [real(db1); imag(db1)]
db1 = u(1:length(u)/2) + 1i*u(length(u)/2+1:end);

if ~arg.rfb
  if any(db1)
    at = arg.b0G * conj(db1);
    bt = arg.a0G * conj(db1);
  else
    at = zeros(size(arg.b0G,1),1);
    bt = zeros(size(arg.a0G,1),1);
  end
else % pulse is parameterized
  if any(db1)
    at = arg.Gat * conj(db1);
    bt = arg.Gbt * conj(db1);
  else
    at = zeros(size(arg.Gat,1),1);
    bt = zeros(size(arg.Gbt,1),1);
  end
end

if arg.mdomswitch
  % dm terms, equilibrium
  dm = 0;dmz = 0;
  if arg.m0(3) ~= 0
    % dm terms, equilibrium
    dm = -2*arg.m0(3)*conj(arg.a0(:).*bt + at.*arg.b0(:));
    % dmz, equilibrium
    dmz = 2*arg.m0(3)*real(arg.a0(:).*conj(at)) - 2*arg.m0(3)*real(arg.b0(:).*conj(bt));
  end
  if any(arg.m0(1:2) ~= 0)
    % dm terms, non-equilibrium
    dm = dm - (arg.m0(1)-1i*arg.m0(2))*conj(2*arg.b0(:).*bt);
    dm = dm + (arg.m0(1)+1i*arg.m0(2))*conj(2*arg.a0(:).*at);
    % dmz, non-equilibrium
    dmz = dmz + 2*real((arg.m0(1)+1i*arg.m0(2))*(conj(arg.a0(:)).*bt + conj(at).*arg.b0(:)));
  end
  if ~arg.exc
    v = [arg.aux1.*dm;arg.aux2.*dmz;u];
  else
    v = [arg.aux1.*dm;u];
  end
else
  if ~arg.bsq
    v = [arg.aux1.*at;arg.aux2.*bt;u];  
  else
    v = [2*arg.aux1.*arg.a0(:).*at;2*arg.aux2.*arg.b0(:).*bt;u];
  end
end


% 
% Gocrf_back(): x = G' * y
%
function v = Gocrf_back(arg,u)

v = zeros(size(arg.a0G,2),1);
Ns = size(arg.a0G,1);
tmp1 = 0; tmp2 = 0;

% tmp1 contains terms that will multiply by system matrix for at
if arg.mdomswitch
  if arg.m0(3) ~= 0
    if ~arg.exc
      tmp1 = arg.m0(3)*2*conj(-arg.b0(:).*u(1:Ns)+conj(arg.a0(:)).*u(Ns+1:2*Ns));
    else
      tmp1 = arg.m0(3)*2*conj(-arg.b0(:).*u(1:Ns));
    end
  end
  if any(arg.m0(1:2) ~= 0)
    % additional alphat terms from Mxy, if non-equilibrium ic
    tmp1 = tmp1 + 2*conj((arg.m0(1)-1i*arg.m0(2))*arg.a0(:).*u(1:Ns)); 
    % additional alphat terms from Mz, if non-equilibrium ic
    if ~arg.exc
      tmp1 = tmp1 + 2*(arg.m0(1)+1i*arg.m0(2))*arg.b0(:).*conj(u(Ns+1:2*Ns));
    end
  end
  % tmp2 contains terms that will multiply by system matrix for betat
  if arg.m0(3) ~= 0
    if ~arg.exc
      tmp2 = -arg.m0(3)*2*conj(arg.a0(:).*u(1:Ns)+conj(arg.b0(:)).*u(Ns+1:2*Ns));
    else
      tmp2 = -arg.m0(3)*2*conj(arg.a0(:).*u(1:Ns));
    end
  end
  if any(arg.m0(1:2) ~= 0)
    % additional betat terms from Mxy, if non-equilibrium ic
    tmp2 = tmp2 - 2*conj((arg.m0(1)+1i*arg.m0(2))*arg.b0(:).*u(1:Ns)); 
    % additional betat terms from Mz, if non-equilibrium ic
    if ~arg.exc
      tmp2 = tmp2 + 2*conj((arg.m0(1)+1i*arg.m0(2))*conj(arg.a0(:)).*u(Ns+1:2*Ns));
    end
  end
else

  if ~arg.bsq
    tmp1 = u(1:Ns);
    tmp2 = u(Ns+1:2*Ns);
  else
    tmp1 = 2*conj(arg.a0(:)).*u(1:Ns);
    tmp2 = 2*conj(arg.b0(:)).*u(Ns+1:2*Ns);
  end
end

if ~arg.rfb
  v = conj(arg.a0G' * (conj(arg.aux2).*tmp2));
  v = v + conj(arg.b0G' * (conj(arg.aux1).*tmp1));
else
  v = conj(arg.Gbt' * (conj(arg.aux2).*tmp2));
  v = v + conj(arg.Gat' * (conj(arg.aux1).*tmp1));
end

if ~arg.exc
  v = [real(v);imag(v)] + u(2*Ns+1:end);
else
  v = [real(v);imag(v)] + u(Ns+1:end);
end


% 
% Gocrf_subindcalc(): Determine subset of points to simulate for 
%                     SVD evaluation
function arg = Gocrf_subindcalc(arg,d,indmask)

% get bin centers
if min(d(:)) == max(d(:))
  fa = min(d(:))*ones(arg.nsubbins+1,1);
else
  fa = min(d(:)):(max(d(:))-min(d(:)))/arg.nsubbins:max(d(:));
end

dind = 1:prod(size(d));

% arrange indices into bins
indices{1} = dind(col(logical(indmask.*(d >= fa(1) & d <= fa(2)))));
for ii = 2:length(fa)-1
  indices{ii} = dind(col(logical(indmask.*(d > fa(ii) & d <= fa(ii+1)))));
end
  
% loop through bins again and again, pulling one from each bin
% each time around at random (use rand*size(bin)), until 
% index vector is full. This works much better than my
% previous method, which failed for filter-based desired patterns
indvec = [];
indind = 0;
rand('twister',1000);
while length(indvec) < arg.subindlen
  if ~isempty(indices{indind+1})
    binind = ceil(length(indices{indind+1})*rand(1,1));
    indvec = [indvec indices{indind+1}(binind)];
    indices{indind+1} = [indices{indind+1}(1:binind-1) ...
                        indices{indind+1}(binind+1:end)];
  end
  indind = mod(indind + 1,arg.nsubbins);
end

arg.subind = indvec;


%
% Gocrf_update(): Perform Bloch sims to update SVD expansions
% 
function [G,m,mz] = Gocrf_update(G,baseB1)  

% perform reduced Bloch sim
if isfield(G.arg, 'sens') && ~isempty(G.arg.sens) % jf

  if ~isempty(G.arg.we) & sum(abs(G.arg.we(:))) ~= 0
    [foo,fooz,alpha0,beta0] = Gocrf_blochsim(baseB1,G.arg,G.arg.simpos(G.arg.subind,:), ...
                                             G.arg.sens(G.arg.subind,:),G.arg.we(G.arg.subind));
  else
    [foo,fooz,alpha0,beta0] = Gocrf_blochsim(baseB1,G.arg,G.arg.simpos(G.arg.subind,:), ...
                                             G.arg.sens(G.arg.subind,:),G.arg.we);
  end

else

  if ~isempty(G.arg.we) & sum(abs(G.arg.we(:))) ~= 0
    [foo,fooz,alpha0,beta0] = Gocrf_blochsim(baseB1,G.arg,G.arg.simpos(G.arg.subind,:), ...
                                             ones(size(G.arg.subind)),G.arg.we(G.arg.subind));
  else
    [foo,fooz,alpha0,beta0] = Gocrf_blochsim(baseB1,G.arg,G.arg.simpos(G.arg.subind,:), ...
                                             ones(size(G.arg.subind)),G.arg.we);
  end

end

% calculate rotating frame transformation
ktmp = -[G.arg.k - repmat(G.arg.k(1,:),[size(G.arg.k,1) 1])];
Asub = exp(1i*2*pi/2*G.arg.simpos(G.arg.subind,:)*ktmp');
if ~isempty(G.arg.fmod)
  t = 0:G.arg.dt:G.arg.dt*(size(G.arg.k,1)-1);
  Asub = Asub.*exp(-1i*2*pi/2*repmat(G.arg.fmod*t,[length(G.arg.subind) 1]));
end

if ~isempty(G.arg.we) & sum(abs(G.arg.we(:))) ~= 0
  t = 0:G.arg.dt:G.arg.dt*(size(G.arg.k,1)-1);
  Asub = Asub.*exp(-1i*2*pi/2*G.arg.we(G.arg.subind).'*t);
end

% get temporal basis funcs via svd
[foo1,foo2,G.arg.a0B] = svd(conj(alpha0)./Asub,'econ');
[foo1,foo2,G.arg.b0B] = svd(conj(beta0)./Asub,'econ');

% truncate basis
G.arg.a0B = G.arg.a0B(:,1:G.arg.Nsvd);
G.arg.b0B = G.arg.b0B(:,1:G.arg.Nsvd);

% re-simulate to get coeffs at all spatial locs
if ~isempty(G.arg.sens)
  [m,mz,G.arg.a0,G.arg.b0,G.arg.a0Ct,G.arg.b0Ct] = Gocrf_blochsim(baseB1,G.arg,G.arg.simpos,G.arg.sens,G.arg.we,G.arg.simmask);
else
  [m,mz,G.arg.a0,G.arg.b0,G.arg.a0Ct,G.arg.b0Ct] = Gocrf_blochsim(baseB1,G.arg,G.arg.simpos,ones(size(G.arg.simpos,1),1),G.arg.we,G.arg.simmask);
end

m = reshape(m,G.arg.N);
mz = reshape(mz,G.arg.N);
G.arg.a0 = reshape(G.arg.a0,G.arg.N);
G.arg.b0 = reshape(G.arg.b0,G.arg.N);

% recalculate subindices based on Bloch sim results, in
% anticipation of the next update
%G.arg = Gocrf_subindcalc(G.arg,real(asin(abs(m))),G.arg.indmask);
jf_tmp = Gocrf_subindcalc(G.arg,real(asin(abs(m))),G.arg.indmask);
G = subsasgn_trick(G, 'arg', jf_tmp); % jf kludge for now

% then stick them into the Gmri objects
G = Gocrf_new_B_Ct(G,G.arg.a0B,G.arg.a0Ct,G.arg.b0B,G.arg.b0Ct,G.arg.a0,G.arg.b0);

if ~isempty(G.arg.RFbasis)
  % evaluate perturbations for basis functions
  for ll = 1:G.arg.Nbasis
    G.arg.Gat(:,ll) = G.arg.b0G * conj(G.arg.RFbasis(:,ll));
    G.arg.Gbt(:,ll) = G.arg.a0G * conj(G.arg.RFbasis(:,ll));
  end  
end


% 
% Gocrf_blochsim(): Perform a Bloch simulation
%
function [m,mz,a,b,a0Ct,b0Ct] = Gocrf_blochsim(B1,arg,pos,sens,we,simmask)

gambar = 42570000;               % gamma/2pi in Hz/T
gam = gambar*2*pi/10000;         % gamma in radians/g

if nargin == 6
  if ~isempty(simmask)
    % filter out un-sim'd spatial locs
    pos = pos(logical(simmask(:)),:);
    sens = sens(logical(simmask(:)),:);
    if ~isempty(we)
      we = we(logical(simmask));
    end
  end
end

Ns = size(pos,1);
Nt = size(arg.g,1);

statea = ones(Ns,1);
stateb = zeros(Ns,1);

% sum up RF over coils
B1 = reshape(B1,[size(arg.g,1) size(sens,2)]).';
bxy = sens * B1;

% sum up gradient over channels
bz = pos * arg.g';

% add off-resonance
if ~isempty(we) & sum(abs(we(:))) ~= 0
  bz = bz + repmat(we(:)/gam*2*pi,1,Nt);
end

if ~isempty(arg.fmod)
  bz = bz + arg.fmod/gam*2*pi;
end

tmp = zeros(2*Ns,1);

if nargout == 4
  % if we are only calling it with four output args,
  % assume we are doing the subspace sim and return all
  % states. Otherwise, we only return the final state.
  returnallstate = 1;
  a = zeros(Ns,Nt);
  b = zeros(Ns,Nt);
else; returnallstate = 0; end

% check if we need to do the running inner product to get spatial coeffs
if nargout == 6; 
  inprod = 1; 
  a0Ct = 0;b0Ct = 0;
  tvec = 0:arg.dt:(Nt-1)*arg.dt;
else; inprod = 0; end

ktmp = -[arg.k - repmat(arg.k(1,:),[size(arg.k,1) 1])];

for tt = 1:Nt
  
  phi = arg.dt*gam*(abs(bxy(:,tt)).^2+bz(:,tt).^2).^0.5;
  normfact = arg.dt*gam*(phi.^(-1));normfact(~isfinite(normfact)) = 0;
  nxy = normfact.*bxy(:,tt);nxy(~isfinite(nxy)) = 0;
  nz = normfact.*bz(:,tt);nz(~isfinite(nz)) = 0;
  cp = cos(phi/2);
  sp = sin(phi/2);
  alpha = cp+1i*nz.*sp;
  beta = 1i*conj(nxy).*sp;
  
  tmpa = alpha.*statea + beta.*stateb;
  tmpb = -conj(beta).*statea + conj(alpha).*stateb;
  
  statea = tmpa;stateb = tmpb;
  if any(~isfinite(statea)); keyboard; end;
  if any(~isfinite(stateb)); keyboard; end;
  if returnallstate
    a(:,tt) = statea;
    b(:,tt) = -conj(stateb); 
  end

  if inprod
    % build drf transform vector
    Asub = exp(-1i*2*pi/2*pos*ktmp(tt,:).');
    if ~isempty(we) & sum(abs(we(:))) ~= 0
      Asub = Asub.*exp(1i*2*pi/2*we(:)*tvec(tt));
    end
    if ~isempty(arg.fmod)
      Asub = Asub.*exp(1i*2*pi/2*arg.fmod*ones(size(pos,1),1)*tvec(tt));
    end
    % get B matrix coeffs
    a0Ct = a0Ct + conj((conj(statea).*Asub)*arg.a0B(tt,:));
    b0Ct = b0Ct - conj((stateb.*Asub)*arg.b0B(tt,:));
  end
  
end

% return final alpha, beta if not returning the whole 
% state progression
if ~returnallstate
  a = statea;
  b = -conj(stateb);
end

% calculate final magnetization state
mxy0 = arg.m0(1)+1i*arg.m0(2);
mz0 = arg.m0(3);
m = mz0*2*conj(statea).*stateb;
m = m + mxy0*conj(statea).^2;
m = m - conj(mxy0)*stateb.^2;
mz = mz0*(statea.*conj(statea) - stateb.*conj(stateb));
mz = mz + 2*real(mxy0*conj(statea).*-conj(stateb));

% embed results, if we did a masked sim
if nargin == 6
  if ~isempty(simmask)
    tmp = zeros(size(simmask(:)));tmp(logical(simmask)) = a;a = tmp;
    tmp = zeros(size(simmask(:)));tmp(logical(simmask)) = b;b = tmp;
    tmp = zeros(size(simmask(:)));tmp(logical(simmask)) = m;m = tmp;
    tmp = zeros(size(simmask(:)));tmp(logical(simmask)) = mz;mz = tmp;
    tmpa0Ct = zeros(length(simmask(:)),size(arg.a0B,2));
    tmpb0Ct = zeros(length(simmask(:)),size(arg.b0B,2));
    for ii = 1:size(arg.a0B,2)
      tmpa0Ct(logical(simmask),ii) = a0Ct(:,ii);
      tmpb0Ct(logical(simmask),ii) = b0Ct(:,ii);
    end
    a0Ct = tmpa0Ct;b0Ct = tmpb0Ct;
  end
end
