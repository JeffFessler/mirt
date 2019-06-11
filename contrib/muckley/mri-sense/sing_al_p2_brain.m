% set values, load data
nx = 144;
ny = 256;
nc = 8;
niter = 1500;
beta = 2^8.5;
nlevels = 3;

load braindat;
load kmask;
load brainxinf;

% generate and sample data
A = (1/sqrt(nx*ny))*Gdft('mask', true(nx,ny), 'fftshift', 1, ...
    'ifftshift', 1);

dat = A*reshape(coilim, [nx*ny nc]); clear coilim;

kmask = kmask((256-nx)/2+1:256-(256-nx)/2,:);
dat = col(dat(col(kmask),:));

n = kmask; b = 1:nx*ny;
n = b(n); clear b;

% don't shrink approximation coeffs
beta = beta*ones(nx,ny);
beta(1:nx/2^nlevels,1:ny/2^nlevels) = 0;
beta = col(beta);

% initial estimate
A = Apsm('knownfn', 'time', 'v', 1, 'n', n(:), 'smap', smap, 'immask', ...
    true(nx,ny), 'nk', nx*ny);
x = A'*dat; x = x./ col(sum(abs(smap).^2,3));

% build the system matrices
Q = cell(nc,1);
for i=1:nc
    Q{i} = (1/sqrt(nx*ny)) * ...
        Gdft('mask', true(nx,ny), 'fftshift', 1, 'ifftshift', 1);
end
Q = block_fatrix(Q, 'type', 'diag');
F = cell(nc,1);
for i=1:nc
	F{i} = (1/sqrt(nx*ny)) * ...
        Gdft('samp', kmask, 'fftshift', 1, 'ifftshift', 1);
end
F = block_fatrix(F, 'type', 'diag');
S = cell(nc,1);
for i=1:nc
    S{i} = Gdiag(col(smap(:,:,i)));
end
S = block_fatrix(S, 'type', 'col');
W = Godwt1(true(nx,ny), 'level', nlevels, 'wname', 'haar');

% set AL parameters
sos = sum(abs(smap).^2,3);
mu = 1/23;
if 0.9*max(col(sos))/min(col(sos)) < 12
    nu2 = (0.1*max(col(sos))*min(col(sos))) / ...
        (0.9*max(col(sos)) - min(col(sos)));
else
    nu2 = (max(col(sos)) - 12*min(col(sos)))/11;
end
nu1 = nu2;

% initializations
d = F'*dat(:);
q = S*x;
z = W*x;
w = x;
eta_q = zeros(size(q));
eta_z = zeros(size(z));
eta_w = zeros(size(w));

% build the Hessians
Hnu2 = Gdiag(1./(col(sos) + nu2));
kmaskrep = repmat(col(kmask), [nc 1]);
Hu = Gdiag(1./(kmaskrep + mu));
Hnu1nu2 = Gdiag(1./(ones(size(x)) + nu2/nu1));

% algorithm book-keeping
shrink = @(t, a) (t-a .* sign(t)) .* (abs(t) > a);
i = 0;
thetime(1) = 0;
xdist(1) = norm(col(x) - col(xinf));

% go
tic;
while i < niter
    printm('iteration %d of %d', i+1, niter);
    
    q = Q'*(Hu*(Q*(d + mu*(S*x + eta_q))));
    
    z = shrink(W*w + eta_z, beta./(mu*nu1));
    
    w = Hnu1nu2*(W'*(z - eta_z) + (nu2/nu1)*(x + eta_w));

    x = Hnu2*(S'*(q - eta_q) + nu2*(w - eta_w));
    
    eta_q = eta_q - (q - S*x);
    eta_z = eta_z - (z - W*w);
    eta_w = eta_w - (w - x);
        
    thetime(i+2) = toc;
    xdist(i+2) = norm(col(x) - col(xinf));
    i = i+1;
end
toc;
x = reshape(x, [nx ny]);
im(x);
