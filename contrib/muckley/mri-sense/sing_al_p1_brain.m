% set values, load data
nx = 144;
ny = 256;
nc = 8;
niter = 250;
mu = 2^-3;
nlevels = 3;
beta = 2^8.5;

load braindat;
load kmask;
load brainxinf;
% load braindb4xinf;
mask = true(nx,ny);

% mask = abs(xinf) > 0.1*max(col(abs(xinf)));

% generate and sample data
A = (1/sqrt(nx*ny))*Gdft('mask', true(nx,ny), 'fftshift', 1, ...
    'ifftshift', 1);

dat = A*reshape(coilim, [nx*ny nc]); clear coilim;

kmask = kmask((256-nx)/2+1:256-(256-nx)/2,:);
dat = col(dat(col(kmask),:));

n = kmask; b = 1:nx*ny;
n = b(n); clear b;

% build system matrices
A = Apsm('knownfn', 'time', 'v', 1, 'n', n(:), 'smap', smap, 'immask', ...
    true(nx,ny), 'nk', nx*ny);
W = Godwt1(mask, 'level', nlevels, 'wname', 'haar');
P = Gdiag(1./(col(sum(abs(smap).^2,3)) + mu));

% don't regularize approx coeffs
beta = beta*ones(nx,ny);
beta(1:nx/2^nlevels,1:ny/2^nlevels) = 0;
beta = col(beta);

% initializations
x = A'*dat; x = x./ col(sum(abs(smap).^2,3));
x = x(mask);
z = W*x;
eta = zeros(size(z));

% system matrix
A = Apsm('knownfn', 'time', 'v', 1, 'n', n(:), 'smap', smap, 'immask', ...
    mask, 'nk', nx*ny);

% algorithm book-keeping
shrink = @(t, a) (t-a .* sign(t)) .* (abs(t) > a);
i = 0;
thetime(1) = 0;
xdist(1) = norm(col(x) - col(xinf(mask)));

% go
tic;
while i < niter
    printm('iteration %d of %d', i+1, niter);
    
    z = shrink(W*x + eta, beta./mu);
    x = qpwls_pcg1(x, [A; sqrt(mu)*W], 1, [dat(:); sqrt(mu)*(z - eta)], ...
        0, 'niter', 5, 'precon', P);
    eta = eta - (z - W*x);
        
    thetime(i+2) = toc;
    xdist(i+2) = norm(col(x) - col(xinf(mask)));
    i = i+1;
end
toc;
x = embed(x, mask);
im(x)
