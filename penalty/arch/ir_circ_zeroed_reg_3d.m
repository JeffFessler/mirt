 function [C,B,phi,varargout] = ir_circ_zeroed_reg_3d(order, maskR, sparseR)
%function [C,B,phi,varargout] = ir_circ_zeroed_reg_3d(order, maskR, sparseR)
% Input: order - order of regularizer
%		maskR - [nx,ny] - mask of regularizer (changes C,B,and phiR)
%		sparseR - binary flag for a sparse output (varargout) instead of fatrix2.
%
% Output: C - [M,nx ny (or) K] - circulant matrix (either fatrix2 of sparse depending on flag)
%			B - [M,1] - binary mask to remove edge wrapping
%			phi - [nx ny,1] - convolution kernal of C'C (e.g., C'C = IFFT2 phi FFT2)
%			vargout - [M,nx ny] - R = BC in fatrix2 form.
%
% I used ifftshift in computing phi/phiR since you can have odd images.

nx = size(maskR,1);
ny = size(maskR,2);
nz = size(maskR,3);

if sparseR == 1
	printm('Using sparse matrix to create R (no C,B,phi returned)')

	argsR = {maskR, 'beta', 1, 'order', order, 'distance_power', 2, ...
		'type_diff', 'spmat', 'type_penal', 'mat'};
	R = Reg1(argsR{:});

	% set B, C, phi to null so that people are not confused.
	C = {};
	B = {};
	phi = {};

else
	if order == 2 %likely make this auto using convolution of basic parts

        psfH = zeros(3,3,3);
        psfH(:,:,2) = [0 0 0;-1 2 -1;0 0 0];
        psfV = zeros(3,3,3);
        psfV(:,:,2) =  [0 -1 0;0  2 0; 0 -1 0];
        psfA = zeros(3,3,3);
        psfA(2,2,:) = [-1 2 -1];
        
	else
		error('only coded arg.order == 2 -- see comment');
	end

	CH = Gblur(maskR, 'psf', psfH, 'type', 'imfilter,circ', 'class', 'fatrix2');
	CV = Gblur(maskR, 'psf', psfV, 'type', 'imfilter,circ', 'class', 'fatrix2');
	CA = Gblur(maskR, 'psf', psfA, 'type', 'imfilter,circ', 'class', 'fatrix2');
	C = [CV;CH;CA]; %clear CH; clear CV;

	% design B for Rs (can be done with only Rweights)
	B = ones(nx,ny,nz);
	BH = B; BH(:,[1,end],:) = 0;
	BV = B; BV([1,end],:,:) = 0;
	BA = B; BA(:,:,[1,end]) = 0;
	BH = BH(:); BV = BV(:); BA = BA(:);
	B = [BV;BH;BA];

	Rw = Rweights(maskR,[1 nx nx*ny],'type_wt','array','edge_type', ...
				'tight','beta',1,'order',2,'distance_power',0);

	B = B .* Rw(:); % if full maskR then B = Rw(:).

	%% compute phi when maskR is all ones.
	%it is faster to just the first column C'C and not use ifftshift
	% but use Fessler's approach to be consistent with toolkit.
	if isempty(find(maskR == 0))
		ig = image_geom_mri('nx',nx,'ny',ny,'nz',nz,'dx',1);
		phi = fftn(ifftshift(reshape(C' * (C * double(ig.unitv(:))), size(maskR))));
		phi = reale(phi);
		phi = max(phi,0);% not really needed for order 2
		phi = phi(:);
	else
		warning('No phi is given when maskR is not all ones');
		phi = [];
	end

end

if nargout == 4
	if ~sparseR
		B = reshape(B, [size(maskR) 3]);
		Bd = Gdiag(B, 'class', 'fatrix2');
		varargout{1} = Bd * C; % implementation is slower than B * where B is a vector and not a fatrix2.
	else
		varargout{1} = R.C;
	end
end


%% compare to old method -- note previous output
if 0 %not perfect, error on order of e-7.

Bm = [BV;BH;BD;BD];
if nargout ~= 4
	B = Gdiag(B, 'class', 'fatrix2');
	varargout{1} = B * C;
end

%old regularizer
mask = true(nx,ny);
args = {mask, 'beta', 1, 'order', order, 'distance_power', 2};
R = Reg1(args{:});
C1 = R.C;
omega = fft2(fftshift(reshape(R.cgrad(R,ig.unitv(:)), [nx ny]))); % R.cgrad = C'C with no added weights
omega = reale(omega);

% check reg ability
tmp = zeros(nx,ny); %rand(nx,ny); %zeros(nx,ny);
tmp(5,5) = 1;
tmp = tmp(:);
CC = varargout{1};
fd = varargout{1} * tmp;
fi = CC' * fd;
fd = reshape(fd, [nx ny 4]); %B .* (C *tmp); fd = reshape(fd, [nx
									%ny 4]);
fi = reshape(fi, [nx ny]);
figure; im(fd(1:10,1:10,:)); %im(fd);
fd1 = C1 * tmp;
fi1 = C1' * fd1;
fd1 = reshape(fd1, [nx ny 4]);
fi1 = reshape(fi1, [nx ny]);
figure; im(fd1(1:10,1:10,:)); %im(fd1);
figure; im(fd(1:10,1:10,:) - fd1(1:10,1:10,:)); %im(fd - fd1);

figure; im(fi(1:10,1:10)); %im(fd);
figure; im(fi1(1:10,1:10)); %im(fd1);
figure; im(fi(1:10,1:10) - fi1(1:10,1:10)); %im(fd - fd1);

%timing
varargout{1} * tmp;
tic;
for i = 1:100
varargout{1} * tmp;
end
disp(['Double fatrix2 time for 100 -- ' num2str(toc)]);
Bm .* (C * tmp);
tic;
for i = 1:100
Bm .* (C * tmp);
end
disp(['.* and fatrix2 time for 100 -- ' num2str(toc)]);
C1 * tmp;
tic;
for i = 1:100
C1 * tmp;
end
disp(['Original single fatrix2 time for 100 -- ' num2str(toc)]);


%test to see if C is masked or not
if 0 % it does not appear to be
	notS = [];
	tmps = find(maskR == 1);
	for i = 1:length(tmps)
		if (norm(C(:,i) - Cf(:,tmps(i))) > eps)
			notS = [notS,i];
		end
	end
end

end
