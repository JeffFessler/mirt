function G = Grid2GlinearROI(Grid,Grid_dec,npxdtc,ctrIdx,block,interp,nblk,vmsk,nnzmsk)
% Grid2Glinear
% input Grid/s generated in fbkp_geometry
% -Grid is a matrix of integers that are coordinate of the detector pixel
% -Grid_dec is [] is the 'nearest neighbour' interp option is used
% -block (integer & scalar) sets the projection undergoing backprojection
% -nblk number of blocks (i.e. projections) used to make up the (rows) sparse matrix
% i.e. the the large block goes from block (i.e. projection) to block + nblk
% -vmsk is the mask as a vector (logical)
% -npxdtc number of detector pixels ( nb in Glinear)
% output one sparse matrix G as by Glinear (see Fessler code)
%
% Gianni Schena March  2003

G=[];
block
nargin

Nx=size(Grid(:,:,1),1); Ny=size(Grid(:,:,1),2); % get the size of the problem
nm=Nx*Ny
%

if nargin < 8
	vmsk=logical(ones(1,nm)); nnzmsk=nm;
	% define default mask - if the mask is not an input parameter
end

if nargin < 7, nblk=1, end

%SPALLOC(M,N,NZMAX) creates an M-by-N all zero sparse matrix with room to eventually hold NZMAX nonzeros.

if strcmp(interp, 'nearest neighbor')
	G = spalloc(nblk*npxdtc,nnzmsk, nblk*nnzmsk);
else
	G = spalloc(nblk*npxdtc,nnzmsk, 2*nblk*nnzmsk); % pre-allocation
end

jc=[1:1:nnzmsk]; % index of columns in the matrix G

for ib=1:nblk % for nblk projections and starting from block

	Gv=Grid(:,:,block+ib-1) ; Gv=Gv(:) ; % starts from block
	ir = double(Gv) + ctrIdx  ; % the coordinate (saved as 'round') + centre
	ir=[ir+npxdtc*(ib-1)]; ir=ir(vmsk)';

	% S = SPARSE(i,j,s,m,n,nzmax) uses the rows of [i,j,s] to generate an
	% m-by-n sparse matrix with space allocated for nzmax nonzeros.

	if strcmp(interp, 'nearest neighbor')
		v=1.;
		T=sparse(ir,jc,double(v),npxdtc*nblk,nnzmsk);

	else % i.e. interp == linear
		dec=double(Grid_dec(:,:,block+ib-1))/100; % use also the decimal saved as integer
		vc = dec(:); % weight for (floor +1) == weight of 'ceil' term
		vc=vc(vmsk);
		T= [sparse(ir,jc,(1-vc),npxdtc*(nblk),nnzmsk) + sparse(ir+1,jc,vc,npxdtc*(nblk),nnzmsk)];
		%T= [sparse(ir,jc',(1-vc)) + sparse([ir+1],jc',vc)];
	end

	G=[G+T];
end % end ib
