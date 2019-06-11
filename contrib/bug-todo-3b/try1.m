function TestPointBaseProjector(ProjectorType,bBlockMatrix)
%|function TestPointBaseProjector(ProjectorType,bBlockMatrix)
%|
%| Test Gtomo3 object, for 3D ESRI image reconstruction.
%| Test 3d forward and backprojection using mex files.
%|
%| in
%|	ProjectorType	'b'|'a'      default: 'b'
%|      'b' is 3b improved projector
%|      'a' is 3a the first version projector
%|	bBlockMatrix	false|true   default: false
%|      true - designed to work with ordered-subsets algorithms

if (nargin==0)
	ProjectorType='b';
	bBlockMatrix=false;
end

%Xtrue init
ig3 = image_geom('nx', 32, 'ny', 34, 'nz', 6, 'dx', 1, 'dy', 1);
xtrue = [
	0 0 0	ig3.nx/2-4 ig3.ny/2-4 inf 0 0 1;
	5 0 0	ig3.nx/15*[1 1 1] 0 0 1;
	4 6 2	ig3.nx/13*[1 1 1] 0 0 1;
	-3 0 0.5 ig3.nx/17*[1 1 1] 0 0 1];
xtrue = ellipsoid_im(ig3, xtrue, 'oversample', 2);
im plc 3 3
im(1, xtrue, 'true images'), cbar
ig3.mask = ig3.circ(ig3.nx/2-2, ig3.ny/2-2, 0, 0);
im(2, ig3.mask, 'support mask'), cbar

%Two type of projectors
if (ProjectorType=='b')
	f.sys_type=FldWrite3b();
else
	f.sys_type=FldWrite3a();
end

%Build the system matrix
f.chat = 1; % Verbosity on
f.nthread = 2;
G = Gtomo3(f.sys_type, ig3.mask, ig3.nx, ig3.ny, ig3.nz, ...
        'nthread', f.nthread, 'chat', f.chat);

% fld_write('x.fld', xtrue)
% G

%Test the forward and backward projector
proj = G*xtrue;
im(3, proj)

XEstim = G'*proj;
im(4, XEstim)
return

f.nblock = 20;
if (bBlockMatrix==true)%
	G = Gblock(G, f.nblock, 0);
end

return

%Reconstructed using pcg
wi = ones(size(proj));
W = diag_sp(wi(:));
R = Reg1(ig3.mask, 'beta', 2^7);
xinit=ones(size(xtrue));
BinaryMask3d=logical(ig3.mask);
[xs, info] = pwls_pcg1(xinit(BinaryMask3d), Gb, W, proj(:),R, ...
    'niter',  3, 'isave',  [1:3],...
    'precon', 1);

xs = ig3.embed(xs);

im(3, xs, 'Reconstructed'), cbar
end


function CmdStr=FldWrite3b(FileName,NView)
Dir=test_dir;
if(nargin==0)
    FileName='viewgeom.fld';
    NView=10;
end
FileName=[Dir FileName];
nu=140;
%nv=1;
nv=8; % jf trying this
nview=NView;
StrEnd=['@-' num2str(nu) ',' num2str(nv) ',' num2str(nview)];
StrStart=['@' FileName];
Strb='3b@-';
CmdStr=[Strb StrStart StrEnd];
%Row vectors:
%%Theta=linspace(0,2*pi,NView); % jf says: we should have -pi/2 < theta < pi/2
Theta = zeros(1,NView);
Phi=linspace(0,2*pi,NView);
su=ones(1,NView);
sv=ones(1,NView);
cu=zeros(1,NView);
cv=zeros(1,NView);
iv_min=zeros(1,NView);
iv_max=nv*ones(1,NView);
%Put all the vectors in 8xNView data
Data=[Theta; Phi; su; sv; cu; cv; iv_min; iv_max];
fld_write(FileName, Data);
end


function CmdStr=FldWrite3a(FileName,NView)
Dir=test_dir;
if(nargin==0)
    FileName='angles.fld';
    NView=10;
end
FileName=[Dir FileName];
nu=140;
nv=1;
su=1;sv=1;sz=1;

StrEnd=['@' FileName];
StrStart=['@' '-'];
Strb=['3a@' num2str(nu) ',' num2str(nv) ',' num2str(su) ',' num2str(sv) ',' num2str(sz) ];
CmdStr=[Strb StrStart StrEnd];

%%Theta=linspace(0,2*pi,NView);
Theta = zeros(1,Nview);
Phi=linspace(0,2*pi,NView);
Data=[Theta;Phi];
fld_write(FileName, Data)
end
