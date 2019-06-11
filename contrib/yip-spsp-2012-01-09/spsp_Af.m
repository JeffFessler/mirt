function ob = spsp_Af(z,f,kz,kf,tsp,fovf,fovz)
%function ob = spsps_Af(z,f,kz,kf,tsp,fovf,fovz)
%
% input 
% z : Nz by 1 vector  (cm)
% f : Nf by 1 vector  (Hz)
% kz, kf : M by 1 vector
% tsp : sampling time of RF pulse in sec
% Construct spsp_Af object
%
%
%This one calls the fast projections!
%
% Copyright, Sangwoo Lee, University of Michigan, 2005

%  keyboard;

arg.z=z(:);
arg.f=f(:);
arg.kz=kz(:);
arg.kf=kf(:);
arg.tsp=tsp;
arg.Nz=length(arg.z);
arg.Nf=length(arg.f);
arg.M=length(arg.kz);
arg.gamma=26751; %rad/sec/g
arg.fovf=fovf; % Hz
arg.fovz=fovz; % cm



nufft_args= {[arg.Nz arg.Nf], [8 8], 2*[arg.Nz arg.Nf], [arg.Nz arg.Nf]/2, 'table', 2^10, 'minmax:kb'};

mask=logical(ones(arg.Nz,arg.Nf));
arg.A=Gmri([kz(:)*fovz/arg.Nz,kf(:)*fovf/arg.Nf],mask,'fov',[arg.Nz arg.Nf], 'nufft',nufft_args);;

if arg.M ~= length(arg.kf)
    error('the length of kz and kf doesn''t match');
end;
arg.dim=[arg.Nz*arg.Nf arg.M];
ob = Fatrix(arg.dim, arg, 'caller', mfilename, ...
        'forw', @spsp_Af_forw, 'back', @spsp_Af_back);

 
function y = spsp_Af_forw(arg, x)

x=complexify(x);
y=i*arg.gamma*arg.tsp*(arg.A'*x);


function x = spsp_Af_back(arg, y)
y=complexify(y);
x = -i*arg.gamma*arg.tsp*(arg.A*y);   
