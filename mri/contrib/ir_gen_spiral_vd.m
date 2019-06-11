 function [Gx, Gy, kxvec, kyvec, sx, sy] = ir_gen_spiral_vd(D, N, nl, gamp, gslew, dsamp)
%function [Gx, Gy, kxvec, kyvec, sx, sy] = ir_gens_piral_vd(D, N, nl, gamp, gslew, dsamp)
%| Subroutine for generating slew-rate-limited spirals.
%| Returns length of resulting gradients, after adding downramp.
%| Craig Meyer, Copyright Leland Stanford Junior University, 1996.
%|
%| If you put 1 in for nl, you just get a standard full density spiral.
%|
%| variable density added by D. Noll, 10/25/01
%| 2015-06-09 cosmetic changes by JF, including built-in test, renaming

if nargin < 1, help(mfilename), error(mfilename), end
if streq(D, 'test'), ir_gen_spiral_vd_test, return, end

MAXDECRATIO = 32; % maximum allowed decimation of input ts

% double A, risetime,GAM,OM,S,absk,targetk;
% short delx, dely;
% int m, n, npts, dnpts, loop, res;
% int dentrans, den1, kdenrad;
% double decratio, om, s;
% double g0, thetan_1, theta, deltheta, taun_1, taun, tauhat;
% double absg,gtilde,B,t1,t2,t3,ac,tgx,tgy;
% double kx, ky, oldkx, oldky;
% double OMF, omf, denrad, scoffset, denoffset, scthat,fractrans, realn, ksv;

%%%%%%%%%%%%%% Initialize scanner variables %%%%%%%%%%%%%%%%%%%%%
GRESMAX = 21000;
if ~exist('nl','var')
	nl=1; % factor of undersampling in the outer regions of kspace
end
if ~exist('gamp','var')
	gamp=2.2; %3.50; % 2.2 for both 1.5 T and 3 T data
end
if ~exist('gslew','var')
	gslew=180; % 200 % 180 for 3T data and 120 (150) for 1.5 T data
end
nramp=0;
% nramp=100;
MAX_PG_WAMP=32766;

gts = 4e-06;

Tmax = GRESMAX*gts;

dts = gts;
opfov = D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GAM = 4257.0;
A = MAX_PG_WAMP;
risetime = gamp/gslew*10000;
OM = 2.0*pi/nl*D/(1/(GAM*gamp*dts));
S = (dts/1e-6)*A/risetime;
absk = 0.;
targetk = N/2;
OMF = OM*nl;
dentrans = dsamp/2;
S0=gslew*100;


ac = A;

loop = 1;
decratio = 1;

while (loop)
	loop = 0;
	om = OM/decratio;
	s = S/decratio;
	omf = OMF/decratio;
	den1 = 0;
	g0 = 0;
	Gx(1) = g0;
	Gy(1) = 0;
	absg = abs(g0); % absg = hypot(g0,0);
	oldkx = 0;
	oldky = 0;
	kx = Gx(1);
	ky = Gy(1);
	thetan_1 = 0;
	taun = 0;
	n = 1;
	while (absk < targetk) %keep generating waveform until reach desired k-space distance
		taun_1 = taun;
		taun = sqrt(kx.^2+ky.^2)/A; %taun = hypot(kx,ky)/A;
		tauhat = taun;
		realn = (n-1)/decratio; %realn = n/decratio
		if (realn > dsamp) %undersampled region
			if (den1 == 0)
				ksv = taun;
				den1 = 1;
			end
			if (realn > (dsamp+dentrans))
				scoffset = scthat;
				denoffset = taun_1;
				scthat = scoffset + om*(tauhat - denoffset);
				fractrans = 1;
			else
				scoffset = scthat;
				denoffset = taun_1;
				fractrans = (realn - dsamp)/( dentrans);
				fractrans = 1 - ( (fractrans-1)*(fractrans-1));
				scthat = scoffset + (omf + (om-omf)*fractrans)*(tauhat - denoffset);
			end
		else %fully sampled region
			scthat = omf*tauhat;
			fractrans = 0;
		end

		theta = atan2(scthat,1.0)+scthat;
		if (absg < ac)
			deltheta = theta-thetan_1;
			B = 1.0/(1.0+tan(deltheta)*tan(deltheta));
			gtilde = absg;
			t1 = s*s;
			t2 = gtilde*gtilde*(1-B);
			if (t2 > t1)
				decratio = decratio * 2.0;
				if (decratio > MAXDECRATIO)
					error('Genspiral failed');
					return;
				end
				loop = 1;
				break;
			end
			t3 = sqrt(t1-t2);
			absg = sqrt(B)*gtilde+t3;
			if (absg > ac)
				absg = ac;
			end
		end

		tgx = absg*cos(theta);
		tgy = absg*sin(theta);
		kx = kx + tgx;
		ky = ky + tgy;
		thetan_1=theta;
		if ~mod(n,round(decratio)) %IRINT(decratio)
			m = n/round(decratio);
			absk = nl*om*taun/(2*pi);
			if (absk > targetk)
				break;
			end
			Gx(m+1) = round((kx-oldkx)/decratio); %&0xfffe;
			Gy(m+1) = round((ky-oldky)/decratio); %&0xfffe;
			oldkx = kx;
			oldky = ky;
		end
		n=n+1;
	end
end

gres1 = m;
npts = m+1;

g = (Gx + i.*Gy)./(MAX_PG_WAMP/gamp); %slew rate vector
s = diff(g)./(gts*1000); % grad vector
Kx = cumsum([0 Gx])*gts*opfov*GAM./(MAX_PG_WAMP/gamp);
Ky = cumsum([0 Gy])*gts*opfov*GAM./(MAX_PG_WAMP/gamp);
k = Kx + 1i.*Ky; % kspace vector
if im
	pr max(abs(k))*2 % matrix
	pr max(abs(g)) % max gradient
	pr max(abs(s)) % max slew rate
end

kxvec = real(k);
kyvec = imag(k);
sx = real(s);
sy = imag(s);

end % ir_gen_spiral_vd


function ir_gen_spiral_vd_test
im plc 1 2
ax = [-1 1 -1 1] * 33;

% generate a fully sampled acquisition for the first 900 samples
% and then transition over the next 450 samples
% to a 3-fold reduction in the number of samples.
[Gx, Gy, kxvec, kyvec, sx, sy] = ir_gen_spiral_vd(24, 64, 3, 2.2, 180, 900);
im subplot 1
plot(kxvec, kyvec, '.'), axis(ax), axis square

[Gx, Gy, kxvec, kyvec, sx, sy] = ir_gen_spiral_vd(24, 64, 1, 2.2, 180, 900);
im subplot 2
plot(kxvec, kyvec, '.'), axis(ax), axis square

end % ir_gen_spiral_vd_test

%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE!
%%%%%%%%%%%%%%%%%%%%%%

%{
res = ceil(2*nl*om*taun/(2*pi));
kdenrad = ceil(2*nl*om*ksv/(2*pi));
sprintf('resolution = %d, high dens region = %d\n', res,kdenrad);

%/* add downramp */
t1 = atan2(Gy(npts-1), Gx(npts-1))+pi;
delx = round(S*cos(t1));%rint(S*cos(t1));
dely = round(S*sin(t1));%rint(S*sin(t1));
if (delx < 0)
	delx = (delx+1); % & 0xfffe;
else
	delx = delx; % & 0xfffe;
end
if (dely < 0)
	dely = (dely+1);% & 0xfffe;
else
	dely = dely; % & 0xfffe;
end
m = abs(Gx(npts-1)/delx); %ramp length needed for x
n = abs(Gy(npts-1)/dely); %ramp length needed for y
if(m>n)
	dnpts = npts+n;
else
	dnpts = npts+m;
end
tgx = Gx(npts-1);
tgy = Gy(npts-1);
%add the ramp points to the end of the gradients
for n=npts:dnpts
	tgx = tgx + delx;
	tgy = tgy + dely;
	Gx(n) = tgx;
	Gy(n) = tgy;
end
Gx(n) = 0;
Gy(n) = 0;

%Check if last point is zero.  If not make it zero or ramp down to zero of
%not already very close.  (hmmm, seems redundant... you guys)
% while ((Gx(n) ~= 0)|(Gy(n) ~= 0))
%
%	if (abs(Gx(n)) < abs(delx))
%		Gx(n+1) = 0;
%	else
%		Gx(n+1) = Gx(n)+delx;
%	end
%	if (abs(Gy(n)) < abs(dely))
%		Gy(n+1) = 0;
%	else
%		Gy(n+1) = Gy(n)+dely;
%	end
%	n=n+1;
% end
gres = n+1;


totx=0;
toty=0;
for j=1:gres-1
	totx=totx+Gx(j)/(MAX_PG_WAMP/gamp);
	toty=toty+Gy(j)/(MAX_PG_WAMP/gamp);
	kxvec(j)=totx;
	kyvec(j)=toty;
end

%******************** now add rephaser ****************************/

tgx = 0;
tgy = 0;
t1 = atan2( kyvec(gres-1), kxvec(gres-1))+ pi;
tempdelx = (.7*S0*cos(t1)); %?? is .7 to limit eddy current?
tempdely = (.7*S0*sin(t1));

nr = ceil(0.7*gamp/(.7*S0*dts)); %number of points for ramp of rephaser

adkx = (abs(kxvec(gres-1)) - abs(tempdelx*dts*nr*nr));
adky = (abs(kyvec(gres-1)) - abs(tempdely*dts*nr*nr));


if ((adkx < 0.0) | (adky < 0.0))
	m = round(sqrt(abs(kxvec(gres-1)/(tempdelx*dts))));
	n = round(sqrt(abs(kyvec(gres-1)/(tempdely*dts))));

	if (m > n)
		nr=m;
	else
		nr = n;
	end

	apareax=abs(tempdelx*dts*nr*nr);
	ratiox=kxvec(gres-1)/apareax;
	apareay=abs(tempdely*dts*nr*nr);
	ratioy=kyvec(gres-1)/apareay;

	tempdelx=tempdelx*(abs(ratiox));
	tempdely=tempdely*(abs(ratioy));

	if (tempdelx < 0)
		delx = (( tempdelx)+1) ; % & 0xfffe;
	else
		delx = ( tempdelx) ; % & 0xfffe;
	end
	if (tempdely < 0)
		dely = (( tempdely)+1) ; % & 0xfffe;
	else
		dely = (( tempdely)) ; % & 0xfffe;
	end

	sprintf('Tri rephaser %d\n',2*nr);

	for n=gres:(gres+nr-1)
		tgx = tgx+delx*dts;
		tgy = tgy+dely*dts;
		Gx(n) = tgx*(MAX_PG_WAMP/gamp);
		Gy(n) = tgy*(MAX_PG_WAMP/gamp);
		totx = totx+Gx(n)/(MAX_PG_WAMP/gamp);
		kxvec(n)= totx;
		toty = toty+Gy(n)/(MAX_PG_WAMP/gamp);
		kyvec(n)= toty;

	end

	for n = gres+nr:(gres+2*nr-1)
		tgx = tgx-delx*dts;
		tgy = tgy-dely*dts;
		Gx(n) = tgx*(MAX_PG_WAMP/gamp);
		Gy(n) = tgy*(MAX_PG_WAMP/gamp);

		totx = totx+Gx(n)/(MAX_PG_WAMP/gamp);
		kxvec(n)= totx;
		toty = toty+Gy(n)/(MAX_PG_WAMP/gamp);
		kyvec(n)= toty;

	end

	Gx(gres+2*nr) = 0;
	Gy(gres+2*nr) = 0;

	kxvec(gres+2*nr)= kxvec(gres+2*nr-1);
	kyvec(gres+2*nr)= kyvec(gres+2*nr-1);

	% if */

else
	m = round(abs(adkx/(tempdelx*dts*nr)));
	n = round(abs(adky/(tempdely*dts*nr)));

	if (m > n)
		nc = m;
	else
		nc = n;
	end

	apareax=abs(tempdelx*dts*nr*(nr +nc));
	ratiox=kxvec(gres-1)/apareax;
	apareay=abs(tempdely*dts*nr*(nr +nc));
	ratioy=kyvec(gres-1)/apareay;

	tempdelx=tempdelx*(abs(ratiox));
	tempdely=tempdely*(abs(ratioy));

	if (tempdelx < 0)
		delx = (( (tempdelx))+1) ; % & 0xfffe;
	else
		delx = ( (tempdelx)) ; % & 0xfffe;
	end
	if (tempdely < 0)
		dely = (( (tempdely)) +1) ; % & 0xfffe;
	else
		dely = ( (tempdely)) ; % & 0xfffe;
	end


	for n = gres:(gres+nr-1)
		tgx = tgx+delx*dts;
		tgy = tgy+dely*dts;
		gxi=tgx*MAX_PG_WAMP/gamp;
		gyi=tgy*MAX_PG_WAMP/gamp;
		Gx(n) = 2*(gxi/2);
		Gy(n) = 2*(gyi/2);

		totx = totx+tgx;
		kxvec(n)= totx;
		toty = toty+tgy;
		kyvec(n)= toty;
	end

	for n = (gres+nr):(gres+nr+nc-1)
		gxi=tgx*MAX_PG_WAMP/gamp;
		gyi=tgy*MAX_PG_WAMP/gamp;
		Gx(n) = 2*(gxi/2);
		Gy(n) = 2*(gyi/2);


		totx = totx+tgx;
		kxvec(n) = totx;
		toty = toty+tgy;
		kyvec(n) = toty;
	end


	ctx=tgx;
	cty=tgy;


	for n = (gres+nc+nr):(gres+nc+2*nr-1)
		tgx = tgx-delx*dts;
		tgy = tgy-dely*dts;
		gxi=tgx*MAX_PG_WAMP/gamp;
		gyi=tgy*MAX_PG_WAMP/gamp;
		Gx(n) = 2*(gxi/2);
		Gy(n) = 2*(gyi/2);

		totx = totx+tgx;
		kxvec(n)= totx;
		toty = toty+tgy;
		kyvec(n)= toty;
	end

	Gx(gres-1+nc+2*nr) = 0;
	Gy(gres-1+nc+2*nr) = 0;

	kxvec(gres-1+nc+2*nr)=kxvec(gres-1+nc+2*nr-1);
	kyvec(gres-1+nc+2*nr)=kyvec(gres-1+nc+2*nr-1);

end

	gres=n;
	Gxr=Gx(end:-1:1);
	Gyr=Gy(end:-1:1);

%[Gx,Gy,Gxr,Gyr,kx,ky]=genspiralvd_rev(22,64,2,2.2,180,300); figure, subplot(121),plot(Gx/MAX_PG_WAMP*gamp),subplot(122),plot(diff(Gx/MAX_PG_WAMP*gamp/10*1000*25));

%}
