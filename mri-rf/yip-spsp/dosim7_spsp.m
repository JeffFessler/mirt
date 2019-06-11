 function [mallsum] = dosim7_spsp(kp,rfp,iop)
%function [mallsum] = dosim7_spsp(kp,rfp,iop)
%|
%|This script put the B1 and gradient waveforms from files
%|(.mag, .ph, .gx, .gy, .gz) into the bloch simulator, and then
%|reads the output of the bloch simulator.
%|
%| in
%| kp	k-space trajectory parameter structure
%| rfp	RF pulse computation parameter structure
%| iop	input-output parameter structure containing Bloch simulation parameters
%|
%| out
%| mallsum	complex-valued Bloch simulation result;
%|		contains only the transverse plane component of magnetization
%|
%| Chun-yu Yip, 4/2/2009

nshots = 1; % variable not used any more (but needed)
pulsen = 1; % variable not used any more (but needed)
signalscale = 2^14-1;
filescale = 32767;

%com = 'blochmulti_spsp_03Dec08';
% jf version to support multiple OS types
arch = mexext; arch = arch(4:end);
dir = path_find_dir('mri-rf/yip-spsp');
com = sprintf('%s/bloch/blochmulti_spsp-%s', dir, arch);

printm('Doing Bloch simulation with %s', com)

exestr = sprintf('%s %s %s %d %d %d %f %d %d %d %d %d %d %f %d %f %1.8f', ...
	com, ...
	[iop.waveformfilespath iop.infname], ...
	[iop.waveformfilespath iop.outfname_from_sim], ...
	rfp.flipangle, iop.bdimf, iop.bdimz, iop.bfovz, kp.npnts, ...
	iop.bfovf, rfp.dfovf, rfp.dfovz, nshots, pulsen, ...
	(kp.npnts-1)*kp.pointtime, kp.npnts, kp.gmax, kp.pointtime)

%[status,originalpath] = system('pwd');
%cd(iop.waveformfilespath);
unix(exestr); % Run the Bloch simulation code

m=zeros([iop.bdimf iop.bdimz nshots]);
ph=zeros([iop.bdimf iop.bdimz nshots]);
mag=zeros([iop.bdimf iop.bdimz nshots]);
mall=zeros([iop.bdimf iop.bdimz nshots]);
magall=zeros([iop.bdimf iop.bdimz nshots]);
phall=zeros([iop.bdimf iop.bdimz nshots]);

% Read in images

outfname_full = [iop.waveformfilespath iop.outfname_from_sim '.' num2str(iop.bdimf) '.' num2str(iop.bdimz) '.' num2str(nshots) '.dat'];
fid = fopen(outfname_full,'r');
temp = fread(fid,2*3*iop.bdimf*iop.bdimz*nshots,'short') / filescale;
fclose(fid);
mtemp = reshape(temp,[3 iop.bdimf iop.bdimz nshots]);

%cd(strtrim(originalpath));

% Build output files
for sindex=1:nshots
	temp=squeeze(mtemp(1,:,:,sindex)+1i*mtemp(2,:,:,sindex));
	m(:,:,sindex)=temp;
	mag(:,:,sindex)=abs(temp);
	ph(:,:,sindex)=atan2(imag(temp),real(temp));
	mall(:,:,sindex)=reshape(squeeze(m(:,:,sindex)),[iop.bdimf iop.bdimz]);
	magall(:,:,sindex)=reshape(squeeze(mag(:,:,sindex)),[iop.bdimf iop.bdimz]);
	phall(:,:,sindex)=reshape(squeeze(ph(:,:,sindex)),[iop.bdimf iop.bdimz]);

end


% Add the shots together
if (nshots>1)
	msum=squeeze(sum(m,3));
	mallsum=squeeze(sum(mall,3));
else
	msum=squeeze(m);
	mallsum=squeeze(mall); % matrix holds the mag results
end
m=squeeze(m);

phsum=atan2(imag(msum),real(msum));
phallsum=atan2(imag(mallsum),real(mallsum)); % matrix holds the phase results

if iop.show_blochsimresults
	bresz = iop.bfovz/iop.bdimz; %cm
	bresf = iop.bfovf/iop.bdimf; %Hz
	z = [-iop.bfovz/2:bresz:iop.bfovz/2-bresz];
	f = [-iop.bfovf/2:bresf:iop.bfovf/2-bresf];
end

if 1
	pl = @(i) subplot(3,4,i);
	pl(9)
	plot(z,abs(mallsum(iop.bdimf/2,:)))
	title('Water slice profile (magnitude)')
	xlabel('z position (cm)')
	ylabel('Normalized excitation magnitude')
	axis([z(1) z(end) 0 1])
	grid

	pl(10)
	plot(z,angle(mallsum(iop.bdimf/2,:)),'r')
	xlabel('z position (cm)')
	ylabel('radians')
	title('Water slice profile (phase)')
	ytick([0 0.5 1])
	grid

	pl(4)
	imagesc(f,z,abs(mallsum).')
	ylabel('z (cm)')
	xlabel('Frequency offset (Hz)')
	colormap default
	title('Bloch simulated SPSP pattern (magnitude)')
	colorbar

	pl(8)
	imagesc(f,z,angle(mallsum).')
	ylabel('z (cm)')
	xlabel('Frequency offset (Hz)')
	title('Bloch simulated SPSP pattern (phase)')
	colormap(gca, hsv)
	colorbar
end
