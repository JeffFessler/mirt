 function write2files(kp,rfp,iop,gz,b)
%function write2files(kp,rfp,iop,gz,b)
%|
%| Function to generate waveform files based on waveforms gz and b, as
%| specified by parameter structure iop.
%|
%| in
%| kp	k-space trajectory parameter structure
%| rfp	RF waveform design parameter structure
%| iop	input-output parameter structure
%| gz	z gradient waveform vector
%| b	complex-valued SPSP ERF waveform vector
%|
%| Chun-yu Yip, 4/1/2009

npnts = length(gz);	%Number of pulse samples
nshots = 1;		%variable not used any more (existence needed)
pulsen = 1;		%variable not used any more (existence needed)
signalscale = 2^14-1;	%constant for scaling waveforms properly
eosloc = npnts*nshots*pulsen; %End Of Sequence LOCation

gx = zeros(size(gz));
gy = zeros(size(gz));
g = [gx gy gz];

% Waveforms on GE Signa (3T) must have even-length
if mod(npnts,2)==1
	b = [b;0];
	g = [g; [0 0 0]];
end

magb1 = abs(b)/max(abs(b));     %Normalize of RF pulse; in scanner
                                %experiment, proper scaling of
                                %RF pulse for the desired flip angle should
                                %be done on the pulse sequence level.
magb1_4sim = abs(b);            %No normalization for simulation
phb1 = atan2(imag((b)),real((b))); %phase in radian [-pi, pi]

% Scale for signa, and put eos bit:
%All non-EOS samples has last bit as 0
%EOC sample has last bit as 1
magb1_final=round(magb1*signalscale)*2;
magb1_final(eosloc)=magb1_final(eosloc)-1;
magb1_4sim_final=round(magb1_4sim*signalscale)*2;
magb1_4sim_final(eosloc)=magb1_4sim_final(eosloc)-1;
phb1_final=round(phb1/pi*signalscale)*2;
phb1_final(eosloc)=phb1_final(eosloc)-1;
gx_final=round(gx*signalscale/kp.gmax)*2;
gx_final(eosloc)=gx_final(eosloc)-1;
gy_final=round(gy*signalscale/kp.gmax)*2;
gy_final(eosloc)=gy_final(eosloc)-1;
gz_final=round(gz*signalscale/kp.gmax)*2;
gz_final(eosloc)=gz_final(eosloc)-1;


if iop.writetofile_sim
	printm('Writing to waveform files for simulation');

	infname_full = [iop.waveformfilespath filesep iop.infname, '.', ...
		num2str(kp.npnts), '.', ...
		int2str(floor(rfp.dfovf)), '.', ...
		int2str(floor(rfp.dfovz)), '.', ...
		num2str(nshots), '.', num2str(pulsen)];

	ir_fwrite([infname_full, '.mag'], magb1_4sim_final, 'short')
	ir_fwrite([infname_full, '.ph'], phb1_final, 'short')
	ir_fwrite([infname_full, '.gx'], gx_final, 'short')
	ir_fwrite([infname_full, '.gy'], gy_final, 'short')
	ir_fwrite([infname_full, '.gz'], gz_final, 'short')
end

if iop.writetofile_scanner

	printm('Writing to waveform files for the scanner...');
	fid=fopen([iop.infname,'.rho'],'w','b');
	fwrite(fid,magb1_final,'short');
	fclose(fid);
	fid=fopen([iop.infname,'.ph'],'w','b');
	fwrite(fid,phb1_final,'short');
	fclose(fid);
	fid=fopen([iop.infname,'.gx'],'w','b');
	fwrite(fid,gx_final,'short');
	fclose(fid);
	fid=fopen([iop.infname,'.gy'],'w','b');
	fwrite(fid,gy_final,'short');
	fclose(fid);
	fid=fopen([iop.infname,'.gz'],'w','b');
	fwrite(fid,gz_final,'short');
	fclose(fid);

	str_transfer = ['mv ',iop.infname,'.* ',iop.waveformfilespath];
	system(str_transfer);
end
