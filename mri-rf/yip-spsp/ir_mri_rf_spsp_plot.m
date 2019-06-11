 function ir_mri_rf_spsp_plot(b, gz, d, f, z, mm, t)
%function ir_mri_rf_spsp_plot(b, gz, d, f, z, mm, t)
%|
%| plot SPSP pulses etc.
%| code extracted from compute_rf_spsp_mgh

%printm('Displaying SPSP pulse design results...')

pl = @(i) subplot(340+i);

pl(1)
plot(t*1000,gz)
axis tight
grid
xlabel('Time (ms)')
ylabel('g/cm')
title('z gradient')

pl(5)
plot(t*1000,real(b),'b-', t*1000,imag(b),'r-')
axis tight
grid
xlabel('Time (ms)')
ylabel('g')
title('SPSP pulse')
legend('Real','Imaginary')

%Display desired pattern of the SPSP recovery pulse.
pl(2)
imagesc(f,z,abs(d))
colormap default;colorbar
xlabel('Frequency offset (Hz)')
ylabel('z (cm)')
title('Desired SPSP pattern (magnitude)')

pl(6)
imagesc(f,z,angle(d).*abs(d),[-pi,pi]);colorbar
colormap(gca, hsv)
xlabel('Frequency offset (Hz)')
ylabel('z (cm)')
title('Desired SPSP pattern (phase)')

pl(3)
imagesc(f,z,abs(mm))
xlabel('frequency (Hz)')
ylabel('z (cm)')
title('Resulting pattern (linear regime prediction)')

pl(7)
imagesc(f,z,angle(mm))
colormap(gca, hsv)
xlabel('frequency (Hz)')
ylabel('z (cm)')
title('Resulting phase pattern (linear regime prediction)')
