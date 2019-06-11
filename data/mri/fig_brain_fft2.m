% fig_brain_fft2
x = imread('brainweb_t1.gif')';
f = fftshift(fft2(x));

im plc 3 4
im(1, x)
axis off, title ''
%colormap(gca, hot)
tmp = abs(f); tmp = tmp / max(tmp(:));
im(2, log(0.0001 + tmp))
axis off, title ''
cmap = colormap(jet);
cmap(1,:) = 0;
colormap(cmap)

im subplot 3
plot(200*tmp(:,end/2+[-10 1 10]), 'r')
axis([0 256 0 200])
xtick([1 256])
ytick([0 200])
axis off
