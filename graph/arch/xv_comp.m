function xv_comp(i1, i2, clim)
% compare two images using xv to toggle back and forth

i1 = (i1 - clim(1)) / (clim(2) - clim(1)) * 255;
i2 = (i2 - clim(1)) / (clim(2) - clim(1)) * 255;
clim = [0 255];

i1 = max(i1, clim(1));
i1 = min(i1, clim(2));
i2 = max(i2, clim(1));
i2 = min(i2, clim(2));
i1 = uint8(i1);
i2 = uint8(i2);

im([i1; i2])
imwrite(i1, 'tmp1.pgm')
imwrite(i2, 'tmp2.pgm')
os_run('xv tmp1.pgm tmp2.pgm')
