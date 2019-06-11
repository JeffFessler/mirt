% runs after parallel_feldkamp_example to create images

close all
im plc 2 2

im(1,xtrue(:,:,:),[0,.02],'True')
im(2,tfdk(:,:,:),[0,.02],'T-FDK')
im(3,pfdk(:,:,:),[0,.02],'P-FDK')
im(4,xfdk(:,:,:),[0,.02],'FDK')
