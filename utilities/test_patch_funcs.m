
%%
n=80;
test_im = cat(3,randn(n,n),phantom(n,n)*3,phantom(n,n)*8);
stride = 5;
for ps = 16:24
    test_im_patches = createpatches(test_im,ps,stride);
    sew_test = sewpatches(test_im_patches,[n,n],stride);
%     vuOnePaneViewer(cat(2,test_im,sew_test)); 
    figure; imagesc(cat(2,test_im(:,:,1),sew_test(:,:,1)))
end
