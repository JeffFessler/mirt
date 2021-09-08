%% test_patch_funcs.m
%  function to test patch creation and recombination (or sewing) with
%  various image, patch, and step sizes
%
%  Melissa Haskell, University of Michigan, 2021-09-08

%% set sizes to test (good to have some even and odd of each)
n_all = 85:88;
stride = 1:5;
ps = 9:15;

nn = numel(n_all); ns = numel(stride); nps = numel(ps);
ntestcases = nn*ns*nps; ind = 0;

%% loop through various image, patch, and step sizes to test
fprintf('Testing patch creation and recombination '); nbytes = 0;
for ii = 1:nn
    for jj = 1:ns
        for kk = 1:nps
            
            ind = ind+1;
            fprintf(repmat('\b',1,nbytes))
            nbytes = fprintf('%d of %d', ind, ntestcases);

            n = n_all(ii);
            s = stride(jj);
            patch_sz = ps(kk);
            
            % create 3d test image and pad
            test_im = cat(3,phantom(n,n),phantom(n,n)*3,phantom(n,n)*8);
            test_im = padarray(test_im,[patch_sz,patch_sz],0,'both');
            n2 = size(test_im,1);
            
            % create patches and sew back together
            test_im_patches = createpatches(test_im,patch_sz,s);
            sew_test = sewpatches(test_im_patches,[n2,n2],s);
            
            % see if there is any significant error
            tmperr = norm(test_im(:)-sew_test(:))/norm(test_im(:));
            if tmperr > 1e-10
                warning(['patch test failed for image size ', num2str(n2),...
                    ', patch size ', num2str(patch_sz), ...
                    ', and stride ', num2str(s)])
            end
        end
    end
end
disp(' ')
disp('Patch test complete.')




