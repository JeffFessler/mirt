function patches = img2patch_index(img, patchSize, stride)

m = size(img, 1); n = size(img, 2);
a = patchSize(1); b = patchSize(2); c = stride;

if stride > 1
	aCut = mod(m, a); bCut = mod(n, b);
	img = img(1:end-aCut, 1:end-bCut);
end

% Extended image from the original image
imgExt = [ img img(:, 1:b-c); img(1:a-c, :) img(1:a-c, 1:b-c) ];

if stride == 1
	patches = im2col(imgExt, patchSize, 'sliding');
return
end

iMat = zeros(size(img), class(img));
iMat( 1:c:end, 1:c:end ) = 1; % Take patches in distances of 'stride'
iPatch = find(iMat);
[i, j] = ind2sub(size(img), iPatch);
patches = zeros(prod(patchSize), length(iPatch));
for k = 1:length(iPatch)
	patch = imgExt(i(k):i(k)+a-1, j(k):j(k)+b-1);
	patches(:, k) = patch(:);
end

end % img2patch_index
