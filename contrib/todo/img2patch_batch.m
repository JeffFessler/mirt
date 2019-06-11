function patches = img2patch_batch(img, patchSize, stride)

a = patchSize(1); b = patchSize(2); c = stride;

if stride > 1
	cut = mod(size(img), patchSize);
	img = img(1:end-cut(1), 1:end-cut(2));
end

nPatches = prod(size(img)/c);

% Extended image from the original image
img = [ img img(:, 1:b-c); img(1:a-c, :) img(1:a-c, 1:b-c) ];
imgSize = size(img);

% Linear indices for the 1st patch
stride = [1 imgSize(1)];
limit = (patchSize-1).*stride;
ind = 1;
for dim = 1:numel(imgSize)
	ind = bsxfun (@plus, ind(:), 0:stride(dim):limit(dim));
end

% Linear indices for all patches
slides = imgSize - patchSize;
limit = slides.*stride;
stride = stride*c;
for dim = 1:numel(imgSize)
	ind = bsxfun( @plus, ind(:), 0:stride(dim):limit(dim) );
end
patches = reshape(img(ind(:)), prod(patchSize), nPatches);

end % img2patch_batch
