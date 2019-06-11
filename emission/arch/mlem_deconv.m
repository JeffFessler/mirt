From fessler@eecs.umich.edu Thu Jun  1 16:56 EDT 2000
To: millsk@engin.umich.edu
Subject: ml-em
Cc: fessler@eecs.umich.edu
Content-Type: text
Content-Length: 539


function f = mlem(f, g, h)

if nargin < 3
	ftrue = zeros(11); ftrue(3:5,4:7) = 1; ftrue(4,5) = 0;
	h = ones(3,3);
	g = conv2(ftrue, h, 'same');
	g = g + 0.1 * randn(size(g));
	g = max(g, 0);

	f = g;
	imagesc(f), drawnow
	for ii=1:100
		f = mlem(f, g, h);
		imagesc(f), drawnow
	end
return
end

if any(g < 0), error 'need nonnegative data', end
gp = conv2(f, h, 'same');
if any(g > 0 & ~gp), error 'model mismatch', end
ratio = g ./ (gp + eps * (gp==0));

sens = conv2(ones(size(f)), h, 'same');
f = f .* conv2(ratio, h, 'same') ./ sens;

