 function om = nufft_samples(type, Nd)
%function om = nufft_samples(type, Nd)
%|
%| default simple sampling patterns

if nargin < 2, ir_usage, end

switch type

case 'epi' % blipped echo-planar cartesian samples
	if length(Nd) == 1
		om = 2*pi*[-Nd/2:Nd/2-1]';
	elseif length(Nd) == 2
		o1 = 2*pi*[-Nd(1)/2:Nd(1)/2-1]' / Nd(1);
		o2 = 2*pi*[-Nd(2)/2:Nd(2)/2-1]' / Nd(2);
		[o1 o2] = ndgrid(o1, o2);
		o1(:,2:2:end) = flipdim(o1(:,2:2:end),1);
		om = [o1(:) o2(:)];
	else
		fail 'only 1d and 2d "epi" implemented'
	end

otherwise
	fail('unknown type "%s"', type)
end
