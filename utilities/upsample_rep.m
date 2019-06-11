 function y = upsample_rep(x, m)
%function y = upsample_rep(x, m)
%|
%| upsample a 2D image a factor of m by simple replication
%| 2018-06-29, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), upsample_rep_test, return, end
if nargin == 1 && streq(x, 'time'), upsample_rep_time, return, end
if nargin < 1, ir_usage, end
if nargin < 2, m = [2 2]; end

if numel(m) > 2, fail 'only 2D done', end

y = upsample2_rep_loop(x, m);


% 1d upsampling of each column using a loop
function y = upsample1_rep_loop(x, m)
[n1, n2] = size(x);
y = zeros(m*n1,n2);
for ii=1:m
	y(ii+m*[0:n1-1],:) = x;
end


% 1d upsampling along first dimension using repmat 
function y = upsample1_rep_repmat(x, m)
M = size(x, 1);
y = x(:); % [MN 1]
y = repmat(y, [1 m]); % [MN m]
y = y.'; % [m MN]
y = reshape(y, m*M, []); % [mM N]


% 2d
function y = upsample2_rep_loop(x, m)
y = upsample1_rep_loop(x, m(1));
y = upsample1_rep_loop(y.', m(end)).';


% 2d
function y = upsample2_rep_repmat(x, m)
y = upsample1_rep_repmat(x, m(1));
y = upsample1_rep_repmat(y.', m(end)).';


% test
function upsample_rep_test
x = reshape(1:24, [4 6]);
if 1
	y1 = upsample2_rep_loop(x, 2);
	y2 = upsample2_rep_repmat(x, 2);
	jf_equal(y1, y2)
	y1 = upsample2_rep_loop(x, [2 3]);
	y2 = upsample2_rep_repmat(x, [2 3]);
	jf_equal(y1, y2)
end
y = upsample_rep(x, 2);
for i1=1:2
	for i2=1:2
		jf_equal(x, y(i1:2:end,i2:2:end))
	end
end


% time
function upsample_rep_time
M = 2^11;
N = 2^11+2;
x = reshape(1:M*N, M, N);
m = 2;
upsample2_rep_loop(ones(m), m); % warm-up
upsample2_rep_repmat(ones(m), m); % warm-up
cpu etic
y1 = upsample2_rep_loop(x, m);
cpu etoc :loop
cpu etic
y2 = upsample2_rep_repmat(x, m);
cpu etoc :repmat
jf_equal(y1, y2)
