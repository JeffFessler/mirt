 function yy = ir_imfill1(xx)
%function yy = ir_imfill1(xx)
%|
%| 1D version of imfill 'holes' that works along 1st dimension of input
%|
%| in
%| xx [N (L)]	logical array
%| out
%| yy [N (L)]	logical array
%|
%| 2015-08-13, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end
if streq(xx, 'test'), ir_imfill1_test, return, end

if ~islogical(xx), fail('input must be logical'), end

yy = cummax(xx,1) & flipdim(cummax(flipdim(xx,1),1),1);


function ir_imfill1_test
x = ellipse_im(64) > 0;
x(end/2,end/4+[1:5]) = 0;
y1 = ir_imfill1(x);
y2 = ir_imfill1(x')';
im plc 2 2
im(1, x)
im(2, y1)
im(3, y2)
im(4, y2-y1)
