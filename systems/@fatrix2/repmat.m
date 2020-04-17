function ob = repmat(A, dims)
%|function ob = repmat(A, dims)
%|
%| called for ob = repmat(A, dims) where A is a fatrix to be
%| "repmatted" along dimensions dim, like the MATALB version
%|  of the function
%|
%| Inputs:
%|  A:      input fatrix
%|  dims:   number of replicated along each dimension
%|
%| Outputs:
%|  ob:     replicated fatrix of A
%|
%| Example:
%|   repmat(A,[3 1]) = [A; A; A]
%|
%| numel(dims) > 2 is not supported

if nargin < 2, ir_usage, end

if numel(dims) > 2 || (prod(dims) > max(dims))
    fail('unsupported for 2d repmat, change to dims = [N 1] or [1 N]');
end

if dims(1) > 1  % vertical cat
    N = dims(1);
    B = cell(N,1);
    for nn = 1:N
        B{nn} = A;
    end
    
    ob = cat(1,B{:});
else  % horizontal cat
    N = dims(2);
    B = cell(1,N);
    for nn = 1:N
        B{nn} = A;
    end
    
    ob = cat(2,B{:});
end

end
