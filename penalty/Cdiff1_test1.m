  function Cdiff1_test1(C)
%|function Cdiff1_test1(C)
%| private helper for Cdiff1_test
%| tests abs and power and adjoints

Cf = C(:,:);

% abs
Ca = abs(C); Caf = Ca(:,:);
jf_equal(Caf, abs(Cf))

% squared
C2 = C.^2; C2f = C2(:,:);
jf_equal(C2f, Cf.^2)

test_adjoint(C);
test_adjoint(Ca);
test_adjoint(C2);
