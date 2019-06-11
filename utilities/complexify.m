 function y = complexify(x)
%function y = complexify(x)
% make y complex if x is not.
% matlab's built in "complex" command chokes on complex inputs!

if isreal(x)
	y = complex(x);
else
	y = x;
end
