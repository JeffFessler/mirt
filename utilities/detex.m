 function y = detex(x);
%function y = detex(x);
% hide all special characters in a string so that tex interpretation
% will not mess it up, like underscores.
% not finished...
y = strrep(x, '_', '\_');
