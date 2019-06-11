function pout=delete_path(pin, paths)

pout=pin;
for ii=1:length(paths)
  pat = ['[^:]*' paths{ii} '[^:]*:'];
  pout = regexprep(pout, pat, '');
  pat = ['[^:]*' paths{ii} '[^:]*$'];
  pout = regexprep(pout, pat, '');
end

