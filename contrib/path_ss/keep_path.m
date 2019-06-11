function pout=keep_path(pin, paths)

pout='';
for ii=1:length(paths)

  % try the all but the last path directory
  pat = ['[^:]*' paths{ii} '[^:]*:'];
  pos = regexp(pin, pat, 'match');
  length(pos)
  if length(pos)~=0
    pos1 = strcat(pos{:});
    pout = strcat(pout, pos1);
  end

  % try the last path directory
  pat = ['[^:]*' paths{ii} '[^:]*$'];
  pos = regexp(pin, pat, 'match');
  length(pos)
  if length(pos)~=0
    pos1 = strcat(pos{:});
    pout = strcat(pout, pos1);
  end
 
end
pout = regexprep(pout,':$','');
