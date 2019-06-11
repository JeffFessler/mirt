function print_path(p)

% print all but the last dir
dirs = regexp(p, '[^:]*:', 'match');
for ii=1:length(dirs)
  fprintf('%s\n',dirs{ii});
end

%print the last dir
lastdir = regexp(p, '[^:]*$', 'match');
if length(lastdir)~=0
  fprintf('%s\n',lastdir{1});
end
