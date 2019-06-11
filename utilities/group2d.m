  function groups = group2d(nx, ny, ngv, mask, chat)
%|function groups = group2d(nx, ny, ngv, mask, chat)
%| build groups for GCA for 2d case
%| ngv [2]	[ngx ngy]
%|
%| Copyright Mar. 1999, Jeff Fessler

if 1 ~= exist('chat'), chat = 1; end
if 1 ~= exist('mask') || isequal(mask, 1)
	mask = ones(nx,ny);
end

	ngx = ngv(1);
	ngy = ngv(2);

	np = nx*ny;

	ng = ngx*ngy;
	groups	= logical(zeros(np, ng));
	for igx = 1:ngx
		for igy = 1:ngy
			ig = igx + (igy-1)*ngx;
			gx = [igx:ngx:nx]';
			gy = [igy:ngy:ny]';
			ggx = gx * ones(1,length(gy));
			ggy = ones(length(gx),1) * gy';
			gg = ggx + (ggy-1) * nx;
			gg = gg(:);
			groups(gg,ig) = ones(size(gg));
		end
	end

	groups = groups(find(mask(:)),:);

	% look at groups
	if chat
		im(groups', 'groups')
		for ig=1:ng
			prompt
			im(embed(groups(:,ig), mask), sprintf('group %d', ig))
		end
		prompt
		im(embed(groups * [1:ng]',mask), 'All Groups')
	end

	if any(sum(groups')) > 1, error 'overlap groups', end
