% ir_pet_blank_unmash.m
% unmash a 1-mashed PET blank scan
% using another unmashed blank scan as a reference for relative weighting

if ~isvar('b2')
f.mashin = '/users/fessler/data/921,fdg,wbdy/raw/muglia_11f7_6022_bl2.s';
%f.mashin = '/users/fessler/data/921,fdg,wbdy/raw/blank,view.s';

f.recent = '/n/ir6/z/fessler/921,71/whe*/whe*_bl.S';
f.recent = '/n/ir6/z/fessler/921,71/got*/bl.s'

	% fix: this will probably not work anymore!
	% need to convert to .fld file, then use fld_read.m instead!
	b1 = ir_read_op(f.mashin); % fix: use fld_read
	b2 = ir_read_op(f.recent); % fix: use fld_read
	im(squeeze(b1(:,:,1:3)))
end

	[nb,na,nz] = size(b2);

% unmash by using b2 for relative weighting
if ~isvar('bu')
	bu = zeros(size(b2));
	for iz=1:47
		disp(iz)
		for ia=1:na/2
			iu = 2*ia-1;
			r0 = double(b1(:,ia,iz));
			r1 = double(b2(:,iu,iz));
			r2 = double(b2(:,iu+1,iz));
			bu(:,iu  ,iz) = r0 .* r1 ./ (r1 + r2);
			bu(:,iu+1,iz) = r0 .* r2 ./ (r1 + r2);
		end
	end
end

fld_write('raw/blank,unmash.fld', bu)

%	t = bu(:,1:2:end,:) + bu(:,2:2:end,:);	% check at end
%	im([b1(:,:,10); t(:,:,10)])

im(bu(:,:,1:3))
im([3*bu(:,:,11); b2(:,:,11)])
