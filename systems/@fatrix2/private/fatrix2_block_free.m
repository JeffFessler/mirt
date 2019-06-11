function fatrix2_block_free(arg)
%printm 'freeing fatrix2_block object static memory'
for mm=1:length(arg.blocks)
	try
		free(blocks{mm})
	catch
	end
end

