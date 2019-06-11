% ir_read_pet_trans_2d_sino_phantom1
% read 2D PET transmission scan data of phantom from raw files

fid1 = fopen('phan,trans.raw', 'r', 'ieee-be');
f1 = fread(fid1, [160 192], 'int16'); 
fclose(fid1);

im(f1)
