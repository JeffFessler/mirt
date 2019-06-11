  function jf_dcm_write(data, dcm_name, WindowCenter, WindowWidth, varargin)
%|function jf_dcm_write(data, dcm_name, WindowCenter, WindowWidth, varargin)
%| write DICOM (.dcm) filr
%| for 2d CT images only
%|
%| based on fld2dcm by Somesh Srivastava 2007-12-19

if (nargin == 1 && streq(data, 'test')) % run test routine
	data = single(1000 + 30*peaks(512));
	im plc 1 2
	im(1, data, [800 1200])
	% a dcm viewer should show an image similar to the one displayed above
	tmp = '/tmp/peaks.dcm';
	jf_dcm_write(data, tmp, 1000, 400)
	tmp = dicomread(tmp)';
	im(2, tmp)
return
end

if (nargin < 4), ir_usage, return, end

% the last 4 options below are required, but seem not to be important
% for display purposes; they can be anything for now.
arg.RescaleIntercept = 0;
arg.RescaleSlope = 1;
arg.PixelSpacing = [1 1] * 700.3125/512;
arg.ImageOrientationPatient = [1 0 0 0 1 0];
arg.ImagePositionPatient = [-250.0000 -250.0000 -57.8750];
arg.ImageType = 'ORIGINAL\SECONDARY\AXIAL';
arg = vararg_pair(arg, varargin);

tmp_file = '/tmp/fld2dcm_temp.dcm';

data = int16(data);
data = permute(data, [2 1]); % transpose in 2d

% write out a dummy dcm file. this lacks required fields.
status = dicomwrite(data, tmp_file, 'ObjectType', 'CT Image Storage');

if (length(status.MissingData) ~= 0)
%	printm 'inserting missing metadata ...'

	metadata = dicominfo(tmp_file);
	delete(tmp_file);

	metadata.WindowCenter	= WindowCenter;
	metadata.WindowWidth	= WindowWidth;
	metadata.RescaleIntercept = arg.RescaleIntercept;
	metadata.RescaleSlope	= arg.RescaleSlope;
	metadata.PixelSpacing	= arg.PixelSpacing;
	metadata.ImageOrientationPatient = arg.ImageOrientationPatient;
	metadata.ImagePositionPatient	= arg.ImagePositionPatient;
	metadata.ImageType	= arg.ImageType;

	status = dicomwrite(data, dcm_name, metadata, ...
		'ObjectType', 'CT Image Storage');
	if (length(status.MissingData) ~= 0)
		fail('failed. quitting')
	end
%	printm 'done'
else
	delete(tmp_file)
	fail('succeded unexpectedly. quitting')
end
