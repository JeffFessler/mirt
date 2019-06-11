 function out = ir_is_live_script
%function out = ir_is_live_script
%|
%| determine if the calling script is running as Matlab "Live Script"
%| 2016-11-01, Jeff Fessler, University of Michigan

out = false;

if ir_is_octave || isfreemat || ~exist('dbstack', 'builtin')
	return
end

% code based on caller_name
st = dbstack;
tmp = {st(:).name};
tmp = strcmp(tmp,  'LiveEditorEvaluationHelper');
out = any(tmp);

% EvaluationOutputsService.evalRegions is top-level caller
