  function out = full(ob, varargin)
% function out = full(ob, varargin)
%| full(ob) = ob(:,:)
%|
%| note: alternative is ob * eye(np) but this may yield the
%| wrong precision because eye() is double but ob may be single.
%|
%| optional argument(s) passed to fatrix2_subsref_colon()
%| to facilitate test_adjoint
%|
%| Copyright 2010-12-05, Jeff Fessler, University of Michigan

out = fatrix2_subsref_colon(ob, ':', varargin{:}); % ob(:,:)
