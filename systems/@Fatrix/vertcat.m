 function c = vertcat(varargin)
%function c = vertcat(varargin)
% called for c = [a; b; ...] where any of them is a fatrix

c = block_fatrix(varargin, 'type', 'col');
