function str=indexstr(S)
%Given a subscripting structure array S, such as used in indexing methods
%SUBSREF and SUBSASGN, a single string STR is returned bearing a subscript 
%expression equivalent to that expressed by S. This is useful when combined 
%with EVAL() to implement overloaded SUBSASGN and SUBSREF methods based on 
%already-defined subscript operations.
%
%EXAMPLE:
%
%   M=substruct('{}',{1:2,':'},'.','fld','()',{1:5}); 
%   str=indexstr(M)
%   
%returns
%
% str = {M(1).subs{:}}.fld(M(3).subs{:})
%
%NOTE: When it exists, the INPUTNAME of S is used in generating STR. Otherwise,
%      this name defaults to 'S'
%
%
%See also SUBSTRUCT, INPUTNAME


name=inputname(1);

%DEFAULT
if isempty(name), name='S'; end


%%
str='';

for ii=1:length(S)

 switch S(ii).type

 case '.'



   str=[str '.' S(ii).subs];

 case '()'

   insert=[name '(' num2str(ii) ')'];

   str=[str '(' insert '.subs{:})'];

 case '{}'

   insert=[name '(' num2str(ii) ')'];

   str=[str '{' insert '.subs{:}}'];

 end

end
