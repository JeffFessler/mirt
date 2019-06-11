  function yn = isfield(st, name)
%|function yn = isfield(st, name)

yn = isfield(st.data, name) || isfield(st.meth, name);
