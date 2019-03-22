function stateArray = functionHandle(stateArray, property, argumentArray)
  % functionHandle (have to be writen)
  
  if length(argumentArray) ~= 1
    [stackTrace, ~] = dbstack;
    error('Wrong number of arguments when evaluating %s function. Check input file', stackTrace(1).name);
  end
  for state = stateArray
    state.(property) = argumentArray;
  end
  
end
