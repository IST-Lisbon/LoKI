% LoKI-B solves a time and space independent form of the two-term 
% electron Boltzmann equation (EBE), for non-magnetised non-equilibrium 
% low-temperature plasmas excited by DC/HF electric fields from 
% different gases or gas mixtures.
% Copyright (C) 2018 A. Tejero-del-Caz, V. Guerra, D. Goncalves, 
% M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
% L. L. Alves
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function stateArray = rotationalDegeneracy_H2(stateArray, property, ...
  argumentArray)
  % rotationalDegeneracy_H2 (have to be writen)
  
  [stackTrace, ~] = dbstack;
  if ~strcmp(property, 'statisticalWeight')
    error(['Trying to use %s function to set up property %s. Check input '...
      'file'], stackTrace(1).name, property);
  end
  if ~isempty(argumentArray)
    error(['Wrong number of arguments when evaluating %s function. Check '...
      'input file'], stackTrace(1).name)
  end
  for state = stateArray
    if ~strcmp(state.type, 'rot')
      error(['Trying to asign rotational degeneracy to non rotational state '...
        '%s. Check input file', state.name]);
    end
    J = str2double(state.rotLevel);
    state.statisticalWeight = (2-(-1)^J)*(2*J+1);
  end
  
end
