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

function stateArray = treanorPopulation(stateArray, property, argumentArray)
  % treanor (have to be writen)
  
  [stackTrace, ~] = dbstack;
  if ~strcmp(property, 'population')
    error(['Trying to use %s function to set up property %s. Check input '...
      'file'], stackTrace(1).name, property);
  end
  if length(argumentArray) ~= 2
    error(['Wrong number of arguments when evaluating %s function. Check '...
      'input file'], stackTrace(1).name)
  end
  temp0 = argumentArray(1);
  temp1 = argumentArray(2);
  norm = 0;
  groundEnergy = [];
  firstEnergy = [];
  for state = stateArray
    if ~strcmp(state.type, 'vib')
      error(['Trying to asign treanor population to non vibrational state '...
        '%s. Check input file', state.name]);
    end
    if strcmp(state.vibLevel, '0')
      groundEnergy = state.energy;
    elseif strcmp(state.vibLevel, '1')
      firstEnergy = state.energy;
    end
  end
  if isempty(groundEnergy) || isempty(firstEnergy)
    error(['Unable to find groundEnergy or firstEnergy to populate state %s'...
      ' and its siblings with function %s.\nCheck input file'], ...
      stateArray(end).name, stackTrace(1).name);
  end
  for state = stateArray
    if isempty(state.energy)
      error(['Unable to find %s energy for the evaluation of %s function.\n'...
        'Check input file'], state.name, stackTrace(1).name);
    elseif isempty(state.statisticalWeight)
      error(['Unable to find %s statistical weight for the evaluation of %s '...
        'function.\nCheck input file'], state.name, stackTrace(1).name);
    end
    vibLevel = str2double(state.vibLevel);
    state.population = state.statisticalWeight*exp(...
      -(1./Constant.boltzmannInEV)*(vibLevel*(firstEnergy-groundEnergy)*...
      (1/temp1-1/temp0)+(state.energy-groundEnergy)/temp0));
    norm = norm + state.population;
  end
  for state = stateArray
    state.population = state.population/norm;
  end
  
end
