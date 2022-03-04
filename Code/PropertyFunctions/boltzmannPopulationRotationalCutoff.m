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

function population = boltzmannPopulationRotationalCutoff(state, argumentArray, workCond)
  % boltzmannPopulationRotationalCutoff is a property function that evaluates a Boltzmann distribution with a certain
  % temperature (specified as the first argument provided to the function in the setup file) and up to a certain 
  % rotational level (specified as the second argument provided to the function in the setup file). This function
  % should be used ONLY to provide the population of rotational distributions functions, for other distributions check
  % boltzmannPopulation and boltzmannPopulationVibrationalCutoff.
  
  % obtain temperature of the distribution (either prescribed, i.e. numeric, or found in the working conditions)
  temperature = argumentArray{1};
  if ~isnumeric(temperature)
    switch temperature
      case 'gasTemperature'
        temperature = workCond.gasTemperature;
      case 'electronTemperature'
        temperature = workCond.electronTemperature/Constant.boltzmannInEV;
      otherwise
        error(['Error found when evaluating population of state %s.\nTemperature ''%s'' not defined in the ' ...
          'working conditions.\nPlease, fix the problem and run the code again.'], state.name, temperature);
    end
  end
  
  % obtain the cutoff rotational level for the Boltzmann distribution (this level is included in the distribution)
  jMax = argumentArray{2};
  if ~isnumeric(jMax)
    error(['Error found when evaluating population of state %s.\nMaximum rotational level ''%s'' should be a ' ...
      'numerical value.\nPlease, fix the problem and run the code again.'], state.name, argumentArray{2});
  end
  
  % error checking (find if state is a rotational level)
  if ~strcmp(state.type, 'rot')
    error(['Trying to asign ''boltzmannPopulationRotationalCutoff'' population to non rotational state %s. ' ...
      'Check input file'], state.name);
  end
  
  % evaluate Boltzmann distribution for the state and its siblings
  norm = 0;
  for stateAux = [state state.siblingArray]
    if isempty(stateAux.energy)
      error(['Unable to find %s energy for the evaluation of ''boltzmannPopulationRotationalCutoff'' function.\n'...
        'Check input file'], stateAux.name);
    elseif isempty(stateAux.statisticalWeight)
      error(['Unable to find %s statistical weight for the evaluation of ''boltzmannPopulationRotationalCutoff'' '...
        'function.\nCheck input file'], stateAux.name);
    end
    if str2double(stateAux.rotLevel) > jMax
      stateAux.population = 0;
    else
      stateAux.population = stateAux.statisticalWeight*exp(-stateAux.energy/(Constant.boltzmannInEV*temperature));
      norm = norm + stateAux.population;
    end
  end
  for stateAux = [state state.siblingArray]
    stateAux.population = stateAux.population/norm;
  end
  
  % return population of the current state
  population = state.population;
  
end
