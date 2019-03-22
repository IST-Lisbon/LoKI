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

classdef Gas < handle
  %Gas Class that stores the information of a certain gas
  %   Class that stores the information of a certain gas present in the
  %   mixture of gases included in the discharge. The electron collisions
  %   and the populations of the different states are to be used in a
  %   Boltzmann solver to obtain the EEDF.

  properties

    ID = -1;                          % ID that identifies the gas in the gas array
    name = '';                        % name of the gas
    mass = [];                        % mass of the gas
    harmonicFrequency = [];           % harmonic oscillator frequency (molecules)
    anharmonicFrequency = [];         % anharmonic oscillator frequency (molecules)
    rotationalConstant = [];          % rotational constant (molecules)
    lennardJonesDistance = [];        % sigma parameter of Lennard-Jones potential
    lennardJonesDepth = [];           % epsilon parameter of Lennard-Jones potential
    electricDipolarMoment = [];       % electric dipolar moment (molecules)
    electricQuadrupoleMoment = [];    % electric cuadrupole moment (molecules)
    polarizability = [];              % polarizability (molecules)
    fraction = 0;                     % fraction in the gas mixture
    stateArray;                       % State array with all the states of the gas (initialized in the subclass)

  end

  methods

    function disp(gas)
      for field = fieldnames(gas)'
        if ~isempty(gas.(field{1}))
          if isnumeric(gas.(field{1}))
            fprintf('%s: ', field{1});
            fprintf('%g ', gas.(field{1}));
            fprintf('\n');
          elseif ischar(gas.(field{1}))
            fprintf('%s: %s\n', field{1}, gas.(field{1}));
          end
        end
      end
    end

    function dispStates(gas)
      fprintf('%s states:\n', gas.name);
      for i = 1:length(gas.stateArray)
        fprintf('\t %d.- %s\n', i, gas.stateArray(i).name);
      end
    end

    function dispFamilyTree(gas)
      fprintf('%s states: (family tree)\n', gas.name);
      fprintf('%s\n', gas.name);
      for eleState = gas.stateArray
        if strcmp(eleState.type, 'ele')
          fprintf('\t|%s\n', eleState.name);
          for vibState = eleState.childArray
            fprintf('\t|\t|%s\n', vibState.name);
            for rotState = vibState.childArray
              fprintf('\t|\t|\t|%s\n', rotState.name);
            end
          end
        end
      end
      for ionState = gas.stateArray
        if strcmp(ionState.type, 'ion')
          fprintf('%s\n', ionState.name);
        end
      end
    end

    function checkPopulationNorms(gas)
    % checkPopulationNorms checks for the population of the different states of the gas to be properly normalised,
    % i. e. the populations of all sibling states should add to one.

      % avoid dummy gases
      if gas.fraction == 0
        return;
      end

      % check norm of electronic/ionic states
      gasNorm = 0;
      electronicStatesToBeChecked = true;
      ionicStatesToBeChecked = true;
      for state = gas.stateArray
        % norm of electronic states
        if strcmp(state.type, 'ele') && electronicStatesToBeChecked
          for eleState = [state state.siblingArray]
            if eleState.population ~= 0
              gasNorm = gasNorm + eleState.population;
              % check norm of vibrational states (if they exist)
              if ~isempty(eleState.childArray)
                vibNorm = 0;
                for vibState = eleState.childArray
                  if vibState.population ~= 0
                    vibNorm = vibNorm + vibState.population;
                    % check norm of rotational states (if they exist)
                    if ~isempty(vibState.childArray)
                      rotNorm = 0;
                      for rotState = vibState.childArray
                        if rotState.population ~= 0
                          rotNorm = rotNorm + rotState.population;
                        end
                      end
                      if abs(rotNorm-1) > 10*eps(1)
                        stateDistName = vibState.name;
                        stateDistName = [stateDistName(1:end-1) ',J=*)'];
                        error('Rotational distribution %s is not properly normalised. (Error = %e)\n', ...
                          stateDistName, rotNorm-1);
                      end
                    end
                  end
                end
                if abs(vibNorm-1) > 10*eps(1)
                  stateDistName = eleState.name;
                  stateDistName = [stateDistName(1:end-1) ',v=*)'];
                  error('Vibrational distribution %s is not properly normalised. (Error = %e)\n', ...
                    stateDistName, vibNorm-1);
                end
              end
            end
          end
          electronicStatesToBeChecked = false;
        end
        if strcmp(state.type, 'ion') && ionicStatesToBeChecked
          for ionState = [state state.siblingArray]
            if ionState.population ~= 0
              gasNorm = gasNorm + ionState.population;
            end
          end
          ionicStatesToBeChecked = false;
        end
      end
      if abs(gasNorm-1) > 10*eps(1)
        stateDistName = [gas.name '(*)'];
        error('Electronic/ionic distribution %s is not properly normalised. (Error = %e)\n', stateDistName, gasNorm-1);
      end

    end

    function renormalizeWithDensities(gas)
    % renormalizeWithDensities is a function that obtain new "populations" of the states of a gas once their
    % "densities" have been updated. Note that the "populations" of the states are just the "densities" normalized
    % to the sibling states and that the "densities" are the absolute densites (in m^-3) normalized to the total
    % gas density.

      % loop over electronic/ionic states to obtain the new gas fraction (sum of electronic/ionic densities)
      gasFraction = 0;
      eleStatesNeedToBeChecked = true;
      ionStatesNeedToBeChecked = true;
      for state = gas.stateArray
        % check electronic states
        if strcmp(state.type, 'ele') && eleStatesNeedToBeChecked
          for eleState = [state.siblingArray state]
            gasFraction = gasFraction + eleState.density;
          end
          eleStatesNeedToBeChecked = false;
        end
        % check ionic states
        if strcmp(state.type, 'ion') && ionStatesNeedToBeChecked
          for ionState = [state.siblingArray state]
            gasFraction = gasFraction + ionState.density;
          end
          ionStatesNeedToBeChecked = false;
        end
        % if everything is checked exit the loop and save the new gas fraction
        if ~eleStatesNeedToBeChecked && ~ionStatesNeedToBeChecked
          gas.fraction = gasFraction;
          break;
        end
      end

      % loop over electronic and ionic states to obtain the new populations (for states with a chemistry equivalent)
      eleStatesNeedToBeChecked = true;
      ionStatesNeedToBeChecked = true;
      for state = gas.stateArray
        % check electronic states
        if strcmp(state.type, 'ele') && eleStatesNeedToBeChecked
          for eleState = [state.siblingArray state]
            if isempty(eleState.chemEquivalent)
              continue
            end
            eleState.population = eleState.density/gasFraction;
            % obtain the new populations of the vibrational distribution (in case it exists)
            if ~isempty(eleState.childArray)
              for vibState = eleState.childArray
                if isempty(vibState.chemEquivalent)
                  continue
                end
                vibState.population = vibState.density/eleState.density;
                % obtain the new populations of the rotational distribution (in case it exists)
                if ~isempty(vibState.childArray)
                  for rotState = vibState.childArray
                    if isempty(rotState.chemEquivalent)
                      continue
                    end
                    rotState.population = rotState.density/vibState.density;
                  end
                end
              end
            end
          end
          eleStatesNeedToBeChecked = false;
        end
        % check ionic states
        if strcmp(state.type, 'ion') && ionStatesNeedToBeChecked
          for ionState = [state.siblingArray state]
            if isempty(ionState.chemEquivalent)
              continue
            end
            ionState.population = ionState.density/gasFraction;
          end
          ionStatesNeedToBeChecked = false;
        end
        % exit the loop if everything is checked
        if ~eleStatesNeedToBeChecked && ~ionStatesNeedToBeChecked
          break;
        end
      end

      % evaluate density of states with the updated populations
      gas.evaluateDensities();

    end
    
    function evaluateDensities(gas)
    % evaluateDensities is a function that evaluates the densities of all the states of the gas from their
    % populations and the gas fraction. Note that the densities of the states are normalized to the total gas
    % density.

      % loop over all the states of the gas
      for state = gas.stateArray
        state.evaluateDensity();
      end

    end

  end

  methods (Static)

    function [gasArray, gasID] = add(gasName, gasArray)
      % addGas looks for gasName in Gas array gas. If a Gas with the same
      % name is found the function returns its ID. If no Gas is found with
      % the same name, a new Gas is added to the array and its ID is returned.

      gasID = Gas.find(gasName, gasArray);
      if gasID ~= -1
        return
      end
      switch class(gasArray)
        case 'EedfGas'
          gasArray(end+1) = EedfGas(gasName);
        case 'ChemGas'
          gasArray(end+1) = ChemGas(gasName);
        otherwise
          error('Gas class ''%s'' not recognised when trying to create gas %s\n.', class(gas), gas.name);
      end
      gasID = gasArray(end).ID;

    end

    function gasID = find(gasName, gasArray)

      for gas = gasArray
        if strcmp(gasName, gas.name)
          gasID = gas.ID;
          return
        end
      end
      gasID = -1;

    end

    function checkFractionNorm(gasArray)
    % checkFractionsNorms checks for the fractions of the different gases
    % in gasArray to add to one. In case the fractions of the gases are not
    % properly normalised, an error message is thrown.
    %
    % Note: the function avoid dummy gases (gases created for the sake of a
    % pretty output)

      norm = 0;
      for gas = gasArray
        if isempty(gas.fraction)
          continue;
        end
        norm = norm + gas.fraction;
      end
      if abs(norm-1) > 10*eps(1)
        error(['Gases fractions are not properly normalised (Error = %e).\n'...
          'Please, check input file.'], norm-1);
      end

    end

  end

end
