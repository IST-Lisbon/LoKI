% LoKI-GM comprises two modules, that can run self-consistently coupled 
% or as standalone tools.
% 
% LoKI-B, which solves the space independent form of the two-term 
% electron Boltzmann equation (EBE) to calculate the isotropic and the 
% anisotropic parts of the electron distribution function, 
% and the associated electron macroscopic parameters. 
% LoKI-B applies to non-magnetised non-equilibrium LTPs, excited by 
% DC/HF electric fields or time-dependent (non-oscillatory) electric fields 
% from different gases or gas mixtures. 
% The tool uses a stationary description for DC fields, 
% a Fourier time-expansion description for HF fields, 
% and a time-dependent description for time-varying fields.
% 
% LoKI-C, which solves the system of zero-dimensional (volume average) 
% rate balance equations for the most relevant 
% charged and neutral species in the plasma. 
% LoKI-C receives as input data the kinetic schemes for the gas/plasma/
% surface system under study, via an intuitive csv-like input file, 
% and gives as output the particle densities of the different gas/plasma/
% surface species, the corresponding creation/destruction reaction rates, 
% and the reduced electric field (and any related quantity, such as 
% the discharge current or the discharge power-density).
% The tool uses several modules to describe the mechanisms 
% (collisional, radiative and transport) controlling the
% creation/destruction of species, namely various transport models 
% for the charged particles and for the neutral particles. 
% LoKI-C includes also a gas/plasma thermal model, for the self-consistent 
% calculation of the gas temperature, and supports multicomponent 
% mean-field microkinetic mesoscopic models to handle surface kinetics 
% in a fully coupled way with volume kinetics.
%
% Copyright (C) 2018 A. Tejero-del-Caz, V. Guerra, D. Goncalves, 
% M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
% L. L. Alves
% 
% Copyright (C) 2026 L. L. Alves, A. Tejero-del-Caz, T. C. Dias, 
% A. Gonçalves, L. Marques, P. Pereira, N. Pinhão, C. D. Pintassilgo, 
% T. Silva, P. Viegas and V. Guerra
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

classdef Grid < handle
  %Grid Class that defines a grid and its interactions with other objects
  %   Class that defines a grid and its interactions with other objects. A
  %   grid object stores three properties: the values of the grid at node 
  %   positions, the values of the grid at cell positions (middle point) and 
  %   the step of the grid. 
  %
  %   Node values -> |     |     |     |     |     |     |     |
  %   Grid        -> o-----o-----o-----o-----o-----o-----o-----o
  %   Cell values ->    |     |     |     |     |     |     |
  
  properties

    maxEnergy = [];         % maximum energy of the grid (in eV)
    cellNumber = [];        % number of cells in the energy grid
    node = [];              % values of the grid at node positions
    cell = [];              % values of the grid at cell positions (i.e. between two consecutive nodes)
    energyStep = [];        % values of the energy steps, as difference between the energies at consecutive nodes 
    step = [];              % value of the first energy step (for backward compatibility with uniform energy grid)

    variableGrid = [];      % features of variable energy grid
    inelCollisions = [];    % activation of inelastic collisions - TO BE DELETED
    ratio = 1.1;            % ratio 'a' for geometric progression (for variable energy grid)
    firstEnergyStep = [];   % first energy step (for variable energy grid)

    isSmart = false;        % smart properties of the energy grid (deactivated by default)
    minEedfDecay = [];      % minimum number of decades of decay for the EEDF
    maxEedfDecay = [];      % maximum number of decades of decay for the EEDF
    updateFactor = [];      % percentage factor to update the maximum value of the energy grid
    
  end
  
  events
    updatedMaxEnergy1
    updatedMaxEnergy2
  end

  methods

    function grid = Grid(gridProperties)
      grid.cellNumber = gridProperties.cellNumber;
      grid.maxEnergy = gridProperties.maxEnergy;

      % Check if variableGrid and inelCollisions are properties
      if isfield(gridProperties, 'variableGrid')
          grid.variableGrid = gridProperties.variableGrid;
      end
      if isfield(gridProperties, 'inelCollisions')
          grid.inelCollisions = gridProperties.inelCollisions;
      end

      % Configure the energy grid
      if grid.variableGrid

          % If variableGrid is true, firstEnergyStep must be provided
          if ~isfield(gridProperties, 'firstEnergyStep')
              error('firstEnergyStep must be provided when variableGrid is true.');
          end
          grid.firstEnergyStep = gridProperties.firstEnergyStep;

          % Perform self-diagnostic test
          grid.selfDiagnosticTest();
          % Compute the progression ratio 'a'
          grid.computeProgressionRatio();
          % Compute the nodes for variable energy grid (geometric progression)
          grid.computeGeometricGrid();

      else

          % Compute the nodes for uniform energy grid
          % grid.step = gridProperties.maxEnergy/grid.cellNumber;
          % grid.node = (0:grid.cellNumber)*grid.step;
          % grid.cell = ((1:grid.cellNumber)-0.5)*grid.step;
          grid.computeUniformGrid();

      end

      % Set smart grid properties if provided
      if isfield(gridProperties, 'smartGrid')
        grid.isSmart = true;
        grid.minEedfDecay = gridProperties.smartGrid.minEedfDecay;
        grid.maxEedfDecay = gridProperties.smartGrid.maxEedfDecay;
        grid.updateFactor = gridProperties.smartGrid.updateFactor;
      end

    end
    
    function computeProgressionRatio(grid)
        % Define the function F(a) for solving progression ratio 'a'
        F = @(a) exp(log((grid.maxEnergy / grid.firstEnergyStep) * (a - 1) + 1) / (grid.cellNumber)) - a;

        % Use fzero to solve for progression ratio 'a'
        options = optimset('Display', 'iter'); % show iterations
        solution = fzero(F, [1+1E-15, 3], options); %default: 1+1E-15

        grid.ratio = solution;

        % Check if the solution is valid
        if grid.ratio <= 1
            error('Invalid progression ratio "a". Ensure that maxEnergy and firstEnergyStep are consistent.');
        end
    end

    function selfDiagnosticTest(grid)
        % Self-diagnostic test to ensure the progression ratio 'a' is valid
        if grid.variableGrid
            % Check if the first energy step allows for increasing energy steps in geometric progression
            if grid.firstEnergyStep * grid.cellNumber >= grid.maxEnergy
                error(['Invalid progression ratio "a" (>= 1). Ensure that maxEnergy and firstEnergyStep are consistent. ' ...
                    'Consider using a smaller value for firstEnergyStep.']); 
            end
        end
    end

    function computeGeometricGrid(grid)
        % Compute node, cell and energy-step values for a variable energy grid (geometric progression)
        grid.node = zeros(1, grid.cellNumber + 1); % Initialize node array
        grid.energyStep = zeros(1, grid.cellNumber); % Initialize step sizes

        % Define the value of the first energy-step (u1)
        u1 = grid.firstEnergyStep;

        % Compute nodes and energy steps using geometric progression
        grid.node(1) = 0; % Starting energy
        for i = 1:grid.cellNumber
            grid.energyStep(i) = u1 * grid.ratio^(i - 1);
            grid.node(i + 1) = grid.node(i) + grid.energyStep(i);
        end

        % Compute cell values at midpoints between nodes
        grid.cell = 0.5 * (grid.node(1:end-1) + grid.node(2:end));

        % Set step property as the value of the first energy-step 
        %  for backward compatibility with uniform energy grid
        grid.step = grid.energyStep(1);
    end

    function computeUniformGrid(grid)
        % Compute node, cell and energy-step values for an uniform energy grid
        grid.node = linspace(0, grid.maxEnergy, grid.cellNumber + 1); % Uniformly spaced nodes
        grid.energyStep = diff(grid.node); % Constant step sizes

        % Compute cell values at midpoints between nodes
        grid.cell = 0.5 * (grid.node(1:end-1) + grid.node(2:end));

        % Set step property as the value of the first energy-step 
        %  for backward compatibility with uniform energy grid
        grid.step = grid.energyStep(1);

    end

    function updateMaxValue(grid, maxValue)
      % resize the grid with a new maximum value
      % grid.step = maxValue/grid.cellNumber;
      % grid.node = (0:grid.cellNumber)*grid.step;
      % grid.cell = ((1:grid.cellNumber)-0.5)*grid.step; 
      grid.maxEnergy = maxValue;
      if grid.variableGrid
          grid.computeProgressionRatio();
          grid.computeGeometricGrid();
      else
          grid.computeUniformGrid();
      end
      
      % broadcast change in the grid
      notify(grid, 'updatedMaxEnergy1'); % update interpolations of collision cross-sections
      notify(grid, 'updatedMaxEnergy2'); % update operators of the Boltzmann equation
    end

    function cellNumber = findCellNumber(grid, energy)
        % Calculate the cell number for a given energy value

        % Ensure that energy value is within grid bounds

        % TO BE CHECKED FOR VARIABLE ENERGY GRID
        if energy < grid.node(1)
            error('Energy value is below the lower grid bound.');
        elseif energy > grid.node(end)
            % This "trick" allows defining cellNumbers ABOVE the maximum
            % number of grid-cells, to automatically set to zero the rate coefficients 
            % of "large-threshold" cross sections
            cellNumber = floor(energy / (grid.maxEnergy / grid.cellNumber));
        else
            if grid.variableGrid
                % For variable energy grid, use the geometric progression formula
                u1 = grid.energyStep(1);
                a = grid.ratio;
                cellNumber = floor(log(energy / u1 * (a - 1) + 1) / log(a) + 1);
        
            else
                % For uniform energy grid, use linear interpolation
                cellNumber = floor(energy / (grid.maxEnergy / grid.cellNumber));

            end
        end
    end
  end

end
