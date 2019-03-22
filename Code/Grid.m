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

    node = [];              % values of the grid at node positions
    cell = [];              % values of the grid at cell positions (i.e. between two consecutive nodes)
    step = [];              % difference between the values at two consecutive nodes
    cellNumber = [];        % number of cells in the energy grid
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
      grid.step = gridProperties.maxEnergy/grid.cellNumber;
      grid.node = (0:grid.cellNumber)*grid.step;
      grid.cell = ((1:grid.cellNumber)-0.5)*grid.step; 
      if isfield(gridProperties, 'smartGrid')
        grid.isSmart = true;
        grid.minEedfDecay = gridProperties.smartGrid.minEedfDecay;
        grid.maxEedfDecay = gridProperties.smartGrid.maxEedfDecay;
        grid.updateFactor = gridProperties.smartGrid.updateFactor;
      end
    end
    
    function updateMaxValue(grid, maxValue)
      % resize the grid with a new maximum value
      grid.step = maxValue/grid.cellNumber;
      grid.node = (0:grid.cellNumber)*grid.step;
      grid.cell = ((1:grid.cellNumber)-0.5)*grid.step; 
      
      % broadcast change in the grid
      notify(grid, 'updatedMaxEnergy1');
      notify(grid, 'updatedMaxEnergy2');
    end

    function adjustedArray = adjustToGrid(grid, originalArray, threshold, mode)
      if strcmp(mode, 'nodes')
        x = grid.node(:);
      elseif strcmp(mode, 'cells')
        x = grid.cell(:);
      else
        error(' Unrecognised mode ''%s'' while adjusting data to Grid', mode)
      end
      adjustedArray = zeros(1, length(x(:)));
      if threshold < x(end)
        i = 2;
        for j = (ceil(threshold/grid.step)+1):length(x(:))
          x1 = originalArray(1,i-1);
          x2 = originalArray(1,i);
          while (x2 < x(j) && i < length(originalArray(1,:)))
            i = i+1;
            x1 = x2;
            x2 = originalArray(1,i);
          end
          if x2 >= x(j)
            y1 = originalArray(2,i-1);
            y2 = originalArray(2,i);
            adjustedArray(j) = y1 + (x(j)-x1)*(y2-y1)/(x2-x1);
            if adjustedArray(j) < 0
              adjustedArray(j) = 0.0;
            end
          else
            adjustedArray(j) = 0.0;
          end
        end
      end
    end

    function result = integrate(grid, integrand)
      % integrate evaluates the integral of the array 'integrand' with
      % respect to the Grid 'grid' using the trapezoidal rule.
      
      result = 0;
      if size(integrand)==size(grid.cell)
        for i = 1:length(grid.cell)
          result = result + integrand(i)*grid.step;
        end
      elseif size(integrand)==size(grid.node)
        for i = 1:length(grid.node)-1
          result = result + grid.step*0.5*(integrand(i)+integrand(i+1));
        end
      else
        error(' Integrand dimensions do not match with grid dimensions');
      end
    end
    
    function result = maxwellianRateCoeff(grid, crossSection, temperatureInEV)
      % maxwellianRateCoeff evaluates the rate coefficient of a certain
      % collision with cross section 'crossSection' whith a maxwellian eedf
      % characterised by temperature 'temperatureInEV'. The 'crossSection'
      % array must be adjusted to the 'grid', and the 'temperatureInEV' must
      % be in electron volts.
      %         result = sqrt(2e/m_e)*integral(sigma(u)*u*f(u)*du)
      
      
      if size(crossSection)==size(grid.cell)
        f = exp(-grid.cell/temperatureInEV);
        f = f/grid.integrate(f.*sqrt(grid.cell));
        result = grid.integrate(crossSection.*grid.cell.*f);
      elseif size(crossSection)==size(grid.node)
        f = exp(-grid.node/temperatureInEV);
        f = f/grid.integrate(f.*sqrt(grid.node));
        result = grid.integrate(crossSection.*grid.node.*f);
      else
        error(' Integrand dimensions do not match with grid dimensions');
      end
      result = result*sqrt(2.0*Constant.electronCharge/Constant.electronMass);
    end
    
  end

end
