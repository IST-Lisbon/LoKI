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

classdef WorkingConditions < handle
  %WorkingConditions Class that contains the information about the working conditions of the simulation
  %   This class contains the information about the working conditins of the simulation. It is also the class
  %   responsible for handling changes in the working conditions by broadcasting any change in them to the
  %   corresponding objects (e.g. when the reduced field is updated the boltzmann object is notified, so it can
  %   reevaluate the field operator).
  
  properties
    gasPressure = [];
    gasTemperature = [];
    gasDensity = [];
    electronDensity = [];
    electronTemperature = [];
    chamberLength = [];
    chamberRadius = [];
    reducedField = [];
    reducedFieldSI = [];
    excitationFrequency = [];
    reducedExcFreqSI= [];
  end
  
  events
    updatedGasPressure
    updatedGasTemperature1
    updatedGasTemperature2
    updatedGasDensity
    updatedElectronDensity
    updatedElectronTemperature
    updatedChamberLength
    updatedReducedField
    updatedExcitationFrequency
  end
  
  methods (Access = public)
    
    function workCond = WorkingConditions(setup)
      
      for field = fieldnames(setup.info.workingConditions)'
        workCond.(field{1}) = setup.info.workingConditions.(field{1})(1);
      end
      workCond.gasDensity = workCond.gasPressure/(Constant.boltzmann*workCond.gasTemperature);
      workCond.reducedFieldSI = workCond.reducedField*1e-21;
      workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
      
    end
    
    function update(workCond, propertiesToUpdate, newValues)
      
      
      if ischar(propertiesToUpdate)
        propertiesToUpdate = {propertiesToUpdate};
      end
      
      for idx = 1:length(propertiesToUpdate)
        workCond.(propertiesToUpdate{idx}) = newValues(idx);
        
        switch propertiesToUpdate{idx}
          case 'gasPressure'
            workCond.gasDensity = workCond.gasPressure/(Constant.boltzmann*workCond.gasTemperature);
            workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
            notify(workCond, 'updatedGasPressure');
            notify(workCond, 'updatedGasDensity');
            notify(workCond, 'updatedExcitationFrequency');
            
          case 'gasTemperature'
            workCond.gasDensity = workCond.gasPressure/(Constant.boltzmann*workCond.gasTemperature);
            workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
            notify(workCond, 'updatedGasTemperature1');
            notify(workCond, 'updatedGasTemperature2');
            notify(workCond, 'updatedGasDensity');
            notify(workCond, 'updatedExcitationFrequency');
            
          case 'electronDensity'
            notify(workCond, 'updatedElectronDensity');
            
          case 'electronTemperature'
            notify(workCond, 'updatedElectronTemperature');
            
          case 'chamberLength'
            notify(workCond, 'updatedChamberLength');
            
          case 'chamberRadius'
            notify(workCond, 'updatedChamberRadius');
            
          case 'reducedField'
            workCond.reducedFieldSI = workCond.reducedField*1e-21;
            notify(workCond, 'updatedReducedField');
            
          case 'excitationFrequency'
            workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
            notify(workCond, 'updatedExcitationFrequency');
        end
      end
      
    end
    
    function workCondStruct = struct(workCond)
    % struct returns an structure with the properties of the object workCond
    
      workCondStruct = struct('gasPressure', workCond.gasPressure, 'gasTemperature', workCond.gasTemperature, ...
        'gasDensity', workCond.gasDensity, 'electronDensity', workCond.electronDensity, ...
        'electronTemperature', workCond.electronTemperature, 'chamberLength', workCond.chamberLength, ...
        'chamberRadius', workCond.chamberRadius, 'reducedField', workCond.reducedField, ...
        'reducedFieldSI', workCond.reducedFieldSI, 'excitationFrequency', workCond.excitationFrequency, ...
        'reducedExcFreqSI', workCond.reducedExcFreqSI);
      
    end
    
  end
  
end

