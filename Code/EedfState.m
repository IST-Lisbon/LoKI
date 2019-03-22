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

classdef EedfState < State
  
  properties
    
    isTarget = false;                 % true if the state is the target of a collision, false otherwise
    collisionArray = Collision.empty; % handle array to the collisions of which the target is state
    collisionArrayExtra = Collision.empty;

  end
  
  events
    
  end
  
  methods
    
    function state = EedfState(gas, ionCharg, eleLevel, vibLevel, rotLevel)
      persistent lastID;
      if isempty(lastID)
        lastID = 0;
      end
      lastID = lastID + 1;
      state.ID = lastID;
      state.gas = gas;
      state.ionCharg = ionCharg;
      state.eleLevel = eleLevel;
      state.vibLevel = vibLevel;
      state.rotLevel = rotLevel;
      if isempty(ionCharg)
        if isempty(rotLevel)
          if isempty(vibLevel)
            state.type = 'ele';
          else
            state.type = 'vib';
          end
        else
          state.type = 'rot';
        end
      else
        state.type = 'ion';
      end
      state.parent = EedfState.empty;
      state.siblingArray = EedfState.empty;
      state.childArray = EedfState.empty;
      state.addFamily;
      gas.stateArray(end+1) = state;
      state.evaluateName;
    end
    
  end
  
end