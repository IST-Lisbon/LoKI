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

classdef Constant
  
  properties (Constant = true)
    
    % List of fundamental constants and units as obtained from the NIST database
    % (https://physics.nist.gov/cuu/Constants/index.html) on Jule 23rd 2018.
    
    boltzmann = 1.38064852e-23;           % Boltzmann constant in J/K
    electronCharge = 1.6021766208e-19;    % Electron charge in C
    electronMass = 9.10938356e-31;        % Electron mass in Kg
    unifiedAtomicMass = 1.660539040e-27;  % Unified Atomic Mass unit (UAM) in kg
    bohrRadius = 5.2917721067e-11;        % Bohr radius in m
    vacuumPermittivity = 8.854187817e-12; % Vacuum permittivity in F/m
    plank = 6.626070040e-34;              % Plank constant in J s
    speedOfLight = 299792458;             % Speed of light in vacuum in m/s
    atmosphereInPa = 101325;              % Standard atmosphere in Pa
    atmosphereInTorr = 760;               % Standard atmosphere in Torr (not obtained from NIST database)
    
  end
  
  methods (Static)
    
    function kBeV = boltzmannInEV()
      % boltzmannInEV Boltzmann constant in eV/K
      
      kBeV = Constant.boltzmann/Constant.electronCharge;
      
    end
    
    function hBar = plankReduced()
      % plankReduced Reduced Plank constant in J s
      
      hBar = Constant.plank/(2*pi);
      
    end
    
    function hBar = plankReducedInEV()
      % plankReducedInEV Reduced Plank constant in eV s
      
      hBar = Constant.plank/(2*pi*Constant.electronCharge);
      
    end

  end
end
