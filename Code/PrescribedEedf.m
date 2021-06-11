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

classdef PrescribedEedf < handle
  %PrescribedEedf Class that evaluates the EEDF corresponding to a prescribed (user defined) electron temperature.
  %
  %   PrescribedEedf evaluates a generalized EEDF, corresponding to a particular electron temperature, as well as the
  %   corresponding swarm parameters and electron impact rate coefficients. The generalized expression, that has as 
  %   limiting cases the Maxwellian and Druyvesteyn distribution functions, is as follows:
  % 
  %   f(u) = (gamma(5/(2g))^3/2)/(gamma(3/(2g))^5/2)*(2/(3KbTe))^3/2*exp(-(2u*gamma(5/(2g)/(3KbTe*gamma(3/(2g))))^g)
  % 
  %   where gamma is the gamma function (https://en.wikipedia.org/wiki/Gamma_function), Kb the boltzmann constant, 
  %   Te the electron temperature and g is the parameter thar controls the shape of the distribution function (g=1
  %   for Maxwellian and g=2 for Druyvesteyn).

  properties

    gasArray = Gas.empty;                   % handle array to the gas mixture
    energyGrid = Grid.empty;                % handle to the energy grid in which the Boltzmann equation is to be solved
    workCond = WorkingConditions.empty;     % handle to the working conditions of the simulation
                                            
    CARgases = [];            % gases described by the continuous approximation for rotations (CAR)

    totalCrossSection = [];   % total momentum transfer cross section
    elasticCrossSection = []; % total elastic cross section
    
    shapeParameter = [];          % value of the parameter controling the shape of the eedf

    eedf = [];                    % eedf
    power = struct.empty;         % power balance
    swarmParam = struct.empty;    % swarm parameters obtained with the eedf
    rateCoeffAll = struct.empty;  % rate coefficients obtained with the eedf and collisions in gasArray
    rateCoeffExtra = struct.empty;% extra rate coefficients for collisions not taken into account to obtain the eedf

  end

  events
    obtainedNewEedf;
  end

  methods (Access = public)

    function prescribedEedf = PrescribedEedf(setup)

      % store the gas array
      prescribedEedf.gasArray = setup.electronKineticsGasArray;
      
      % store gases for which the CAR is activated (in case there is any)
      if isfield(setup.info.electronKinetics, 'CARgases')
        prescribedEedf.CARgases = setup.info.electronKinetics.CARgases;
      end
      
      % store the energy grid and add corresponding listener
      prescribedEedf.energyGrid = setup.energyGrid;
      addlistener(prescribedEedf.energyGrid, 'updatedMaxEnergy2', @prescribedEedf.evaluateTotalAndElasticCrossSections);

      % store working conditions
      prescribedEedf.workCond = setup.workCond;
      
      % store the value of the parameter that controls the shape of the eedf (1 Maxwellian, 2 Druyvesteyn)
      prescribedEedf.shapeParameter = setup.info.electronKinetics.shapeParameter;

      % evaluate total momentum transfer and total elastic cross sections
      prescribedEedf.evaluateTotalAndElasticCrossSections();

    end

    function solve(prescribedEedf)
            
      % when the smart grid is activated the minimum number of decades of decay in the eedf is ensured
      if prescribedEedf.energyGrid.isSmart
        g = prescribedEedf.shapeParameter;
        Te = prescribedEedf.workCond.electronTemperature;
        decades = 0.5*(prescribedEedf.energyGrid.minEedfDecay + prescribedEedf.energyGrid.maxEedfDecay);
        maxEnergy = (decades/log10(exp(1)))^(1/g)*1.5*Te*gamma(3/(2*g))/gamma(5/(2*g));
        prescribedEedf.energyGrid.updateMaxValue(maxEnergy);
      end
      
      % evaluate EEDF and store it in the properties of the prescribedEedf object
      prescribedEedf.eedf  = prescribedEedf.evaluateEEDF();
      
      % evaluate power balance
      prescribedEedf.evaluatePower();

      % evaluate rate coefficients
      prescribedEedf.evaluateRateCoeff();

      % evaluate transport parameters
      prescribedEedf.evaluateSwarmParameters();

      % bradcast obtention of a solution for the EEDF
      notify(prescribedEedf, 'obtainedNewEedf');

    end
    
    function updateDensityDependencies(prescribedEedf)
    % updateDensityDependencies is a function that evaluates all the species density dependencies of the
    % prescribedEedf object.
      
      prescribedEedf.evaluateTotalAndElasticCrossSections();
      
    end

  end

  methods (Access = private)
    
    function evaluateTotalAndElasticCrossSections(prescribedEedf,~,~)
      
      % reset values to zero
      prescribedEedf.totalCrossSection = zeros(size(prescribedEedf.energyGrid.node));
      prescribedEedf.elasticCrossSection = zeros(size(prescribedEedf.energyGrid.node));
      % loop over each gas in the mixture
      for gas = prescribedEedf.gasArray
        % avoid dummy gasses 
        if isempty(gas.collisionArray)
          continue;
        end
        % evaluation of the mass ratio
        massRatio = Constant.electronMass/gas.mass;
        % loop over each collision with the gas
        for collision = gas.collisionArray
          % avoid effective collisions
          if strcmp(collision.type, 'Effective')
            continue;
          end
          % add collision cross section to the total momentum transfer cross section (also superelastic)
          prescribedEedf.totalCrossSection = prescribedEedf.totalCrossSection + collision.target.density*collision.crossSection;
          if collision.isReverse
            prescribedEedf.totalCrossSection = prescribedEedf.totalCrossSection + ...
              collision.productArray.density*collision.superElasticCrossSection;
          end
          % add elastic collision cross section to the total elastic cross section (weighted by the mass ratio)
          if strcmp(collision.type, 'Elastic')
            prescribedEedf.elasticCrossSection = prescribedEedf.elasticCrossSection + ...
              massRatio*collision.target.density*collision.crossSection;
            continue;
          end
        end
      end
    
    end
    
    function eedf = evaluateEEDF(prescribedEedf)
    
      % evaluate the distribution without normalization
      g = prescribedEedf.shapeParameter;
      Te = prescribedEedf.workCond.electronTemperature;
      energy = prescribedEedf.energyGrid.cell;
      eedf = gamma(5/(2*g))^(1.5)*gamma(3/(2*g))^(-2.5)*(2/(3*Te))^(1.5)*...
        exp(-(energy*gamma(5/(2*g))*2/(gamma(3/(2*g))*3*Te)).^g);

      % renormalising the distribution
      eedf = eedf/(sum(sqrt(prescribedEedf.energyGrid.cell).*eedf)*prescribedEedf.energyGrid.step);
      
    end

    function power = evaluatePower(prescribedEedf)
      
      % initialize power structure
      power = struct('field', 0, 'elasticNet', 0, 'elasticGain', 0, 'elasticLoss', 0, 'carNet', 0, 'carGain', 0, ...
        'carLoss', 0, 'excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
        'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
        'ionizationIne', 0, 'attachmentIne', 0, 'inelastic', 0, 'superelastic', 0, 'eDensGrowth', 0, ...
        'electronElectron', 0, 'gases', '');
      
      % save a local copy of the EEDF because of performance reasons
      eedfLocal = prescribedEedf.eedf;
      
      % save a local copy of the energy grid information
      energyGridLocal = prescribedEedf.energyGrid;    % energyGrid object
      N = energyGridLocal.cellNumber;                 % number of cells in the energy grid
      energyStep = energyGridLocal.step;              % energy step
      energyCell = energyGridLocal.cell;              % value of the energy at the cells
      energyNode = energyGridLocal.node;              % value of the energy at the nodes
      
      % multiplicative constant to obtain the right units
      factor = sqrt(2*Constant.electronCharge/Constant.electronMass);
      
      % evaluate power absorved per electron at unit gas density due to the electric field (auxiliary value)
      % evaluation of the elements of the electric field operator (Boltzmann)
      g_E = energyNode./(3*prescribedEedf.totalCrossSection);
      auxPowerField = factor*sum(eedfLocal.*(g_E(2:end)-g_E(1:end-1)));
      power.field = [];
      
      % auxiliary quantities needed to evaluate the elastic and CAR powers
      kTg = Constant.boltzmannInEV*prescribedEedf.workCond.gasTemperature; 
      aux1 = kTg+energyStep*0.5;
      aux2 = kTg-energyStep*0.5;
      
      % evaluate power absorved per electron at unit gas density due to elastic collisions
      % evaluation of the elements of the elastic operator (Boltzmann)
      g_c = 2*energyNode.^2.*prescribedEedf.elasticCrossSection;
      g_c(1) = 0;
      g_c(end) = 0;
      power.elasticNet = factor*sum(eedfLocal.*(g_c(2:end)*aux2-g_c(1:end-1)*aux1));
      power.elasticGain = factor*kTg*sum(eedfLocal.*(g_c(2:end)-g_c(1:end-1)));
      power.elasticLoss = power.elasticNet-power.elasticGain;
      
      % evaluate power absorved per electron at unit gas density due to rotations CAR
      if ~isempty(prescribedEedf.CARgases)
        % evaluation of the elements of the CAR operator (Boltzmann)
        sigma0B = 0;
        for gasName = prescribedEedf.CARgases
          gasID = Gas.find(gasName, prescribedEedf.gasArray);
          gas = prescribedEedf.gasArray(gasID);
          sigma0B = sigma0B + gas.fraction*gas.electricQuadrupoleMoment*gas.rotationalConstant;
        end
        g_car =4*energyNode.*(8.0*pi*sigma0B/(15.0*Constant.electronCharge));
        power.carNet = factor*sum(eedfLocal.*(g_car(2:end)*aux2-g_car(1:end-1)*aux1));
        power.carGain = factor*kTg*sum(eedfLocal.*(g_car(2:end)-g_car(1:end-1)));
        power.carLoss = power.carNet-power.carGain;
      end
      
      % evaluate power absorved per electron at unit gas density due to the inelastic/super-elastic collisions
      % loop over each gas in the mixture
      for gas = prescribedEedf.gasArray
        gasName = gas.name;
        % initialize power balance information of this gas
        gasPower = struct('excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
          'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
          'ionizationIne', 0, 'attachmentIne', 0);
        % loop over each collision with the gas
        for collision = gas.collisionArray
          % collision type
          collType = collision.type;
          % avoid Effective or Elastic collisions and collisions whos threshold is larger than the maximum energy
          if strcmp(collType, 'Effective') || strcmp(collType, 'Elastic') || collision.threshold > energyNode(end)
            continue;
          end
          % switch to lower case because of aesthetical reasons
          collType = lower(collType);
          % evaluate cross section at cell positions
          cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
          % evaluate departure cell
          lmin=floor(collision.threshold/energyStep);
          % add contribution to the power due to the inelastic collisions
          gasPower.([lower(collType) 'Ine']) = gasPower.([lower(collType) 'Ine']) - factor*collision.target.density*...
            energyStep*energyNode(lmin+1)*sum(eedfLocal(1+lmin:N).*energyCell(1+lmin:N).*cellCrossSection(1+lmin:N));
          % add contribution to the power due to the superelastic collisions
          if collision.isReverse
            statWeightRatio = collision.target.statisticalWeight/collision.productArray.statisticalWeight;
            gasPower.([lower(collType) 'Sup']) = gasPower.([lower(collType) 'Sup']) + factor*statWeightRatio*...
              collision.productArray.density*energyStep*energyNode(lmin+1)*sum(eedfLocal(1:N-lmin).*...
              energyCell(1+lmin:N).*cellCrossSection(1+lmin:N));
          end
        end
        % evaluate net values (for each gas)
        gasPower.excitationNet = gasPower.excitationIne + gasPower.excitationSup;
        gasPower.vibrationalNet = gasPower.vibrationalIne + gasPower.vibrationalSup;
        gasPower.rotationalNet = gasPower.rotationalIne + gasPower.rotationalSup;
        gasPower.inelastic = gasPower.excitationIne + gasPower.vibrationalIne + gasPower.rotationalIne + ...
          gasPower.ionizationIne + gasPower.attachmentIne;
        gasPower.superelastic = gasPower.excitationSup + gasPower.vibrationalSup + gasPower.rotationalSup;
        
        % evaluate net values (for the gas mixture)
        power.excitationIne = power.excitationIne + gasPower.excitationIne;
        power.excitationSup = power.excitationSup + gasPower.excitationSup;
        power.vibrationalIne = power.vibrationalIne + gasPower.vibrationalIne;
        power.vibrationalSup = power.vibrationalSup + gasPower.vibrationalSup;
        power.rotationalIne = power.rotationalIne + gasPower.rotationalIne;
        power.rotationalSup = power.rotationalSup + gasPower.rotationalSup;
        power.ionizationIne = power.ionizationIne + gasPower.ionizationIne;
        power.attachmentIne = power.attachmentIne + gasPower.attachmentIne;
        
        % store power balance information of this gas in the main structure
        power.gases.(gasName) = gasPower;
      end
      power.excitationNet = power.excitationIne + power.excitationSup;
      power.vibrationalNet = power.vibrationalIne + power.vibrationalSup;
      power.rotationalNet = power.rotationalIne + power.rotationalSup;
      power.inelastic = power.excitationIne + power.vibrationalIne + power.rotationalIne + power.ionizationIne + ...
        power.attachmentIne;
      power.superelastic = power.excitationSup + power.vibrationalSup + power.rotationalSup;
      
      % evaluate the electric field by ensuring perfect power balance
      power.field = -(power.elasticNet + power.carNet + power.inelastic + power.superelastic);
      power.balance = 0;
      power.relativeBalance = 0;
      prescribedEedf.workCond.update('reducedField', sqrt(power.field/auxPowerField)/1e-21);
      
      % evaluate reference power (chanel of max power)
      powerTypes = {'field' 'elasticGain' 'carGain' 'excitationSup' 'vibrationalSup' 'rotationalSup'};
      powerValues = [power.field power.elasticGain power.carGain power.excitationSup power.vibrationalSup ...
        power.rotationalSup];
      [maxPower, referencePowerIdx] = max(powerValues);
      power.reference = maxPower;
      power.referenceType = powerTypes{referencePowerIdx};
      
      % store power balance information in the prescribedEedf properties
      prescribedEedf.power = power;

    end

    function swarmParam = evaluateSwarmParameters(prescribedEedf)
      
      % save local copy of total momentum transfer cross section
      totalCrossSectionAux = prescribedEedf.totalCrossSection;
      
      % initialize transport parameters structure
      swarmParam = struct('redDiffCoeff', [], 'redMobility', [], 'redDiffCoeffEnergy', [], 'redMobilityEnergy', [], ...
        'redTownsendCoeff', [], 'redAttCoeff', [], 'meanEnergy', [], 'characEnergy', [], 'Te', [], 'driftVelocity', []);

      % evaluate reduced diffusion coefficient
      swarmParam.redDiffCoeff = (2.0/3.0)*sqrt(2.0*Constant.electronCharge/Constant.electronMass)*...
        sum(prescribedEedf.energyGrid.cell.*prescribedEedf.eedf./(prescribedEedf.totalCrossSection(1:end-1)+...
        totalCrossSectionAux(2:end)))*prescribedEedf.energyGrid.step;

      % evaluate reduced mobility
      swarmParam.redMobility = -sqrt(2.0*Constant.electronCharge/Constant.electronMass)/3.0*...
        sum(prescribedEedf.energyGrid.node(2:end-1).*(prescribedEedf.eedf(2:end)-prescribedEedf.eedf(1:end-1))./...
        totalCrossSectionAux(2:end-1));
      
      % evaluate reduced energy diffusion coefficient
      swarmParam.redDiffCoeffEnergy = (2.0/3.0)*sqrt(2.0*Constant.electronCharge/Constant.electronMass)*...
        sum(prescribedEedf.energyGrid.cell.^2.*prescribedEedf.eedf./(totalCrossSectionAux(1:end-1)+...
        totalCrossSectionAux(2:end)))*prescribedEedf.energyGrid.step;
      
      % evaluate reduced energy mobility
      swarmParam.redMobilityEnergy = -sqrt(2.0*Constant.electronCharge/Constant.electronMass)/3.0*...
        sum(prescribedEedf.energyGrid.node(2:end-1).^2.*(prescribedEedf.eedf(2:end)-prescribedEedf.eedf(1:end-1))./...
        totalCrossSectionAux(2:end-1));
      
      % evaluate drift velocity
      swarmParam.driftVelocity = swarmParam.redMobility*prescribedEedf.workCond.reducedFieldSI;
      
      % evaluate reduced Townsend coefficient
      totalIonRateCoeff = 0;
      for gas = prescribedEedf.gasArray
        for collision = gas.collisionArray
          if strcmp(collision.type, 'Ionization')
            totalIonRateCoeff = totalIonRateCoeff + collision.target.density*collision.ineRateCoeff;
          end
        end
      end
      swarmParam.redTownsendCoeff = totalIonRateCoeff / swarmParam.driftVelocity;

      % evaluate reduced attachment coefficient
      totalAttRateCoeff = 0;
      for gas = prescribedEedf.gasArray
        for collision = gas.collisionArray
          if strcmp(collision.type, 'Attachment')
            totalAttRateCoeff = totalAttRateCoeff + collision.target.density*collision.ineRateCoeff;
          end
        end
      end
      swarmParam.redAttCoeff = totalAttRateCoeff / swarmParam.driftVelocity;

      % evaluate mean energy
      swarmParam.meanEnergy = sum(prescribedEedf.energyGrid.cell(1:end).^(1.5).*prescribedEedf.eedf(1:end))*...
        prescribedEedf.energyGrid.step;

      % evaluate characteristic energy
      swarmParam.characEnergy = swarmParam.redDiffCoeff/swarmParam.redMobility;

      % evaluate electron temperature
      swarmParam.Te = prescribedEedf.workCond.electronTemperature; %swarmParam.meanEnergy*2/3;

      % store swarm parameters information in the prescribedEedf properties
      prescribedEedf.swarmParam = swarmParam;

    end

    function [rateCoeffAll, rateCoeffExtra] = evaluateRateCoeff(prescribedEedf)

      % initialize rateCoeffAll structure
      rateCoeffAll = struct.empty;
      rateCoeffExtra = struct.empty;

      % evaluate rate coefficient for all collision
      for gas = prescribedEedf.gasArray
        % collisions taken into account for solving the eedf
        for collision = gas.collisionArray
          rateCoeffAll(end+1).collID = collision.ID;
          [ineRate, supRate] = collision.evaluateRateCoeff(prescribedEedf.eedf);
          rateCoeffAll(end).value = [ineRate, supRate];
          rateCoeffAll(end).collDescription = collision.description;
        end
        % collisions not taken into account for solving the eedf
        for collision = gas.collisionArrayExtra
          rateCoeffExtra(end+1).collID = collision.ID;
          [ineRate, supRate] = collision.evaluateRateCoeff(prescribedEedf.eedf);
          rateCoeffExtra(end).value = [ineRate, supRate];
          rateCoeffExtra(end).collDescription = collision.description;
        end
      end
      
      % store rate coefficients information in the prescribedEedf properties
      prescribedEedf.rateCoeffAll = rateCoeffAll;
      prescribedEedf.rateCoeffExtra = rateCoeffExtra;

    end

  end

end
