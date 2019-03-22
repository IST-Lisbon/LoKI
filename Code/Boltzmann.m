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

classdef Boltzmann < handle
  %Boltzmann Class that solves the boltzmann equation under certain conditions to
  %obtain an EEDF
  %
  %   Boltzmann creates the different boltzmann matrices corresponding
  %   to the different elements of the boltzmann equation. In order to
  %   create a BoltzmannEq object, a gasArray with the gas mixture
  %   information and an energyGrid must be sent to the constructor.
  
  properties
    
    gasArray = Gas.empty;                   % handle array to the gas mixture
    energyGrid = Grid.empty;                % handle to the energy grid in which the Boltzmann equation is to be solved
    workCond = WorkingConditions.empty;     % handle to the working conditions of the simulation
    CARgases = [];                          % gases described by the continuous approximation for rotations (CAR)
    
    ionCollOpType;            % describes the type of ionization collision operator
    eDensGrowthModel;         % describes the electron density growth model
    includeEECollisions;      % indicates if electron-electrons collisions are accounted for
    
    nonLinearAlgorithm;       % algorithm used to solve the non-linear operators (ionization, attachment or e-e)
    mixingParameter;          % solution mixing fraction used on the iterative scheme of the growth routines
    maxEedfRelError;          % maximum relative difference for the eedf between two consecutive iterations (non-linear routines)
    maxPowerBalanceRelError;  % maximum value for the relative power balance (threshold for the warning message)
    odeOptions;               % user defined options for the ODE solver (used when nonLinearAlgorithm is iterativeSolution)

    totalCrossSection = [];   % total momentum transfer cross section
    elasticCrossSection = []; % total elastic cross section
    
    g_E = [];                 % elements of the field operator
    g_c = [];                 % elements of the elastic operator
    g_car = [];               % elements of the CAR operator
    g_fieldSpatialGrowth = [];% elements of the field operator with spatial electron density growth model
    g_fieldTemporalGrowth = [];% elements of the field operator with temporal electron density growth model
    Afinal = [];              % elements of the energy upflux from electron-electron collisions (used in power balance)
    Bfinal = [];              % elements of the energy downflux from electron-electron collisions (used in power balance)
    
    ionizationThresholdIsSmallerThanUmax;% indicates if at least one ionization threshold is within the energy grid range
    includeNonConservativeIonization = false; % boolean that indicates if growth models are to be used because of ionization
    includeNonConservativeAttachment = false; % boolean that indicates if growth models are to be used because of attachment
    attachmentThresholdIsSmallerThanUmax;% indicates if at least one attachment threshold is within the energy grid range
    
    CIEff;                    % effective ionization rate coefficient
    alphaRedEff;              % reduced first Townsend ionization coefficient
    alphaee;                  % alpha constant of electron-electron collisions
    
    fieldMatrix = [];         % matrix of the electric field operator of the boltzmann equation (continuous term)
    fieldMatrixTempGrowth = [];% matrix of the electric field operator of the boltzmann equation with temporal growth model
    fieldMatrixSpatGrowth = [];% matrix of another part of the electric field operator of the boltzmann equation with spatial growth model
    elasticMatrix = [];       % matrix of the elastic collision operator of the boltzmann equation (continuous term)
    CARMatrix = [];           % matrix of the CAR operator of the boltzmann equation (continuous term)
    continuousMatrix = [];    % matrix of the continuous operator of the boltzmann equation (sum of all continuos terms)
    inelasticMatrix = [];     % matrix of the discrete operator of the boltzmann equation (sum of all discrete terms)
    ionizationConservativeMatrix = []; % matrix of the conservative ionization operator of the boltzmann equation (discrete terms)
    ionizationMatrix = [];    % matrix of the non-conservative ionization operator of the boltzmann equation (discrete terms)
    attachmentConservativeMatrix = []; % matrix of the conservative attachment operator of the boltzmann equation (discrete terms)
    attachmentMatrix = [];    % matrix of the non-conservative attachment operator of the boltzmann equation (discrete terms)
    ionTemporalGrowth= [];    % matrix of the electron-density temporal growth term
    ionSpatialGrowthD = [];   % matrix of the electron-density spatial growth diffusion term
    ionSpatialGrowthU = [];   % matrix of the electron-density spatial growth mobility term
    Aee = [];                 % auxiliary matrix to calculate upflux vector of electron-electron collisions
    Bee = [];                 % auxiliary matrix to calculate downflux vector of electron-electron collisions
    
    eedf = [];                    % eedf obtained as solution to the boltzmann equation
    power = struct.empty;         % power balance of the boltzmann equation
    swarmParam = struct.empty;    % swarm parameters obtained with the eedf
    rateCoeffAll = struct.empty;  % rate coefficients obtained with the eedf and collisions in gasArray
    rateCoeffExtra = struct.empty;% extra rate coefficients for collisions not taken into account to obtain the eedf
    firstAnisotropy = [];         % value of the first anisotropy of the EDF (two term approximation)
    
  end
  
  events
    obtainedNewEedf;
  end
  
  methods (Access = public)
    
    function boltzmann = Boltzmann(setup)
      
      % store the gas array
      boltzmann.gasArray = setup.electronKineticsGasArray;
      
      % store the energy grid and add corresponding listener
      boltzmann.energyGrid = setup.energyGrid;
      addlistener(boltzmann.energyGrid, 'updatedMaxEnergy2', @boltzmann.evaluateMatrix);
      
      % store working conditions and add corresponding listeners
      boltzmann.workCond = setup.workCond;
      %       addlistener(workCond, 'updatedGasTemperature', @boltzmann.);
      %       addlistener(workCond, 'updatedGasDensity', @boltzmann.);
      addlistener(boltzmann.workCond, 'updatedReducedField', @boltzmann.evaluateFieldOperator);
      addlistener(boltzmann.workCond, 'updatedExcitationFrequency', @boltzmann.evaluateFieldOperator);
      
      % store gases for which the CAR is activated (in case there is any)
      if isfield(setup.info.electronKinetics, 'CARgases')
        boltzmann.CARgases = setup.info.electronKinetics.CARgases;
      end
      
      % store electron-impact ionization operator type
      boltzmann.ionCollOpType = setup.info.electronKinetics.ionizationOperatorType;
      
      % store electron density growth model type
      boltzmann.eDensGrowthModel = setup.info.electronKinetics.growthModelType;
      
      % store electron-electron collisions setups
      boltzmann.includeEECollisions = setup.info.electronKinetics.includeEECollisions;
      
      % store information about the numerical details on how to solve the Boltzmann equation
      boltzmann.nonLinearAlgorithm = setup.info.electronKinetics.numerics.nonLinearRoutines.algorithm;
      boltzmann.maxEedfRelError = setup.info.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError;
      if strcmp(boltzmann.nonLinearAlgorithm, 'mixingDirectSolutions') 
        boltzmann.mixingParameter = setup.info.electronKinetics.numerics.nonLinearRoutines.mixingParameter;
      else
        % store configuration of the ODE solver
        options = odeset();
        for parameter = fields(options)'
          if isfield(setup.info.electronKinetics.numerics.nonLinearRoutines, 'odeSetParameters') && ...
              isfield(setup.info.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, parameter{1})
            options.(parameter{1}) = setup.info.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.(parameter{1});
          end
        end
        boltzmann.odeOptions = options;
      end
      boltzmann.maxPowerBalanceRelError = setup.info.electronKinetics.numerics.maxPowerBalanceRelError;
      
      % allocate memory for different properties
%       boltzmann.totalCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
%       boltzmann.elasticCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
%       boltzmann.fieldMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.fieldMatrixTempGrowth = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.elasticMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
      boltzmann.CARMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.g_car = zeros(1,boltzmann.energyGrid.cellNumber+1);
%       boltzmann.continuousMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.discreteMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.ionizationIneMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
      
      % initialize quantities of the ionization and electron-electron routines
%       boltzmann.ionizationMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.Aee = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.Bee = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.ionSpatialGrowthD  = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.ionSpatialGrowthU = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.fieldMatrixSpatGrowth = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.fieldMatrixTempGrowth = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.ionTemporalGrowth = zeros(boltzmann.energyGrid.cellNumber);
      
      %evaluate boltzmann matrix
      boltzmann.evaluateMatrix();
      
    end
    
    function solve(boltzmann)
      
      % save appropiate method for the selectec non-linear algorithm
      switch boltzmann.nonLinearAlgorithm
        case 'mixingDirectSolutions'
          nonLinearSolver = str2func('mixingDirectSolutions');
        case 'iterativeSolution'
          nonLinearSolver = str2func('iterativeSolution');
      end

      % check for the presence of non-linear operators in the Boltzmann equation
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment || ...
          boltzmann.includeEECollisions
        % iterate non-linear solver until convergence
        auxEedf = nonLinearSolver(boltzmann);
      else
        % invert the boltzmann matrix to obtain an eedf (stationary solution)
        auxEedf = boltzmann.invertLinearMatrix;
      end
      
      % when the smart grid is activated the minimum number of decades of decay in the eedf is ensured
      if boltzmann.energyGrid.isSmart
        decades = log10(auxEedf(1))-log10(auxEedf(end));
        while decades < boltzmann.energyGrid.minEedfDecay
          % increase maximum value of the energy grid
          boltzmann.energyGrid.updateMaxValue(boltzmann.energyGrid.node(end)*(1+boltzmann.energyGrid.updateFactor));
          % check for the presence of non-linear operators in the Boltzmann equation
          if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment || ...
              boltzmann.includeEECollisions
            % iterate non-linear solver until convergence
            auxEedf = nonLinearSolver(boltzmann);
          else
            % invert the boltzmann matrix to obtain an eedf (stationary solution)
            auxEedf = boltzmann.invertLinearMatrix;
          end
          % check decades of decay
          decades = log10(auxEedf(1))-log10(auxEedf(end));
        end
        while decades > boltzmann.energyGrid.maxEedfDecay
          % decrease maximum value of the energy grid
          boltzmann.energyGrid.updateMaxValue(boltzmann.energyGrid.node(end)/(1+boltzmann.energyGrid.updateFactor));
          if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment || ...
              boltzmann.includeEECollisions
            % iterate non-linear solver until convergence
            auxEedf = nonLinearSolver(boltzmann);
          else
            % invert the boltzmann matrix to obtain an eedf (stationary solution)
            auxEedf = boltzmann.invertLinearMatrix;
          end
          % check decades of decay
          decades = log10(auxEedf(1))-log10(auxEedf(end));
        end
      end
      
      % evaluate power balance
      boltzmann.evaluatePower(true);
      
      % evaluate rate coefficients
      boltzmann.evaluateRateCoeff();
      
      % evaluate transport parameters
      boltzmann.evaluateSwarmParameters();
      
      % evaluate first anisotropy 
      boltzmann.evaluateFirstAnisotropy();
      
      % bradcast obtention of a solution for the boltzmann equation
      notify(boltzmann, 'obtainedNewEedf');
      
    end
    
    function updateDensityDependencies(boltzmann)
    % updateDensityDependencies is a function that evaluates all the species density dependencies of the boltzmann
    % object.
      
      boltzmann.evaluateMatrix();
      
    end
    
  end
  
  methods (Access = private)
    
    function evaluateMatrix(boltzmann, ~, ~)
      
      % evaluate continuous operators of the Boltzmann equation (field+elastic+CAR operators)
      boltzmann.evaluateContinuousOperators();
      
      % evaluate discrete operators of the Boltzmann equation (inelastic and superelastic operators without including
      % the contributions from ionization nor attachment)
      boltzmann.evaluateInelasticOperators();
      
      % evaluate ionization operator of the Boltzmann equation
      boltzmann.evaluateIonizationOperator();
      
      % evaluate attachment operator of the Boltzmann equation
      boltzmann.evaluateAttachmentOperator();
      
    end
    
    function evaluateContinuousOperators(boltzmann)
      
      % evaluate total and elastic momentum transfer cross sections for continuous operators
      boltzmann.evaluateTotalAndElasticCrossSections();
      
      % evaluate field operator
      boltzmann.evaluateFieldOperator();
      
      % evaluate elastic operator
      boltzmann.evaluateElasticOperator();
      
      % evaluate CAR operator
      if ~isempty(boltzmann.CARgases)
        boltzmann.evaluateCAROperator();
      end
      
    end
    
    function evaluateTotalAndElasticCrossSections(boltzmann)
      
      % reset values to zero
      boltzmann.totalCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
      boltzmann.elasticCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
      % loop over each gas in the mixture
      for gas = boltzmann.gasArray
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
          boltzmann.totalCrossSection = boltzmann.totalCrossSection + collision.target.density*collision.crossSection;
          if collision.isReverse
            boltzmann.totalCrossSection = boltzmann.totalCrossSection + ...
              collision.productArray.density*collision.superElasticCrossSection;
          end
          % add elastic collision cross section to the total elastic cross section (weighted by the mass ratio)
          if strcmp(collision.type, 'Elastic')
            boltzmann.elasticCrossSection = boltzmann.elasticCrossSection + ...
              massRatio*collision.target.density*collision.crossSection;
            continue;
          end
        end
      end
      
    end
    
    function evaluateFieldOperator(boltzmann, ~, ~)
      % evaluateFieldOperator is in charge of the evaluation of the elements of the field operator. The elements of
      % the operator are given by the following expression:
      %
      %              g_E(u) = -N*sqrt(2*e/me)*EoN^2/(sigmaT(u)*(1+WoN^2/(2*e*u*sigmaT^2(u)/me)))
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % definition of intermediate variables
      e = Constant.electronCharge;                % electron charge
      me = Constant.electronMass;                 % electron mass
      EoN = boltzmann.workCond.reducedFieldSI;    % reduced electric field (SI units)
      WoN = boltzmann.workCond.reducedExcFreqSI;  % reduced angular exitation frequency (SI units)
      
      % evaluation of the elements of the operator
      boltzmann.g_E = EoN^2*boltzmann.energyGrid.node./(3*boltzmann.totalCrossSection.*(1+me*WoN^2./...
        (2*e*boltzmann.energyGrid.node.*boltzmann.totalCrossSection.^2)));
      
      % evaluation of the boundary conditions
      boltzmann.g_E(1) = 0;
      boltzmann.g_E(end) = 0;
      
      % loading of the matrix
      if isempty(boltzmann.fieldMatrix)
        boltzmann.fieldMatrix = zeros(boltzmann.energyGrid.cellNumber);
      end
      for k = 1:boltzmann.energyGrid.cellNumber
        boltzmann.fieldMatrix(k,k) = -(boltzmann.g_E(k)+boltzmann.g_E(k+1))/boltzmann.energyGrid.step^2;
        if k>1
          boltzmann.fieldMatrix(k,k-1) = boltzmann.g_E(k)/boltzmann.energyGrid.step^2;
        end
        if k<boltzmann.energyGrid.cellNumber
          boltzmann.fieldMatrix(k,k+1) = boltzmann.g_E(k+1)/boltzmann.energyGrid.step^2;
        end
      end
      
    end
    
    function evaluateElasticOperator(boltzmann)
      % evaluateElasticOperator is in charge of the evaluation of the elements of the elastic collision operator. The
      % elements of the operator are given by the following expression:
      %
      %                               g_c(u) = -N*sqrt(2*e/me)*2*u^2*sigmaC(u)
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % definition of intermediate variables
      kb = Constant.boltzmannInEV;                % boltzmann constant in eV
      Tg = boltzmann.workCond.gasTemperature;     % gas temperature in K
      factor1 = (kb*Tg/boltzmann.energyGrid.step+0.5)/boltzmann.energyGrid.step;
      factor2 = (kb*Tg/boltzmann.energyGrid.step-0.5)/boltzmann.energyGrid.step;
      
      % evaluation of the elements of the operator
      boltzmann.g_c = 2*boltzmann.energyGrid.node.^2.*boltzmann.elasticCrossSection;
      
      % evaluation of the boundary conditions
      boltzmann.g_c(1) = 0;
      boltzmann.g_c(end) = 0;
      
      % loading of the matrix
      if isempty(boltzmann.elasticMatrix)
        boltzmann.elasticMatrix = zeros(boltzmann.energyGrid.cellNumber);
      end
      for k = 1:boltzmann.energyGrid.cellNumber
        boltzmann.elasticMatrix(k,k) = -(boltzmann.g_c(k)*factor1+boltzmann.g_c(k+1)*factor2);
        if k>1
          boltzmann.elasticMatrix(k,k-1) = boltzmann.g_c(k)*factor2;
        end
        if k<boltzmann.energyGrid.cellNumber
          boltzmann.elasticMatrix(k,k+1) = boltzmann.g_c(k+1)*factor1;
        end
      end
      
    end
    
    function evaluateCAROperator(boltzmann)
      % evaluateCAROperator is in charge of the evaluation of the elements of the continuous approximation for rotation
      % (CAR) operator with Chapman-Cowling correction. The operator is developed for gases well described by the
      % rotational cross sections proposed by Gerjuoy and Stein (complete development of the operator can be found in
      % "M A Ridenti, L L Alves, V Guerra and J Amorim, Plasma Sources Sci. Technol. 24 (2015) 035002 (16pp)).
      %
      % The elements of the operator are given by the following expression:
      %
      %                               g_car(u) = -N*sqrt(2*e/me)*4*u*sigma0*B
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % definition of intermediate variables
      kb = Constant.boltzmannInEV;                % boltzmann constant in eV
      Tg = boltzmann.workCond.gasTemperature;     % gas temperature in K
      factor1 = (kb*Tg/boltzmann.energyGrid.step+0.5)/boltzmann.energyGrid.step;
      factor2 = (kb*Tg/boltzmann.energyGrid.step-0.5)/boltzmann.energyGrid.step;
      
      % evaluation of the elements of the operator
      sigma0B = 0;
      for gasName = boltzmann.CARgases
        gasID = Gas.find(gasName, boltzmann.gasArray);
        gas = boltzmann.gasArray(gasID);
        sigma0B = sigma0B + gas.fraction*gas.electricQuadrupoleMoment*gas.rotationalConstant;
      end
      sigma0B = 8.0*pi*sigma0B/(15.0*Constant.electronCharge);
      boltzmann.g_car =4*boltzmann.energyGrid.node.*sigma0B;
      
      % evaluation of the boundary conditions
      boltzmann.g_car(1) = 0;
      boltzmann.g_car(end) = 0;
      
      % loading of the matrix
      if isempty(boltzmann.CARMatrix)
        boltzmann.CARMatrix = zeros(boltzmann.energyGrid.cellNumber);
      end
      for k = 1:boltzmann.energyGrid.cellNumber
        boltzmann.CARMatrix(k,k) = -(boltzmann.g_car(k)*factor1+boltzmann.g_car(k+1)*factor2);
        if k>1
          boltzmann.CARMatrix(k,k-1) = boltzmann.g_car(k)*factor2;
        end
        if k<boltzmann.energyGrid.cellNumber
          boltzmann.CARMatrix(k,k+1) = boltzmann.g_car(k+1)*factor1;
        end
      end
      
    end
    
    function evaluateInelasticOperators(boltzmann)
      % evaluateInelasticOperators is in charge of the evaluation of the elements of the inelastic and superelastic 
      % collision operator (it does not include the contribution due to ionization nor attachment collisions)
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % define local copies of variables used multiple times along the function
      energyNode = boltzmann.energyGrid.node;
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      
      % allocate memory for the matrix
      inelasticMatrixAux = zeros(cellNumber);
      
      % loop over each gas in the mixture
      for gas = boltzmann.gasArray
        % loop over each collision with the gas
        for collision = gas.collisionArray
          collisionType = collision.type;
          threshold = collision.threshold;
          % avoid Effective, Elastic, Ionization and attachment collisions, also avoid collisions which threshold is 
          % larger than the maximum energy or smaller than the energy step
          if strcmp(collisionType, 'Effective') || strcmp(collisionType, 'Elastic') || ...
              strcmp(collisionType, 'Ionization') || strcmp(collisionType, 'Attachment') || ...
              threshold > energyNode(end) || threshold < energyStep
            continue;
          end
          % evaluate numerical threshold
          numThreshold=floor(threshold/energyStep);
          % evaluate cross section at cell positions
          nodeCrossSection = collision.crossSection;
          cellCrossSection = 0.5*(nodeCrossSection(1:end-1)+nodeCrossSection(2:end));
          % store density of the target of the collision
          targetDensity = collision.target.density;
          % load matrix
          for k=1:cellNumber
            if(k<=cellNumber-numThreshold)
              inelasticMatrixAux(k,k+numThreshold) = inelasticMatrixAux(k,k+numThreshold) + ...
                targetDensity*energyCell(k+numThreshold)*cellCrossSection(k+numThreshold);
            end
            inelasticMatrixAux(k,k) = inelasticMatrixAux(k,k) - targetDensity*energyCell(k)*cellCrossSection(k);
          end
          % include superelastic collision
          if collision.isReverse
            statWeightRatio = collision.target.statisticalWeight/collision.productArray.statisticalWeight;
            productDensity = collision.productArray.density;
            for k=1:cellNumber
              if(k>=numThreshold+1)
                inelasticMatrixAux(k,k-numThreshold) = inelasticMatrixAux(k,k-numThreshold) + statWeightRatio*...
                  productDensity*energyCell(k)*cellCrossSection(k);
              end
              if(k<=cellNumber-numThreshold)
                inelasticMatrixAux(k,k) = inelasticMatrixAux(k,k) - statWeightRatio*...
                  productDensity*energyCell(k+numThreshold)*cellCrossSection(k+numThreshold);
              end
            end
          end
        end
      end
      
      % save discrete matrix in the boltzmann object
      boltzmann.inelasticMatrix = inelasticMatrixAux;
      
    end
    
    function evaluateIonizationOperator(boltzmann)
      % evaluateIonizationOperator is in charge of the evaluation of the elements of the ionization
      % conservative or non-conservative collision operator
      
      % define local copies of variables used multiple times along the function
      energyNode = boltzmann.energyGrid.node;
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      ionConservativeMatrix = zeros(cellNumber);
      ionNonConservativeMatrix = zeros(cellNumber);
      thresholdIsSmallerThanUmax = false;
      
      for gas = boltzmann.gasArray
        % loop over each collision with the gas
        for collision = gas.collisionArray
          threshold = collision.threshold;
          if strcmp(collision.type, 'Ionization') && threshold <= energyNode(end)
            % set flag for non-conservative algorithms
            thresholdIsSmallerThanUmax = true;
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1) + collision.crossSection(2:end));
            % local copy of target density
            density = collision.target.density;
            % evaluate numerical threhold
            numThreshold=floor(threshold/energyStep);
            
            % evaluation of nonconservative ionization collisional operator (in case it is needed)
            switch boltzmann.ionCollOpType
              case 'oneTakesAll'
                for k=1:cellNumber
                  if (k<=cellNumber-numThreshold)
                    ionNonConservativeMatrix(k,k+numThreshold) = ionNonConservativeMatrix(k,k+numThreshold) + ...
                      density*energyCell(k+numThreshold)*cellCrossSection(k+numThreshold);
                  end
                  ionNonConservativeMatrix(k,k) =  ionNonConservativeMatrix(k,k) - ...
                    density*energyCell(k)*cellCrossSection(k);
                  ionNonConservativeMatrix(1,k) = ionNonConservativeMatrix(1,k) + ...
                    density*energyCell(k)*cellCrossSection(k);
                end
              case 'equalSharing'
                % using linear indeces in order to improve performance (just in this particular calculation)
                ionNonConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) = ...
                  ionNonConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) - ...
                  density*energyCell(1:cellNumber).*cellCrossSection(1:cellNumber);
                for k=1:(cellNumber-numThreshold)/2
                  ionNonConservativeMatrix(k,2*k+numThreshold) = ionNonConservativeMatrix(k,2*k+numThreshold) + ...
                    4*density*energyCell(2*k+numThreshold)*cellCrossSection(2*k+numThreshold);
                end
              case 'usingSDCS'
                W = gas.OPBParameter;
                if isempty(W)
                  W = threshold;
                end
                for k=1:cellNumber
                  half = floor((k-numThreshold)/2);
                  final = min([2*k+numThreshold cellNumber]);
                  ionNonConservativeMatrix(k,k)= ionNonConservativeMatrix(k,k) - ...
                    energyCell(k)*density*energyStep*cellCrossSection(k)*...
                    sum(1/(W*atan((energyCell(k)-threshold)/(2*W)))./(1+(energyCell(1:half)/W).^2));
                  if k+numThreshold+1<=cellNumber
                    ionNonConservativeMatrix(k,(k+numThreshold+1):final) = ...
                      ionNonConservativeMatrix(k,(k+numThreshold+1):final) + energyCell((k+numThreshold+1):final)*...
                      density*energyStep.*cellCrossSection((k+numThreshold+1):final)./...
                      (W*atan((energyCell((k+numThreshold+1):final)-threshold)/(2*W)))./...
                      (1+(energyCell(((k+numThreshold+1):final)-k-numThreshold)/W).^2);
                  end
                  ionNonConservativeMatrix(k,(2*k+numThreshold):cellNumber) = ...
                    ionNonConservativeMatrix(k,(2*k+numThreshold):cellNumber) + ...
                    energyCell((2*k+numThreshold):cellNumber)*density*energyStep.*...
                    cellCrossSection((2*k+numThreshold):cellNumber)./...
                    (W*atan((energyCell(2*k+numThreshold:cellNumber)-threshold)/(2*W)))./(1+(energyCell(k)/W).^2);
                end
            end
            % avoid writing of conservative collision operator if threshold is smaller than the energy step
            if numThreshold==0
              continue
            end
            % evaluation of conservative ionization collisional operator (always needed)
            for k=1:cellNumber
              if (k<=cellNumber-numThreshold)
                ionConservativeMatrix(k,k+numThreshold) = ionConservativeMatrix(k,k+numThreshold) + ...
                  density*energyCell(k+numThreshold)*cellCrossSection(k+numThreshold);
              end
              ionConservativeMatrix(k,k) = ionConservativeMatrix(k,k) - density*energyCell(k)*cellCrossSection(k);
            end
          end
        end
      end
      
      % save matrixes in as properties of the boltzmann object
      boltzmann.ionizationConservativeMatrix = ionConservativeMatrix;
      boltzmann.ionizationMatrix = ionNonConservativeMatrix;
      boltzmann.ionizationThresholdIsSmallerThanUmax = thresholdIsSmallerThanUmax;
      
      if ~strcmp(boltzmann.ionCollOpType, 'conservative') && thresholdIsSmallerThanUmax
        boltzmann.includeNonConservativeIonization = true;
      end
      
    end
    
    function evaluateAttachmentOperator(boltzmann)
      % evaluateAttachmentOperator is in charge of the evaluation of the elements of the attachment
      % conservative or non-conservative collision operator
      
      % define local copies of variables used multiple times along the function
      energyNode = boltzmann.energyGrid.node;
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      attConservativeMatrix = zeros(cellNumber);
      attNonConservativeMatrix = zeros(cellNumber);
      thresholdIsSmallerThanUmax = false;
      
      for gas = boltzmann.gasArray
        % loop over each collision with the gas
        for collision = gas.collisionArray
          threshold = collision.threshold;
          if strcmp(collision.type, 'Attachment') && threshold <= energyNode(end)
            % set flag for non-conservative algorithms
            thresholdIsSmallerThanUmax = true;
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1) + collision.crossSection(2:end));
            % local copy of target density
            density = collision.target.density;
            % evaluate numerical threshold
            numThreshold=floor(threshold/energyStep);
            
            % evaluation of nonconservative attachment collisional operator (always)
            for k=1:cellNumber
              attNonConservativeMatrix(k,k) = attNonConservativeMatrix(k,k) - ...
                density*energyCell(k)*cellCrossSection(k);
            end
            % avoid writing of conservative operator if threshold is smaller than the energy step
            if numThreshold==0
              continue
            end
            % evaluation of conservative attachment collisional operator (always needed)
            for k=1:cellNumber
              if (k<=cellNumber-numThreshold)
                attConservativeMatrix(k,k+numThreshold) = attConservativeMatrix(k,k+numThreshold) + ...
                  density*energyCell(k+numThreshold)*cellCrossSection(k+numThreshold);
              end
              attConservativeMatrix(k,k) = attConservativeMatrix(k,k) - density*energyCell(k)*cellCrossSection(k);
            end
            
          end
        end
      end
      
      % save matrixes in as properties of the boltzmann object
      boltzmann.attachmentConservativeMatrix = attConservativeMatrix;
      boltzmann.attachmentMatrix = attNonConservativeMatrix;
      boltzmann.attachmentThresholdIsSmallerThanUmax = thresholdIsSmallerThanUmax;
      
      if thresholdIsSmallerThanUmax
        boltzmann.includeNonConservativeAttachment = true;
      end
      
    end
    
    function eedf = invertLinearMatrix(boltzmann)
      % invertLinearMatrix invert the matrix of the discretized Boltzmann equation without considering non-linear
      % operators (non-conservative ionization or attachment, growth models, e-e collisions) in order to obtain an eedf
      
      % local copies of energy grid variables
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      
      % sum all contributions to the boltzmann matrix and reescale it to avoid very small numbers
      matrix = 1.e20*(boltzmann.fieldMatrix + boltzmann.elasticMatrix + boltzmann.CARMatrix + ...
        boltzmann.inelasticMatrix + boltzmann.ionizationConservativeMatrix + boltzmann.attachmentConservativeMatrix);
      
      % include normalization condition for the EEDF in the Boltzmann matrix
%       aux = matrixAux(1,:);
%       matrix(1,:) = matrix(1,:) + energyStep*sqrt(energyCell);
      matrix(1,:) = energyStep*sqrt(energyCell);
      
      % invert the Boltzmann matrix
      eedf = (matrix\([1 zeros(1,cellNumber-1)]'))';
      
      % restore proper nomalization for the EEDF
%       eedf(1) = -sum(aux(2:end).*eedf(2:end))/aux(1);
      
      % renormalize of the EEDF
      eedf = eedf/dot(eedf,sqrt(energyCell)*energyStep);
      
      % store eedf in the properties of the boltzmann object
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = iterativeSolution(boltzmann)
      % eedfTimeIntegration performs a time integration of the Boltzmann equation (to be completed when the function is
      % completely developed)
      
      % initial guess for the eedf from the linear Boltzmann equation (i.e. without non-linear operators)
      eedf = boltzmann.invertLinearMatrix();
      
      % iterative ode solver
      eedfTimeDerivative(0,eedf',boltzmann,true); % flush persistent variables previously stored in function memory
      odeOptionsLocal = boltzmann.odeOptions;
      odeOptionsLocal.NonNegative = 1:length(eedf);
      odeOptionsLocal.Events = @steadyStateEventFcn;
      [~,~,~,eedf,~] = ode15s(@eedfTimeDerivative,[0 inf],eedf',odeOptionsLocal,boltzmann,false);
      
      % select final solution if more than one is found
      if length(eedf(:,1))>1
        eedf = eedf(end,:);
      end
      
      % ensure proper normalization
      norm = sum(eedf.*sqrt(boltzmann.energyGrid.cell))*boltzmann.energyGrid.step;
      eedf = eedf/norm;
      
      % store eedf in the properties of the boltzmann object
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = mixingDirectSolutions(boltzmann)
      % mixingDirectSolutions solves the non-linear Boltzmann equation (i.e. including one of more of the following
      % operators: non-conservative ionization, non-conservative attachment or e-e collision) with an algorithm of
      % mixing direct solutions of a linearized Boltzmann equation.
      
      % local copies of energy grid variables
      cellNumber = boltzmann.energyGrid.cellNumber;
      
      % initial guess for the eedf from the linear Boltzmann equation (i.e. without non-linear operators)
      boltzmann.invertLinearMatrix();
      
      % select the appropiate method depending on the growth model (if non-conservative operators are activated)
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'spatial'
            solveEGrowthModel = str2func('solveSpatialGrowthMatrix');
            % reset spatial growth related matrices
            boltzmann.ionSpatialGrowthD  = zeros(cellNumber);
            boltzmann.ionSpatialGrowthU = zeros(cellNumber);
            boltzmann.fieldMatrixSpatGrowth = zeros(cellNumber);
          case 'temporal'
            solveEGrowthModel = str2func('solveTemporalGrowthMatrix');
            % reset temporal growth related matrices
            boltzmann.fieldMatrixTempGrowth = zeros(cellNumber);
            boltzmann.ionTemporalGrowth = zeros(cellNumber);
        end
      else
        solveEGrowthModel = [];
      end
      
      % solve the non-linear Boltzmann system to obtain an eedf (different cases depending on the activated operators)
      if boltzmann.includeEECollisions
        % reset e-e collisions related variables
        boltzmann.alphaee = 0;
        boltzmann.Aee = zeros(cellNumber);
        boltzmann.Bee = zeros(cellNumber);
        boltzmann.Afinal = zeros(1, cellNumber);
        boltzmann.Bfinal = zeros(1, cellNumber);
        if isempty(solveEGrowthModel)
          % solution with the e-e collisions only
          eedf = boltzmann.solveEEColl();
        else
          % solution with e-e collisions and growth model
          globalIteration = 0;        % counter for the number of iterations of the global cycle
          maxGlobalIterations = 20;   % maximum number of iterations for the global cycle
          while globalIteration < maxGlobalIterations
            eedfOld = solveEGrowthModel(boltzmann);   % cycle updating only the growth model operator
            eedf = boltzmann.solveEEColl();           % cycle updating only the e-e collision operator
            % exit global cycle in case of convergence (maximum relative difference of each value of the eedf < 1e-9)
            if max(abs(eedfOld-eedf)./eedfOld) < boltzmann.maxEedfRelError 
              break;
            end
            globalIteration = globalIteration+1;
            % throw a warning message if the global cycle has not converged
            if globalIteration == maxGlobalIterations
              warning('Global cycle (mixing of direct solutions): EEDF did not converge after %d iterations\n', ...
                maxGlobalIterations);
            end
          end
        end
      else
        % solution with growth model only 
        eedf = solveEGrowthModel(boltzmann);
      end
      
    end
    
    function eedf = solveSpatialGrowthMatrix(boltzmann)
      
      e = Constant.electronCharge;                % electron charge
      me = Constant.electronMass;                 % electron mass
      EoN = boltzmann.workCond.reducedFieldSI;    % reduced electric field (SI units)
      
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      energyNode = boltzmann.energyGrid.node;
      eedf = boltzmann.eedf;
      
      ionizationMatrixAux = boltzmann.ionizationMatrix;
      attachmentMatrixAux = boltzmann.attachmentMatrix;
      MatrixFieldSpatialGrowth = zeros(cellNumber);
      totalCrossSectionAux = boltzmann.totalCrossSection;
      cellTotalCrossSectionAux = 0.5*(totalCrossSectionAux(1:end-1) + totalCrossSectionAux(2:end));
      
      % writing of the Boltzmann matrix without the growth model and the electron-electron collisional operator
      baseMatrix = boltzmann.elasticMatrix + boltzmann.CARMatrix + boltzmann.fieldMatrix + ...
        boltzmann.inelasticMatrix + boltzmann.ionizationMatrix + boltzmann.attachmentMatrix;
      
      % electron-electron collisions terms
      Mee = zeros(cellNumber);
      if boltzmann.includeEECollisions
        alphaEE = boltzmann.alphaee;
        auxA=boltzmann.Aee;
        auxB=boltzmann.Bee;
        A = (alphaEE/energyStep)*(auxA*eedf');
        B = (alphaEE/energyStep)*(auxB*eedf');
        Mee(1:cellNumber+1:cellNumber*cellNumber) = -(A(1:cellNumber)+B(1:cellNumber));
        Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = B(2:cellNumber);
        Mee(2:cellNumber+1:cellNumber*cellNumber) = A(1:cellNumber-1);
      end
      
      mixingParam = boltzmann.mixingParameter;
      
      % calculation of the ionization rate
      integrandCI = sqrt(2*e/me)*energyStep*sum(ionizationMatrixAux+attachmentMatrixAux)';
      CIEffNew = dot(eedf,integrandCI);
      
      % mixing of solutions
      CIEffOld = CIEffNew/3;
      CIEffNew = mixingParam*CIEffNew + (1-mixingParam)*CIEffOld;
      
      % writting of the diffusion and mobility components of the spatial growth terms
      D = zeros(cellNumber);
      U = zeros(cellNumber);
      D0 = energyCell./(3*cellTotalCrossSectionAux);
      U0sup = EoN/(6*energyStep)*[0 energyCell(1:cellNumber-1)./(cellTotalCrossSectionAux(1:cellNumber-1))];
      U0inf = -EoN/(6*energyStep)*[energyCell(2:cellNumber)./(cellTotalCrossSectionAux(2:cellNumber)) 0];
      U0 = U0sup+U0inf;
      
      % calculation of the gas density times the diffusion coefficient using the matrix of the diffusion component
      %of the spatial growth model
      ND = sqrt(2*e/me)*energyStep*dot(D0,eedf);
      
      % calculation of the mobility times the electric field using the matrix of the diffusion component of the
      %spatial growth model
      muE = -sqrt(2*e/me)*energyStep*dot(U0,eedf);
      
      % initial eedf guess may lead to negative root argument, in which case, the reduced townsend coefficient
      %should be calculated assuming that there is no electron density gradient
      if muE^2-4*CIEffNew*ND < 0
        alphaRedEffNew = CIEffNew/muE;
      else
        alphaRedEffNew = (muE - sqrt(muE^2-4*CIEffNew*ND))/(2*ND);
      end
      
      iter=0;
      convergence=0;
      while(convergence==0)
        % writing of MatrixI which refers to the additional electric field terms of the spatial growth model
        g_fieldSpatialGrowthAux = alphaRedEffNew*EoN*energyNode./(6*totalCrossSectionAux);
        g_fieldSpatialGrowthAux(1) = 0;
        g_fieldSpatialGrowthAux(end) = 0;
        for k=1:cellNumber
          MatrixFieldSpatialGrowth(k,k) = -(g_fieldSpatialGrowthAux(k) - g_fieldSpatialGrowthAux(k+1))/energyStep;
          if k>1
            MatrixFieldSpatialGrowth(k,k-1) = -g_fieldSpatialGrowthAux(k)/energyStep;
          end
          if k<cellNumber
            MatrixFieldSpatialGrowth(k,k+1) = g_fieldSpatialGrowthAux(k+1)/energyStep;
          end
        end
        
        % calculation of the diffusion and mobility matrices of the spatial growth model
        D(1:cellNumber+1:cellNumber*cellNumber) = alphaRedEffNew*alphaRedEffNew*D0;
        U(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = alphaRedEffNew*U0sup(2:cellNumber);
        U(2:cellNumber+1:cellNumber*cellNumber) = alphaRedEffNew*U0inf(1:cellNumber-1);
        
        % writting of the Boltzmann matrix (expansion of the following expression because of performance reasons)
        % matrixAux =  1.e20*(MatrixI + D + U + baseMatrix + Mee)
        matrixAux = 1.e20*baseMatrix;
        matrixAux(1:cellNumber+1:cellNumber*cellNumber) = matrixAux(1:cellNumber+1:cellNumber*cellNumber) + ...
          1.e20*(MatrixFieldSpatialGrowth(1:cellNumber+1:cellNumber*cellNumber) + ...
          D(1:cellNumber+1:cellNumber*cellNumber) + Mee(1:cellNumber+1:cellNumber*cellNumber));
        matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = ...
          matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldSpatialGrowth(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + ...
          U(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) );
        matrixAux(2:cellNumber+1:cellNumber*cellNumber) = matrixAux(2:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldSpatialGrowth(2:cellNumber+1:cellNumber*cellNumber) + ...
          U(2:cellNumber+1:cellNumber*cellNumber) + Mee(2:cellNumber+1:cellNumber*cellNumber) );
        
        % include normalization condition for the EEDF in the Boltzmann matrix
        aux = matrixAux(1,:);
        matrixAux(1,:) = matrixAux(1,:) + energyStep*sqrt(energyCell);
        
        % invert the Boltzmann matrix
        eedfNew = eedf;
        eedf = (matrixAux\([1 zeros(1,cellNumber-1)]'))';
        
        % restore proper nomalization for the EEDF
        eedf(1) = -sum(aux(2:end).*eedf(2:end))/aux(1);
        
        % renormalize of the EEDF
        norm = sum(sqrt(energyCell).*eedf*energyStep);
        eedf = eedf/norm;
        
        % calculation of the ionization rate
        CIEffOld = CIEffNew;
        CIEffNew = dot(eedf,integrandCI);
        % mixing of solutions
        CIEffNew = mixingParam*CIEffNew + (1-mixingParam)*CIEffOld;
        
        % calculation of the gas density times diffusion coefficient, and the mobility coefficient times the
        % electric field
        ND = sqrt(2*e/me)*energyStep*dot(D0,eedf);
        muE = -sqrt(2*e/me)*energyStep*dot(U0,eedf);
        
        % calculation of the new effective reduced first Townsend coefficient
        alphaRedEffOld = alphaRedEffNew;
        if muE^2-4*CIEffNew*ND < 0
          alphaRedEffNew = CIEffNew/muE;
        else
          alphaRedEffNew = (muE - sqrt(muE^2-4*CIEffNew*ND))/(2*ND);
        end
        % mixing of solutions
        alphaRedEffNew =(alphaRedEffNew*mixingParam + (1-mixingParam)*alphaRedEffOld);
        
        % testing of convergence criteria
        if (alphaRedEffNew==0 || (alphaRedEffNew-alphaRedEffOld)/alphaRedEffOld <1e-10) &&  ...
            max(abs(eedf - eedfNew)./eedf) < boltzmann.maxEedfRelError || iter>150
          convergence=1;
          if iter > 150 && ~boltzmann.includeEECollisions
            warning('Spatial growth iterative scheme did not converge\n');
          end
        end
        
        iter = iter+1;
      end
      
      % copy of terms used in power balance and swarm parameters calculation
      boltzmann.alphaRedEff = alphaRedEffOld;
      boltzmann.CIEff = CIEffOld;
      boltzmann.g_fieldSpatialGrowth = g_fieldSpatialGrowthAux;
      boltzmann.eedf = eedf;
      
      % copy of spatial growth model terms used on the electron-electron collisions routine
      if boltzmann.includeEECollisions
        boltzmann.ionSpatialGrowthD = D;
        boltzmann.ionSpatialGrowthU = U;
        boltzmann.fieldMatrixSpatGrowth = MatrixFieldSpatialGrowth;
      end
      
      % saving eedf
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = solveTemporalGrowthMatrix(boltzmann)
      
      e = Constant.electronCharge;                % electron charge
      me = Constant.electronMass;                 % electron mass
      EoN = boltzmann.workCond.reducedFieldSI;    % reduced electric field (SI units)
      WoN = boltzmann.workCond.reducedExcFreqSI;  % reduced angular exitation frequency (SI units)
      
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      energyNode = boltzmann.energyGrid.node;
      eedf = boltzmann.eedf;
      
      ionizationMatrixAux = boltzmann.ionizationMatrix;
      attachmentMatrixAux = boltzmann.attachmentMatrix;
      totalCrossSectionAux = boltzmann.totalCrossSection;
      totalCSI = zeros(size(totalCrossSectionAux));
      g_fieldTemporalGrowthAux = zeros(size(boltzmann.g_E));
      MatrixFieldTemporalGrowth = zeros(cellNumber);
      growthMatrix = zeros(cellNumber);
      
      % writing of the Boltzmann matrix without the growth model and the electron-electron collisional operator
      baseMatrix = boltzmann.elasticMatrix + boltzmann.CARMatrix + boltzmann.inelasticMatrix + ...
        ionizationMatrixAux + attachmentMatrixAux;
      
      % electron-electron collisions terms (not updated in the density growth operator)
      Mee = zeros(cellNumber);
      if boltzmann.includeEECollisions
        alphaEE = boltzmann.alphaee;
        auxA=boltzmann.Aee;
        auxB=boltzmann.Bee;
        A = (alphaEE/energyStep)*(auxA*eedf');
        B = (alphaEE/energyStep)*(auxB*eedf');
        Mee(1:cellNumber+1:cellNumber*cellNumber) = -(A(1:cellNumber)+B(1:cellNumber));
        Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = B(2:cellNumber);
        Mee(2:cellNumber+1:cellNumber*cellNumber) = A(1:cellNumber-1);
      end
      
      mixingParam = boltzmann.mixingParameter;
      
      % calculation of the ionization rate
      integrandCI = sqrt(2*e/me)*energyStep*sum(ionizationMatrixAux+attachmentMatrixAux)';
      CIEffNew = dot(eedf,integrandCI);
      
      % mixing of solutions
      CIEffOld = CIEffNew/3;
      CIEffNew = mixingParam*CIEffNew + (1-mixingParam)*CIEffOld;
      
      iter=0;
      convergence=0;
      while(convergence==0)
        % writing the total cross section plus the ionizatino rate divided by the electron velocity
        totalCSI(1)= totalCrossSectionAux(1);
        for j=2:cellNumber+1
          totalCSI(j) =  totalCrossSectionAux(j) + CIEffNew*sqrt(me./(2*e*(j-1)*energyStep));
        end
        
        % writing of the MatrixI which refers to the electric field term of the temporal growth model and
        % growthMatrix that refers to the time variation term (dn/dt) of the temporal grwoth model
        g_fieldTemporalGrowthAux = EoN^2*energyNode./(3*totalCSI.*(1+me*WoN^2./(2*e*energyNode.*totalCSI.*totalCSI)));
        g_fieldTemporalGrowthAux(1) = 0;
        g_fieldTemporalGrowthAux(end) = 0;
        for k = 1:cellNumber
          MatrixFieldTemporalGrowth(k,k) = -(g_fieldTemporalGrowthAux(k)+g_fieldTemporalGrowthAux(k+1))/energyStep^2;
          if k>1
            MatrixFieldTemporalGrowth(k,k-1) = g_fieldTemporalGrowthAux(k)/energyStep^2;
          end
          if k<cellNumber
            MatrixFieldTemporalGrowth(k,k+1) = g_fieldTemporalGrowthAux(k+1)/energyStep^2;
          end
          growthMatrix(k,k) = -CIEffNew*sqrt(me/(2*e))*sqrt(energyCell(k));
        end
        
        % writting of the Boltzmann matrix (expansion of the following expression because of performance reasons)
%         matrixAux =  1.e20*(growthMatrix + MatrixI + baseMatrix + Mee);
        matrixAux = 1.e20*baseMatrix;
        matrixAux(1:cellNumber+1:cellNumber*cellNumber) = matrixAux(1:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldTemporalGrowth(1:cellNumber+1:cellNumber*cellNumber) + ...
          growthMatrix(1:cellNumber+1:cellNumber*cellNumber) + Mee(1:cellNumber+1:cellNumber*cellNumber) );
        matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = ...
          matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldTemporalGrowth(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + ...
          Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) );
        matrixAux(2:cellNumber+1:cellNumber*cellNumber) = matrixAux(2:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldTemporalGrowth(2:cellNumber+1:cellNumber*cellNumber) + ...
          Mee(2:cellNumber+1:cellNumber*cellNumber) );
        
        % include normalization condition for the EEDF in the Boltzmann matrix
        aux = matrixAux(1,:);
        matrixAux(1,:) = matrixAux(1,:) + energyStep*sqrt(energyCell);
        
        % invert the Boltzmann matrix
        eedfNew = eedf;
        eedf = (matrixAux\([1 zeros(1,cellNumber-1)]'))';
        
        % restore proper nomalization for the EEDF
        eedf(1) = -sum(aux(2:end).*eedf(2:end))/aux(1);
        
        % renormalize of the EEDF
        norm = sum(sqrt(energyCell).*eedf*energyStep);
        eedf = eedf/norm;
        
        % calculation of the ionization rate
        CIEffOld = CIEffNew;
        CIEffNew = dot(eedf,integrandCI);
        % mixing of solutions
        CIEffNew = mixingParam*CIEffNew + (1-mixingParam)*CIEffOld;
        
        % testing of convergence criteria
        if (CIEffNew==0 || abs((CIEffNew-CIEffOld)/CIEffOld)<10e-10) &&  ...
            max(abs(eedf-eedfNew)./eedf) < boltzmann.maxEedfRelError || iter>150
          convergence=1;
          if iter > 150 && ~boltzmann.includeEECollisions
            warning('Temporal growth iterative scheme did not converge\n');
          end
        end
        
        iter = iter+1;
      end
      
      % copy of terms used in power balance and swarm parameters calculation
      boltzmann.CIEff = CIEffOld;
      boltzmann.g_fieldTemporalGrowth = g_fieldTemporalGrowthAux;
      boltzmann.eedf = eedf;
      
      % copy of temporal growth model terms used on the electron-electron collisions routine
      if boltzmann.includeEECollisions
        boltzmann.ionTemporalGrowth = growthMatrix;
        boltzmann.fieldMatrixTempGrowth = MatrixFieldTemporalGrowth;
      end
      
      % saving eedf
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = solveEEColl(boltzmann)
      
      e = Constant.electronCharge;                % electron charge
      e0 = Constant.vacuumPermittivity;           % vacuum permitivity (SI units)
      ne = boltzmann.workCond.electronDensity;    % electron density (SI units)
      n0 = boltzmann.workCond.gasDensity;         % gas density (SI units)
      EoN = boltzmann.workCond.reducedField;      % reduced electric field (Td)
      
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      energyNode = boltzmann.energyGrid.node;
      eedf = boltzmann.eedf;
      
      Mee = zeros(cellNumber);
      
      % writting of Boltzmann matrix without the electron-electron collisional operator
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'spatial'
            baseMatrix = boltzmann.ionizationMatrix + boltzmann.attachmentMatrix + boltzmann.elasticMatrix + ...
              boltzmann.CARMatrix + boltzmann.inelasticMatrix + boltzmann.fieldMatrix + boltzmann.ionSpatialGrowthD + ...
              boltzmann.ionSpatialGrowthU +  boltzmann.fieldMatrixSpatGrowth;
          case 'temporal'
            baseMatrix = boltzmann.ionizationMatrix + boltzmann.attachmentMatrix + boltzmann.elasticMatrix + ...
              boltzmann.CARMatrix + boltzmann.inelasticMatrix + boltzmann.fieldMatrixTempGrowth + ...
              boltzmann.ionTemporalGrowth;
        end
      else
        baseMatrix = boltzmann.ionizationConservativeMatrix + boltzmann.attachmentConservativeMatrix + ...
          boltzmann.elasticMatrix + boltzmann.CARMatrix + boltzmann.inelasticMatrix + boltzmann.fieldMatrix ;
      end
      
      % writing of auxiliary matrix A, without constant alpha
      auxA = zeros(cellNumber,cellNumber);
      auxEnergyArray = -(energyStep/2)*sqrt(energyCell)+(2/3)*energyCell.^(3/2);
      % from power conservation, terms on the last row (cellNumber,:) and on the first column (:,1) are zero
      for k=1:cellNumber-1
        auxA(k,2:k) = auxEnergyArray(2:k);
        auxA(k,k+1:cellNumber) = (2/3)*energyNode(k+1)^(3/2);
      end
      
      % detailed balance condition
      for k=1:cellNumber-1
        auxA(k,2:cellNumber) = sqrt(auxA(k,2:cellNumber).*auxA(1:cellNumber-1,k+1)'); 
      end
      
      % writing of auxiliary matrix B, without constant alpha
      auxB = transpose(auxA);
      
      % Calculation of constant alpha
      % mean energy
      meanEnergy = energyStep*dot(energyCell.^(3/2),eedf);
      % electron temperature in eV Te = (2/3)*meanEnergy
      Te = (2/3)*meanEnergy;
      % Coulomb logarithm
      logC = log(12*pi*(e0*Te/e)^(3/2)/sqrt(ne));
      % alpha
      alpha = (ne/n0)*(e^2/(8*pi*e0^2))*logC;
      
      % initial estimate ratio=PE/Pee set to zero
      ratioNew=0;
      % save of eedf calculated without electron-electron collisional operator
      eedfNew=eedf;
      
      % cycle variables set to zero
      conve=0;
      iter=0;
      
      while conve==0
        % calculation of electron-electron collisions vectors of upflux (A) and downflux (B)
        A = (alpha/energyStep)*(auxA*eedf');
        B = (alpha/energyStep)*(auxB*eedf');
        
        % writing of the electron-electron operator
        Mee(1:cellNumber+1:cellNumber*cellNumber) = -(A(1:cellNumber)+B(1:cellNumber));
        Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = B(2:cellNumber);
        Mee(2:cellNumber+1:cellNumber*cellNumber) = A(1:cellNumber-1);
        
        % sum all contributions to the Boltzmann matrix and reescale it to avoid very small numbers
        matrixAux = 1e20*(baseMatrix + Mee);
        
        % include normalization condition for the EEDF in the Boltzmann matrix
        aux = matrixAux(1,:);
        matrixAux(1,:) = matrixAux(1,:) + energyStep*sqrt(energyCell);
        
        % invert the Boltzmann matrix
        eedf = (matrixAux\([1 zeros(1,cellNumber-1)]'))';
        
        % restore proper nomalization for the EEDF
        eedf(1) = -sum(aux(2:end).*eedf(2:end))/aux(1);
        norm = sum(sqrt(energyCell(1:cellNumber)).*eedf(1:cellNumber))*energyStep;
        eedf = eedf/norm;
        
        % calculation of ratio =PE/Pee
        % field power
        boltzmann.eedf = eedf;
        boltzmann.Afinal = A;
        boltzmann.Bfinal = B;
        boltzmann.evaluatePower(false);
        Preference = boltzmann.power.reference;
        Pee = boltzmann.power.electronElectron;
        ratio = abs(Pee/Preference);
        
        % test convergence parameters
        if  max(abs(eedf-eedfNew)./eedf) < boltzmann.maxEedfRelError
          if abs(ratio)<1e-9
            conve = 1;
          end
          if abs(ratio)>1e-9 && iter>200
            fprintf('e-e iterative scheme: EEDF has converged but the abs(ratio)=abs(Pee/Pref)= %.16g > 1e-9\n',ratio);
            fprintf('Ionization degree:%f \t Reduced electric field:%f Td \n',ne/n0,EoN);
            conve = 1;
          end
        elseif iter == 300 && isempty(boltzmann.eDensGrowthModel)
          conve = 1;
          warning('Electron-electron iterative scheme: EEDF did not converge\n');
        end
        
        % acceleration scheme based on the Newton-Raphson
        if iter>0 && conve~=1
          
          % saving quantities calculated from the inversion of the Boltzmann matrix
          ratioOld = ratioNew;
          ratioNew=ratio;
          eedfOld = eedfNew;
          eedfNew = eedf;
          
          % eedf estimate
          eedf = eedfNew - (eedfNew-eedfOld)*ratioNew/(ratioNew-ratioOld);
          
          % correction of the eedf estimate
          for i=1:cellNumber
            if(abs(norm*eedf(i)) < eps)
              eedf(i) = abs(eedf(i));
            elseif((eedf(i)) < 0 && abs(eedf(i)))
              eedf(i) = abs(eedf(i));
            end
          end
        end
        
        % calculation of constant alpha
        meanEnergy = energyStep*dot(energyCell.^(3/2),eedf);
        % electron temperature
        Te = (2/3)*meanEnergy;
        % Coulomb logarithm
        logC = log(12*pi*(e0*Te/e)^(3/2)/sqrt(ne));
        alpha = (ne/n0)*(e^2/(8*pi*e0^2))*logC;
        
        iter = iter+1;
      end
      
      % saving auxiliary matrix to be used on the ionization routine
      if ~strcmp(boltzmann.ionCollOpType,'conservative')
        boltzmann.Aee = auxA;
        boltzmann.Bee = auxB;
        boltzmann.alphaee = alpha;
      end
      
      % saving eedf
      boltzmann.eedf = eedf;
      
    end
    
    function power = evaluatePower(boltzmann, isFinalSolution)
      
      % initialize power structure
      power = struct('field', 0, 'elasticNet', 0, 'elasticGain', 0, 'elasticLoss', 0, 'carNet', 0, 'carGain', 0, ...
        'carLoss', 0, 'excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
        'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
        'ionizationIne', 0, 'attachmentIne', 0, 'inelastic', 0, 'superelastic', 0, 'eDensGrowth', 0, ...
        'electronElectron', 0, 'gases', '');
      
      % save a local copy of the EEDF because of performance reasons
      eedfLocal = boltzmann.eedf;
      
      % save a local copy of the energy grid information
      energyGridLocal = boltzmann.energyGrid;   % energyGrid object
      N = energyGridLocal.cellNumber;           % number of cells in the energy grid
      energyStep = energyGridLocal.step;        % energy step
      energyCell = energyGridLocal.cell;        % value of the energy at the cells
      energyNode = energyGridLocal.node;        % value of the energy at the nodes
      
      % multiplicative constant to obtain the right units
      factor = sqrt(2*Constant.electronCharge/Constant.electronMass);
      
      % auxiliary quantities needed to evaluate the elastic and CAR powers
      kTg = Constant.boltzmannInEV*boltzmann.workCond.gasTemperature;
      aux1 = kTg+energyStep*0.5;
      aux2 = kTg-energyStep*0.5;
      
      % evaluate power absorved per electron at unit gas density due to elastic collisions
      g_cLocal = boltzmann.g_c; % elements of the elastic collision operator (local copy)
      power.elasticNet = factor*sum(eedfLocal.*(g_cLocal(2:end)*aux2-g_cLocal(1:end-1)*aux1));
      power.elasticGain = factor*kTg*sum(eedfLocal.*(g_cLocal(2:end)-g_cLocal(1:end-1)));
      power.elasticLoss = power.elasticNet-power.elasticGain;
      
      % evaluate power absorved per electron at unit gas density due to rotations CAR
      if ~isempty(boltzmann.CARgases)
        g_carLocal = boltzmann.g_car; % elements of the CAR operator (local copy)
        power.carNet = factor*sum(eedfLocal.*(g_carLocal(2:end)*aux2-g_carLocal(1:end-1)*aux1));
        power.carGain = factor*kTg*sum(eedfLocal.*(g_carLocal(2:end)-g_carLocal(1:end-1)));
        power.carLoss = power.carNet-power.carGain;
      end
      
      % evaluate power gained from electric field and lossed due to electron density growth terms
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'temporal'
            g_ELocal = boltzmann.g_fieldTemporalGrowth; % elements of the electric field operator (local copy)
            power.field = factor*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
            
            power.eDensGrowth = - boltzmann.CIEff*energyStep*sum(eedfLocal.*energyCell.*sqrt(energyCell));
          case 'spatial'
            g_ELocal = boltzmann.g_E; % elements of the electric field operator (local copy)
            power.field = factor*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
            % evaluate correction due to the electric field term of the spatial growth model
            g_ELocal = boltzmann.g_fieldSpatialGrowth;
            correction = factor*sum(eedfLocal.*(-g_ELocal(2:end)-g_ELocal(1:end-1)))*energyStep;
            power.field = power.field+correction;
            
            totalCrossSectionLocal = boltzmann.totalCrossSection;
            alphaRedEffLocal = boltzmann.alphaRedEff;
            cellTotalCrossSection =  0.5*(totalCrossSectionLocal(1:N) + totalCrossSectionLocal(2:N+1));
            
            % diffusion contribution
            powerDiffusion = alphaRedEffLocal^2*factor*energyStep/3*sum(energyCell(1:N).^2.*eedfLocal(1:N)./...
              cellTotalCrossSection(1:N));
            
            % mobility contribution
            powerMobility = factor*alphaRedEffLocal*(boltzmann.workCond.reducedFieldSI/6)*(...
              energyCell(1)^2*eedfLocal(2)/cellTotalCrossSection(1) - ...
              energyCell(N)^2*eedfLocal(N-1)/cellTotalCrossSection(N) + ...
              sum(energyCell(2:N-1).^2.*(eedfLocal(3:N)-eedfLocal(1:N-2))./cellTotalCrossSection(2:N-1)));
            
            % power of spatial growth component
            power.eDensGrowth =  powerDiffusion + powerMobility;
        end
      else
        g_ELocal = boltzmann.g_E; % elements of the electric field operator (local copy)
        power.field = factor*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
      end
      
      % power absorved per electron at unit gas density due to electron-electron collisions (to be revised)
      if boltzmann.includeEECollisions
        A = boltzmann.Afinal;
        B = boltzmann.Bfinal;
        power.electronElectron = -factor*sum((A(1:N)-B(1:N)).*eedfLocal(1:N)')*energyStep^2;
      end
      
      % evaluate power absorved per electron at unit gas density due to the inelastic/super-elastic collisions
      % loop over each gas in the mixture
      for gas = boltzmann.gasArray
        gasName = gas.name;
        % initialize power balance information of this gas
        gasPower = struct('excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
          'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
          'ionizationIne', 0, 'attachmentIne', 0);
        % loop over each collision with the gas
        for collision = gas.collisionArray
          % collision type
          collType = collision.type;
          % avoid Effective or Elastic collisions and collisions which threshold is larger than the maximum energy
          if strcmp(collType, 'Effective') || strcmp(collType, 'Elastic') || collision.threshold > energyNode(end)
            continue;
          elseif strcmp(collision.type, 'Ionization') && ~strcmp(boltzmann.ionCollOpType,'conservative')
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
            
            threshold = collision.threshold;
            lmin=floor(threshold/energyStep);
            
            switch boltzmann.ionCollOpType
              case 'equalSharing'
                ionizationIneAux = -factor*collision.target.density*energyStep*(...
                  sum(energyCell(lmin:N).^2.*cellCrossSection(lmin:N).*eedfLocal(lmin:N)) + ...
                  2*energyCell(lmin+1)*sum(energyCell(2+lmin:2:N).*cellCrossSection(2+lmin:2:N).*...
                  eedfLocal(2+lmin:2:N))-2*sum(energyCell(2+lmin:2:N).^2.*cellCrossSection(2+lmin:2:N).*...
                  eedfLocal(2+lmin:2:N)));
                
              case 'oneTakesAll'
                ionizationIneAux = -factor*collision.target.density*energyStep*energyCell(lmin)*...
                  sum(eedfLocal(lmin:N).*energyCell(lmin:N).*cellCrossSection(lmin:N));
                
              case 'usingSDCS'
                %calculation of the total (integrated) ionization cross section using the SDCS
                TICS = zeros(size(eedfLocal));
                W = gas.OPBParameter;
                if isempty(W)
                  W = threshold;
                end
                auxArray = 1./(1+(energyCell(1:N)/W).^2);
                for k=2:N
                  auxArray(k) = auxArray(k)+auxArray(k-1);
                  kmax = floor((k-lmin)/2);
                  if kmax>0
                    TICS(k) = TICS(k) + cellCrossSection(k)/(W*atan((energyCell(k)-threshold)/(2*W)))*auxArray(kmax);
                  end
                end
                TICS = energyStep*TICS;
                ionizationIneAux = -factor*collision.target.density*energyCell(lmin+1)*energyStep*...
                  sum(eedfLocal(1:N).*energyCell(1:N).*TICS(1:N));
            end
            gasPower.ionizationIne = gasPower.ionizationIne+ionizationIneAux;
            continue;
          elseif strcmp(collision.type, 'Attachment') && ~strcmp(boltzmann.ionCollOpType,'conservative')
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
            threshold = collision.threshold;
            lmin=floor(threshold/energyStep);
            gasPower.attachmentIne = gasPower.attachmentIne - factor*collision.target.density*energyStep*...
              sum(eedfLocal(1+lmin:N).*energyCell(1+lmin:N).^2.*cellCrossSection(1+lmin:N));
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
      
      % evaluate power balance
      powerValues = [power.field power.elasticGain power.elasticLoss power.carGain power.carLoss ...
        power.excitationSup power.excitationIne power.vibrationalSup power.vibrationalIne ...
        power.rotationalSup power.rotationalIne power.eDensGrowth power.electronElectron];
      totalGain = 0;
      totalLoss = 0;
      for powerValue = powerValues
        if powerValue > 0 
          totalGain = totalGain+powerValue;
        else
          totalLoss = totalLoss+powerValue;
        end
      end
      power.balance = power.field + power.elasticNet + power.carNet + power.inelastic + power.superelastic + ...
        power.eDensGrowth + power.electronElectron;
      power.relativeBalance = abs(power.balance)/totalGain;
      power.reference = totalGain;
      
      % store power balance information in the boltzmann properties
      boltzmann.power = power;
      
      % check for errors in the power balance for final solution
      if isFinalSolution && power.relativeBalance > boltzmann.maxPowerBalanceRelError
        warning(sprintf(['Relative power balance greater than %e.\n' ...
          'Results may be wrong, please check input/output of the simulation'], boltzmann.maxPowerBalanceRelError));
      end
      
    end
    
    function swarmParam = evaluateSwarmParameters(boltzmann)
      
      % save local copies of different constants
      me = Constant.electronMass;
      e = Constant.electronCharge;
      
      % save a local copy of the energy grid information
      energyNode = boltzmann.energyGrid.node;
      
      % save local copy of total momentum transfer cross section
      totalCrossSectionAux = boltzmann.totalCrossSection;
      if strcmp(boltzmann.eDensGrowthModel,'temporal') && ...
          ( boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment )
        totalCrossSectionAux(2:end) = totalCrossSectionAux(2:end) + boltzmann.CIEff*sqrt(me./(2*e*energyNode(2:end)));
      end
      
      % initialize transport parameters structure
      swarmParam = struct('redDiffCoeff', [], 'redMobCoeff', [], 'redTownsendCoeff', [], 'redAttCoeff', [],...
        'meanEnergy', [], 'characEnergy', [], 'Te', [], 'driftVelocity', []);
      
      % evaluate reduced diffusion coefficient
      swarmParam.redDiffCoeff = (2.0/3.0)*sqrt(2.0*Constant.electronCharge/Constant.electronMass)*...
        sum(boltzmann.energyGrid.cell.*boltzmann.eedf./(totalCrossSectionAux(1:end-1)+...
        totalCrossSectionAux(2:end)))*boltzmann.energyGrid.step;
      
      % evaluate reduced mobility coefficient
      swarmParam.redMobCoeff = -sqrt(2.0*Constant.electronCharge/Constant.electronMass)/3.0*...
        sum(boltzmann.energyGrid.node(2:end-1).*(boltzmann.eedf(2:end)-boltzmann.eedf(1:end-1))./...
        totalCrossSectionAux(2:end-1));
      
      % evaluate drift velocity
      if strcmp(boltzmann.eDensGrowthModel,'spatial') && ...
          ( boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment )
        swarmParam.driftVelocity = - swarmParam.redDiffCoeff*boltzmann.alphaRedEff + ...
          swarmParam.redMobCoeff*boltzmann.workCond.reducedFieldSI;
      else
        swarmParam.driftVelocity = swarmParam.redMobCoeff*boltzmann.workCond.reducedFieldSI;
      end
      
      % evaluate reduced Townsend coefficient
      totalIonRateCoeff = 0;
      for gas = boltzmann.gasArray
        for collision = gas.collisionArray
          if strcmp(collision.type, 'Ionization')
            totalIonRateCoeff = totalIonRateCoeff + collision.target.density*collision.ineRateCoeff;
            break;
          end
        end
      end
      swarmParam.redTownsendCoeff = totalIonRateCoeff / swarmParam.driftVelocity;
      
      % evaluate reduced attachment coefficient
      totalAttRateCoeff = 0;
      for gas = boltzmann.gasArray
        for collision = gas.collisionArray
          if strcmp(collision.type, 'Attachment')
            totalAttRateCoeff = totalAttRateCoeff + collision.target.density*collision.ineRateCoeff;
            break;
          end
        end
      end
      swarmParam.redAttCoeff = totalAttRateCoeff / swarmParam.driftVelocity;
      
      % evaluate mean energy
      swarmParam.meanEnergy = sum(boltzmann.energyGrid.cell(1:end).^(1.5).*boltzmann.eedf(1:end))*...
        boltzmann.energyGrid.step;
      
      % evaluate characteristic energy
      swarmParam.characEnergy = swarmParam.redDiffCoeff/swarmParam.redMobCoeff;
      
      % evaluate electron temperature
      swarmParam.Te = (2./3.)*swarmParam.meanEnergy;
      boltzmann.workCond.update('electronTemperature', swarmParam.Te);
      
      % store swarm parameters information in the boltzmann properties
      boltzmann.swarmParam = swarmParam;
      
    end
    
    function [rateCoeffAll, rateCoeffExtra] = evaluateRateCoeff(boltzmann)
      
      % initialize rateCoeffAll and rateCoeffExtra structures
      rateCoeffAll = struct.empty;
      rateCoeffExtra = struct.empty;
      
      % evaluate rate coefficient for all collision 
      for gas = boltzmann.gasArray
        % collisions taken into account for solving the eedf
        for collision = gas.collisionArray
          rateCoeffAll(end+1).collID = collision.ID;
          [ineRate, supRate] = collision.evaluateRateCoeff(boltzmann.eedf);
          rateCoeffAll(end).value = [ineRate, supRate];
          rateCoeffAll(end).collDescription = collision.description;
        end
        % collisions not taken into account for solving the eedf
        for collision = gas.collisionArrayExtra
          rateCoeffExtra(end+1).collID = collision.ID;
          [ineRate, supRate] = collision.evaluateRateCoeff(boltzmann.eedf);
          rateCoeffExtra(end).value = [ineRate, supRate];
          rateCoeffExtra(end).collDescription = collision.description;
        end
      end
      
      % store rate coefficients information in the boltzmann properties
      boltzmann.rateCoeffAll = rateCoeffAll;
      boltzmann.rateCoeffExtra = rateCoeffExtra;
      
    end
    
    function evaluateFirstAnisotropy(boltzmann)
    % evaluateFirstAnisotropy evaluates the value of the first anisotropy of the electron distribution function in
    % energy space in the framework of the two term approximation.
      
      % local copy of variables
      localEedf = boltzmann.eedf;                             % electron energy distribution function (isotropic)
      energyCell = boltzmann.energyGrid.cell;                 % values of energy at cell position (same as eedf)
      energyStep = boltzmann.energyGrid.step;                 % energy step of the energy grid
      localTotalCrossSection = boltzmann.totalCrossSection;   % total momentum transfer cross section
      e = Constant.electronCharge;                            % electron charge
      me = Constant.electronMass;                             % electron mass
      EoN = boltzmann.workCond.reducedFieldSI;                % reduced electric field (SI units)
      WoN = boltzmann.workCond.reducedExcFreqSI;              % reduced angular exitation frequency (SI units)
      
      % evaluate derivative of the eedf
      eedfDerivative = zeros(size(localEedf));
      eedfDerivative(1) = (localEedf(2)-localEedf(1))/energyStep;                     % 1st order forward approximation
      eedfDerivative(end) = (localEedf(end)-localEedf(end-1))/energyStep;             % 1st order backward approximation
      eedfDerivative(2:end-1) = (localEedf(3:end)-localEedf(1:end-2))/(energyStep*2); % 2nd order centered approximation
      
      % evaluate total momentum transfer cross section at cell positions
      totalCrossSectionCell = (localTotalCrossSection(1:end-1)+localTotalCrossSection(2:end))/2.0;
      
      % evaluate the first anisotropy
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'temporal'
            totalCrossSectionCell = totalCrossSectionCell + boltzmann.CIEff./(sqrt(2*e*energyCell/me));
            if WoN == 0
              boltzmann.firstAnisotropy = -EoN.*eedfDerivative./totalCrossSectionCell;
            else
              boltzmann.firstAnisotropy = -EoN*sqrt(2).*eedfDerivative./(totalCrossSectionCell+WoN^2*me./...
                (2.0*e*energyCell.*totalCrossSectionCell));
            end
          case 'spatial'
            boltzmann.firstAnisotropy = -(boltzmann.alphaRedEff.*localEedf+EoN.*eedfDerivative)./totalCrossSectionCell;
        end
      elseif WoN == 0
        boltzmann.firstAnisotropy = -EoN.*eedfDerivative./totalCrossSectionCell;
      else
        boltzmann.firstAnisotropy = -EoN*sqrt(2).*eedfDerivative./(totalCrossSectionCell+WoN^2*me./...
          (2.0*e*energyCell.*totalCrossSectionCell));
      end
      
    end
    
  end
  
end

function dfdt = eedfTimeDerivative(~,eedf,boltzmann,clearPersistentVars)
% eedfTimeDerivative evaluates the time derivative of the EEDF at a given time t, assuming that the EEDF at this precise
% time is given by eedf. 
  
  % local copy of fundamental constants
  persistent e;
  persistent me;
  persistent e0;
  if isempty(e)
    e = Constant.electronCharge;                % electron charge (SI units)
    me = Constant.electronMass;                 % electron mass (SI units)
    e0 = Constant.vacuumPermittivity;           % vacuum permitivity (SI units)
  end
  
  % local copy of other simulation constants and values
  persistent cellNumber;
  persistent energyStep;
  persistent energyCell;
  persistent energyNode;
  persistent ne;
  persistent N;
  persistent EoN;
  persistent WoN;
  if isempty(cellNumber)
    cellNumber = boltzmann.energyGrid.cellNumber;   % number of cells used in the energy grid
    energyStep = boltzmann.energyGrid.step;         % energy step of the energy grid
    energyCell = boltzmann.energyGrid.cell;         % energy at cells of the energy grid
    energyNode = boltzmann.energyGrid.node;         % energy at nodes of the energy grid
    ne = boltzmann.workCond.electronDensity;        % electron density (SI units)
    N = boltzmann.workCond.gasDensity;              % gas density (SI units)
    EoN = boltzmann.workCond.reducedFieldSI;        % reduced electric field (SI units)
    WoN = boltzmann.workCond.reducedExcFreqSI;      % reduced angular exitation frequency (SI units)
  end
    
  % renormalize EEDF
  eedf = eedf/sum(eedf.*sqrt(energyCell')*energyStep);
  
  % evaluate time idependent elements of the Boltzmann equation (persistent for performance reasons)
  persistent matrix;
  persistent fieldMatrix;
  persistent nonConservativeMatrix;
  persistent CIEffIntegrand;
  persistent totalCrossSectionAux;
  persistent cellTotalCrossSection;
  persistent growthMatrixDiagElements;
  persistent growthMatrixSupElements;
  persistent growthMatrixInfElements;
  persistent g_extraFieldSpatialGrowth;
  persistent eeMatrixAuxA;
  persistent eeMatrixAuxB;
  if isempty(matrix)
    % evaluate basic boltzmann matrix (without field, ionization, attachment or e-e collisions operators)
    matrix = boltzmann.elasticMatrix + boltzmann.CARMatrix + boltzmann.inelasticMatrix;
    % save local copy of the regular field operator (without growth models)
    fieldMatrix = boltzmann.fieldMatrix;
    
    % include ionization operator (either conservative or not)
    nonConservativeMatrix = zeros(cellNumber);
    if boltzmann.includeNonConservativeIonization
      matrix = matrix + boltzmann.ionizationMatrix;
      nonConservativeMatrix = nonConservativeMatrix + boltzmann.ionizationMatrix;
    else
      matrix = matrix + boltzmann.ionizationConservativeMatrix;
    end
    
    % include attachment operator (either conservative or not)
    if boltzmann.includeNonConservativeAttachment
      matrix = matrix + boltzmann.attachmentMatrix;
      nonConservativeMatrix = nonConservativeMatrix + boltzmann.attachmentMatrix;
    else
      matrix = matrix + boltzmann.attachmentConservativeMatrix;
    end
    
    % evaluate time idependent elements of the growth model operators (in case they are activated)
    if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
      % evaluate integrand to calculate the effective ionization rate
      CIEffIntegrand = sum(nonConservativeMatrix)';
      % local copies of total momentum transfer cross sections (at nodes and cells)
      totalCrossSectionAux = boltzmann.totalCrossSection;
      cellTotalCrossSection = 0.5*(totalCrossSectionAux(1:end-1) + totalCrossSectionAux(2:end));
      % choose growth model for the electron density (either spatial or temporal)
      switch boltzmann.eDensGrowthModel
        case 'temporal'
          % evaluate matrix elements of the temporal growth operator (e-density time variation term, dn/dt)
          growthMatrixDiagElements = -sqrt(me/(2*e))*sqrt(energyCell);
        case 'spatial'
          % evaluate the diffusion and mobility components of the spatial growth operator
          growthMatrixDiagElements = ...    % diffusion component
            energyCell./(3*cellTotalCrossSection);
          growthMatrixSupElements = ...     % mobility component (diag sup)
            EoN/(6*energyStep)*[0 energyCell(1:cellNumber-1)./(cellTotalCrossSection(1:cellNumber-1))];
          growthMatrixInfElements = ...     % mobility component (diag inf)
            -EoN/(6*energyStep)*[energyCell(2:cellNumber)./(cellTotalCrossSection(2:cellNumber)) 0];
          % evaluate components of the extra electric field operator of the spatial growth model
          g_extraFieldSpatialGrowth = EoN*energyNode/energyStep./(6*totalCrossSectionAux);
          g_extraFieldSpatialGrowth(1) = 0;
          g_extraFieldSpatialGrowth(end) = 0;
      end
    end
    
    % evaluate time idependent elements of the e-e collision operator (in case that e-e collisions are activated)
    if boltzmann.includeEECollisions
      % evaluating auxiliary matrix used in the evaluation of the ee upflux vector without multiplicative constant
      eeMatrixAuxA = zeros(cellNumber);
      auxEnergyArray = -(energyStep/2)*sqrt(energyCell)+(2/3)*energyCell.^(3/2);
      % because of power conservation, terms on the last row (cellNumber,:) and on the first column (:,1) are zero
      for k=1:cellNumber-1
        eeMatrixAuxA(k,2:k) = auxEnergyArray(2:k);
        eeMatrixAuxA(k,k+1:cellNumber) = (2/3)*energyNode(k+1)^(3/2);
      end
      % detailed balance condition
      for k=1:cellNumber-1
        eeMatrixAuxA(k,2:cellNumber) = sqrt(eeMatrixAuxA(k,2:cellNumber).*eeMatrixAuxA(1:cellNumber-1,k+1)');
      end
      % evaluating auxiliary matrix used in the evaluation of the ee downflux vector without multiplicative constant
      eeMatrixAuxB = transpose(eeMatrixAuxA);
    end
    
  end
  
  % flush persistent memory for a new integration of the Boltzmann equation
  if clearPersistentVars
    vars = whos;
    vars = vars([vars.persistent]);
    varName = {vars.name};
    clear(varName{:});
    dfdt = [];
    return
  end
  
  % evaluate time derivative (time dependent) of each discrete component of the eedf (except e-e collisions operator)
  if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
    % calculation of the current effective ionization rate (needed for both growth models, spatial and temporal)
    CIEff = sqrt(2*e/me)*energyStep*dot(eedf,CIEffIntegrand);
    % copy CIEff into boltzmann properties (needed for power balance and swarm parameters calculation)
    boltzmann.CIEff = CIEff;
    % choose growth model for the electron density (either spatial or temporal)
    switch boltzmann.eDensGrowthModel
      case 'temporal'
        % writing the modified total momentum transfer cross section (total momentum transfer cross section plus the
        % ionization rate divided by the electron velocity) needed for the reevaluation of the field operator
        totalCrossSectionMod = totalCrossSectionAux + CIEff*sqrt(me/2*e./energyNode);
        % writing the electric field operator matrix of the temporal growth model (with totalCrossSectionMod)
        g_fieldTemporalGrowthAux = EoN^2*energyNode./( 3*energyStep^2*totalCrossSectionMod.* ...
          (1+me*WoN^2./(2*e*energyNode.*totalCrossSectionMod.*totalCrossSectionMod)) );
        g_fieldTemporalGrowthAux(1) = 0;
        g_fieldTemporalGrowthAux(end) = 0;
        % copy field terms into boltzmann properties (needed for power balance and swarm parameters calculation)
        boltzmann.g_fieldTemporalGrowth = g_fieldTemporalGrowthAux*energyStep^2;
        % evaluate time derivative of each discrete component of the eedf (case with temporal growth models)
        dfdt = N*sqrt(2*e/me)*( ( ... % multiplicative constant to obtain proper units of time
          matrix*eedf + ...           % full basic boltzmann matrix
          ( -g_fieldTemporalGrowthAux(1:cellNumber)-g_fieldTemporalGrowthAux(2:cellNumber+1) + ... % time dependent diagonal component
          CIEff*growthMatrixDiagElements)'.*eedf + ...
          [ 0 g_fieldTemporalGrowthAux(2:cellNumber) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...  % time dependent inf. diagonal component
          [ g_fieldTemporalGrowthAux(2:cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...      % time dependent sup. diagonal component
          )./sqrt(energyCell'));    % divide by the square root of the energy to obtain derivative of the EEDF
      case 'spatial'
        % evaluation of gas density times diffusion coefficient and mobility coefficient times the electric field
        ND = sqrt(2*e/me)*energyStep*dot(growthMatrixDiagElements,eedf);
        muE = -sqrt(2*e/me)*energyStep*dot(growthMatrixSupElements+growthMatrixInfElements,eedf);
        % calculation of the effective reduced first Townsend coefficient
        if muE^2-4*CIEff*ND < 0
          alphaRedEff = CIEff/muE;
        else
          alphaRedEff = (muE - sqrt(muE^2-4*CIEff*ND))/(2*ND);
        end
        % copy alphaRedEff into boltzmann properties (needed for power balance and swarm parameters calculation)
        boltzmann.alphaRedEff = alphaRedEff;
        % copy field terms into boltzmann properties (needed for power balance and swarm parameters calculation)
        boltzmann.g_fieldSpatialGrowth = alphaRedEff*energyStep*g_extraFieldSpatialGrowth;
        % evaluate time derivative of each discrete component of the eedf (case with spatial growth models)
        dfdt = N*sqrt(2*e/me)*( ( ...   % multiplicative constant to obtain proper units of time
          matrix*eedf+...               % full basic boltzmann matrix (without field nor growth operators)
          (fieldMatrix(1:cellNumber+1:cellNumber*cellNumber) - ...              % time dependent diagonal component
          alphaRedEff*(g_extraFieldSpatialGrowth(1:cellNumber)-g_extraFieldSpatialGrowth(2:cellNumber+1)) + ...
          alphaRedEff^2*growthMatrixDiagElements)'.*eedf + ...
          [ 0 (fieldMatrix(2:cellNumber+1:cellNumber*cellNumber) - ...          % time dependent inf. diagonal component
          alphaRedEff*g_extraFieldSpatialGrowth(2:cellNumber) + ...
          alphaRedEff*growthMatrixInfElements(1:cellNumber-1)) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...
          [ fieldMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + ...  % time dependent sup. diagonal component
          alphaRedEff*g_extraFieldSpatialGrowth(2:cellNumber) + ...
          alphaRedEff*growthMatrixSupElements(2:cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...
          )./sqrt(energyCell'));        % divide by the square root of the energy to obtain derivative of the EEDF
    end
  else
    % evaluate time derivative of each discrete component of the eedf (case without growth models)
    dfdt = N*sqrt(2*e/me)*( ( ... % multiplicative constant to obtain proper units of time
      matrix*eedf + ...           % full basic boltzmann matrix (without field operator contribution)
      (fieldMatrix(1:cellNumber+1:cellNumber*cellNumber))'.*eedf + ...
      [ 0 fieldMatrix(2:cellNumber+1:cellNumber*cellNumber) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...
      [ fieldMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...
      )./sqrt(energyCell'));      % divide by the square root of the energy to obtain derivative of the EEDF
  end
  
  % evaluate e-e contribution to the time derivative (time dependent) of each discrete component of the eedf
  if boltzmann.includeEECollisions
    % electron temperature in eV Te = (2/3)*meanEnergy
    Te = (2/3)*energyStep*dot(energyCell.^(3/2),eedf);
    % Coulomb logarithm
    logC = log(12*pi*(e0*Te/e)^(3/2)/sqrt(ne));
    % multiplicative constant
    eeConstant = (ne/N)*(e^2/(8*pi*e0^2))*logC;
    % calculation of electron-electron collisions vectors of upflux (A) and downflux (B)
    A = (eeConstant/energyStep)*(eeMatrixAuxA*eedf);
    B = (eeConstant/energyStep)*(eeMatrixAuxB*eedf);
    % copy A and B into boltzmann properties (needed for power balance calculations)
    boltzmann.Afinal = A;
    boltzmann.Bfinal = B;
    % add contribution to time derivative of each discrete component of the eedf due to e-e collisions
    dfdt = dfdt + N*sqrt(2*e/me)*( ( ... % multiplicative constant to obtain proper units of time
      (-A(1:cellNumber)-B(1:cellNumber)).*eedf + ...                       % time dependent diagonal component
      [ 0; A(1:cellNumber-1) ].*[ 0; eedf(1:cellNumber-1) ] + ...          % time dependent inf. diagonal component
      [ B(2:cellNumber); 0 ].*[ eedf(2:cellNumber); 0 ] ...                % time dependent sup. diagonal component
      )./sqrt(energyCell'));    % divide by the square root of the energy to obtain derivative of the EEDF
  end
  
end

function [value,isTerminal,direction] = steadyStateEventFcn(t,eedf,boltzmann,~)
  % evaluate time derivative of each discrete component of the eedf

  persistent eedfOld;
  persistent maxEedfDifference;
  if isempty(eedfOld)
    eedfOld = eedf;
    maxEedfDifference = boltzmann.maxEedfRelError;
    value = 1;
  elseif max(abs(eedf-eedfOld)./eedfOld) < maxEedfDifference
    value = 0;
  else
    eedfOld = eedf;
    value = 1;
  end

  isTerminal = 1;   % stop the integration
  direction = 0;    % stop in any direction
  
end
