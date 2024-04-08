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

classdef Output < handle
  
  properties
    folder = '';    % main output folder
    subFolder = ''; % sub folder for output of different jobs
    
    isSimulationHF = false;             % boolean to know if the electron kinetics (Boltzmann only) is HF
    
    logIsToBeSaved = false;             % boolean to know if the log (as written by the CLI) must be saved
    inputsAreToBeSaved = false;         % boolean to know if the input files must be saved
    eedfIsToBeSaved = false;            % boolean to know if the eedf must be saved
    swarmParamsIsToBeSaved = false;     % boolean to know if the swarm parameters info must be saved
    rateCoeffsIsToBeSaved = false;      % boolean to know if the rate coefficients info must be saved
    powerBalanceIsToBeSaved = false;    % boolean to know if the power balance info must be saved
    lookUpTableIsToBeSaved = false;     % boolean to know if look-up tables with results must be saved
    
  end
  
  methods (Access = public)
    
    function output = Output(setup)
      
      % set output folder (if not specified in the setup, a generic folder with a timestamp is created)
      if isfield(setup.info.output, 'folder')
        output.folder = ['Output' filesep setup.info.output.folder];
      else
        output.folder = ['Output' filesep 'Simulation ' datestr(datetime, 'dd mmm yyyy HHMMSS')];
      end
      % create output folder in case it doesn't exist
      if 7 ~= exist(output.folder, 'file')
        mkdir(output.folder);
      end
      
      % set initial output subfolder (in case multiple jobs are to be run)
      outputSubFolder = '';
      if setup.numberOfJobs > 1
        for i = setup.numberOfBatches:-1:1
          outputSubFolder = sprintf('%s%s%s_%g', outputSubFolder, filesep, setup.batches(i).property, ...
            setup.batches(i).value(1));
        end
      end
      % save output sub folder info (folder inside the output.folder folder)
      output.subFolder = outputSubFolder;
      
      % save what information must be saved
      dataFiles = setup.info.output.dataFiles;
      if ischar(dataFiles)
        dataFiles = {dataFiles};
      end
      for dataFile = dataFiles
        switch dataFile{1}
          case 'log'
            output.logIsToBeSaved = true;
            output.initializeLogFile(setup.cli.logStr);
          case 'inputs'
            output.inputsAreToBeSaved = true;
            output.saveInputFiles(setup);
          case 'eedf'
            output.eedfIsToBeSaved = true;
          case 'swarmParameters'
            output.swarmParamsIsToBeSaved = true;
          case 'rateCoefficients'
            output.rateCoeffsIsToBeSaved = true;
          case 'powerBalance'
            output.powerBalanceIsToBeSaved = true;
          case 'lookUpTable'
            output.lookUpTableIsToBeSaved = true;
        end
      end
      
      % save the setup information for reference (always saved)
      output.saveSetupInfo(setup.unparsedInfo);

      % add listener to status messages of the setup object
      addlistener(setup, 'genericStatusMessage', @output.genericStatusMessage);
      % add listener of the working conditions object
      addlistener(setup.workCond, 'genericStatusMessage', @output.genericStatusMessage);
      
      % add listener to status messages of the electron kinetics object
      addlistener(setup.electronKinetics, 'genericStatusMessage', @output.genericStatusMessage);
      % add listener to output results when a new solution for the EEDF is found
      addlistener(setup.electronKinetics, 'obtainedNewEedf', @output.electronKineticsSolution);
      % save the information if the electron kinetics is HF
      if setup.workCond.reducedExcFreqSI>0
        output.isSimulationHF = true;
      end
      
    end
    
  end
  
  methods (Access = private)
    
    function saveSetupInfo(output, setupCellArray)
    % saveSetupInfo saves the setup of the current simulation
    
      fileName = [output.folder filesep 'setup.txt'];
      fileID = fopen(fileName, 'wt');
      
      for cell = setupCellArray
        fprintf(fileID, '%s\n', cell{1});
      end
      
      fclose(fileID);
      
    end

    function saveInputFiles(output, setup)
    % saveInputFiles saves all the input files found in the setup of the simulation inside an Input folder in the Output
    % folder
      
      % find setup file
      files = {['Input' filesep setup.fileName]};

      % find electron kinetics input files
      if setup.enableElectronKinetics
        % find cross-section files (regular)
        for file = setup.info.electronKinetics.LXCatFiles
          files{end+1} = ['Input' filesep file{1}];
        end
        % find cross-section files (extra)
        if isfield('LXCatFilesExtra', setup.info.electronKinetics)
          for file = setup.info.electronKinetics.LXCatFilesExtra
            files{end+1} = ['Input' filesep file{1}];
          end
        end
        % find gas property files
        for field = fieldnames(setup.info.electronKinetics.gasProperties)'
          entries = setup.info.electronKinetics.gasProperties.(field{1});
          if ischar(entries)
            entries = {entries};
          end
          for entry = entries
            file = ['Input' filesep entry{1}];
            if isfile(file)
              files{end+1} = file;
            end
          end
        end
        % find state property files
        for field = fieldnames(setup.info.electronKinetics.stateProperties)'
          entries = setup.info.electronKinetics.stateProperties.(field{1});
          if ischar(entries)
            entries = {entries};
          end
          for entry = entries
            file = ['Input' filesep entry{1}];
            if isfile(file)
              files{end+1} = file;
            end
          end
        end
      end

      % create Input folder inside the current output folder
      inputOutputFolder = [output.folder filesep 'Input'];
      if ~isfolder(inputOutputFolder)
        mkdir(inputOutputFolder);
      end
      % copy input files to output folder
      for file = files
        finalFile = [output.folder filesep file{1}];
        [finalFolder, ~, ~] = fileparts(finalFile);
        if ~isfolder(finalFolder)
          mkdir(finalFolder);
        end
        copyfile(file{1}, finalFile);
      end

    end
    
    function initializeLogFile(output, logCellArray)
    % initializeLogFile initialized the output file containing the log of the simulation and writes previous messages of
    % the log produced before the creation of the output object
    
      fileName = [output.folder filesep 'log.txt'];
      fileID = fopen(fileName, 'wt');
      
      for cell = logCellArray
        fprintf(fileID, '%s\n', cell{1});
      end
      
      fclose(fileID);
      
    end
    
    function genericStatusMessage(output, ~, statusEventData)
      
      if output.logIsToBeSaved
        fileName = [output.folder filesep 'log.txt'];
        fileID = fopen(fileName, 'at');
        fprintf(fileID, statusEventData.message);
        fclose(fileID);
      end

    end

    function electronKineticsSolution(output, electronKinetics, ~)
    
      % create subfolder name in case of time-dependent boltzmann calculations
      if isa(electronKinetics, 'Boltzmann') && electronKinetics.isTimeDependent
        output.subFolder = sprintf('%stime_%e', filesep, electronKinetics.workCond.currentTime);
      end
      % create subfolder in case it is needed (when performing runs of simmulations or in time-dependent Boltzmann)
      if ~isempty(output.subFolder) && (output.eedfIsToBeSaved || output.powerBalanceIsToBeSaved || ...
          output.swarmParamsIsToBeSaved || output.rateCoeffsIsToBeSaved )
        if 7 ~= exist([output.folder output.subFolder], 'file')
          mkdir([output.folder output.subFolder]);
        end
      end
      
      % save selected results of the electron kinetics
      if output.eedfIsToBeSaved
        if isa(electronKinetics, 'Boltzmann')
          output.saveEedf(electronKinetics.eedf, electronKinetics.firstAnisotropy, electronKinetics.energyGrid.cell);
        else
          output.saveEedf(electronKinetics.eedf, [], electronKinetics.energyGrid.cell);
        end
      end
      if output.swarmParamsIsToBeSaved
        output.saveSwarm(electronKinetics.swarmParam, electronKinetics.workCond.reducedField);
      end
      if output.rateCoeffsIsToBeSaved
        output.saveRateCoefficients(electronKinetics.rateCoeffAll, electronKinetics.rateCoeffExtra, []);
      end
      if output.powerBalanceIsToBeSaved
        output.savePower(electronKinetics.power);
      end
      if output.lookUpTableIsToBeSaved
        output.saveLookUpTable(electronKinetics);
      end
      
    end
    
    function saveEedf(output, eedf, firstAnisotropy, energy)
    % saveEedf saves the eedf information of the current simulation
      
      % create file name
      fileName = [output.folder output.subFolder filesep 'eedf.txt'];
      
      % open file
      fileID = fopen(fileName, 'wt');
      
      % save information into the file
      if isempty(firstAnisotropy)
        fprintf(fileID, 'Energy(eV)           EEDF(eV^-(3/2))\n');
        values(2:2:2*length(eedf)) = eedf;
        values(1:2:2*length(eedf)) = energy;
        fprintf(fileID, '%#.14e %#.14e \n', values);
      else
        fprintf(fileID, 'Energy(eV)           EEDF(eV^-(3/2))      Anisotropy(eV^-(3/2))\n');
        values(3:3:3*length(eedf)) = firstAnisotropy;
        values(2:3:3*length(eedf)) = eedf;
        values(1:3:3*length(eedf)) = energy;
        fprintf(fileID, '%#.14e %#.14e %#.14e \n', values);
      end
      
      % close file
      fclose(fileID);
      
    end
    
    function saveSwarm(output, swarmParam, reducedField)
    % saveSwarm saves the swarm parameters information of the current simulation
      
      % create file name
      fileName = [output.folder output.subFolder filesep 'swarmParameters.txt'];
      
      % open file
      fileID = fopen(fileName, 'wt');
      
      % save information into the file
      fprintf(fileID, '               Reduced electric field = %#.14e (Td)\n', reducedField);
      fprintf(fileID, '        Reduced diffusion coefficient = %#.14e ((ms)^-1)\n', swarmParam.redDiffCoeff);
      fprintf(fileID, '                     Reduced mobility = %#.14e ((msV)^-1)\n', swarmParam.redMobility);
      if output.isSimulationHF 
        fprintf(fileID, '                  Reduced mobility HF = %#.14e%+#.14ei ((msV)^-1)\n', ...
          real(swarmParam.redMobilityHF), imag(swarmParam.redMobilityHF));
      else
        fprintf(fileID, '                       Drift velocity = %#.14e (ms^-1)\n', swarmParam.driftVelocity);
        fprintf(fileID, '         Reduced Townsend coefficient = %#.14e (m^2)\n', swarmParam.redTownsendCoeff);
        fprintf(fileID, '       Reduced attachment coefficient = %#.14e (m^2)\n', swarmParam.redAttCoeff);
      end
      fprintf(fileID, ' Reduced energy diffusion coefficient = %#.14e (eV(ms)^-1)\n', swarmParam.redDiffCoeffEnergy);
      fprintf(fileID, '              Reduced energy mobility = %#.14e (eV(msV)^-1)\n', swarmParam.redMobilityEnergy);
      fprintf(fileID, '                          Mean energy = %#.14e (eV)\n', swarmParam.meanEnergy);
      fprintf(fileID, '                Characteristic energy = %#.14e (eV)\n', swarmParam.characEnergy);
      fprintf(fileID, '                 Electron temperature = %#.14e (eV)\n', swarmParam.Te);
      
      % close file
      fclose(fileID);
      
    end
    
    function saveRateCoefficients(output, eKineticsRateCoeffs, eKineticsRateCoeffsExtra, reactionsInfo)
    % saveRateCoefficients saves the rate coefficients obtained in the current simulation
      
      % create file name
      fileName = [output.folder output.subFolder filesep 'rateCoefficients.txt'];
      
      % open file
      fileID = fopen(fileName, 'wt');
      
      % save information into the file
      if ~isempty(eKineticsRateCoeffs)
        fprintf(fileID, '%s\n*    e-Kinetics Rate Coefficients    *\n%s\n\n', repmat('*', 1,38), repmat('*', 1,38));
        fprintf(fileID, 'ID   Ine.R.Coeff.(m^3s^-1) Sup.R.Coeff.(m^3s^-1) Threshold(eV)         Description\n');
        for rateCoeff = eKineticsRateCoeffs
          if length(rateCoeff.value) == 1
            fprintf(fileID, '%4d %20.14e  (N/A)                 %20.14e  %s\n', rateCoeff.collID, rateCoeff.value, ...
              rateCoeff.energy, rateCoeff.collDescription);
          else
            fprintf(fileID, '%4d %20.14e  %20.14e  %20.14e  %s\n', rateCoeff.collID, rateCoeff.value(1), ...
              rateCoeff.value(2), rateCoeff.energy, rateCoeff.collDescription);
          end
        end
      end
      if ~isempty(eKineticsRateCoeffsExtra)
        fprintf(fileID, '\n%s\n* e-Kinetics Extra Rate Coefficients *\n%s\n\n', repmat('*', 1,38), repmat('*', 1,38));
        fprintf(fileID, 'ID   Ine.R.Coeff.(m^3s^-1) Sup.R.Coeff.(m^3s^-1) Threshold(eV)         Description\n');
        for rateCoeff = eKineticsRateCoeffsExtra
          if length(rateCoeff.value) == 1
            fprintf(fileID, '%4d %20.14e  (N/A)                 %20.14e  %s\n', rateCoeff.collID, rateCoeff.value, ...
              rateCoeff.energy, rateCoeff.collDescription);
          else
            fprintf(fileID, '%4d %20.14e  %20.14e  %20.14e  %s\n', rateCoeff.collID, rateCoeff.value(1), ...
              rateCoeff.value(2), rateCoeff.energy, rateCoeff.collDescription);
          end
        end
      end
      if ~isempty(reactionsInfo)
        fprintf(fileID, '\n%s\n*     Chemistry Rate Coefficients    *\n%s\n\n', repmat('*', 1,38), repmat('*', 1,38));
        fprintf(fileID, ['ID   Dir.R.Coeff.(S.I.)    Inv.R.Coeff.(S.I.)    Enthalpy(eV)          ' ...
          'Net.Reac.Rate(m^-3s^-1) Description\n']);
        for reaction = reactionsInfo
          if length(reaction.rateCoeff) == 1
            fprintf(fileID, '%4d %20.14e  (N/A)                 %+20.14e %+20.14e   %s\n', reaction.reactID, ...
              reaction.rateCoeff, reaction.energy, reaction.netRate, reaction.description);
          else
            fprintf(fileID, '%4d %20.14e  %20.14e  %+20.14e %+20.14e   %s\n', reaction.reactID, ...
              reaction.rateCoeff(1), reaction.rateCoeff(2), reaction.energy, reaction.netRate, reaction.description);
          end
        end
      end

      % close file
      fclose(fileID);
      
    end
    
    function savePower(output, power)
    % savePower saves the power balance information of the current simulation
      
      % create file name
      fileName = [output.folder output.subFolder filesep 'powerBalance.txt'];
      
      % open file
      fileID = fopen(fileName, 'wt');
      
      % save information into the file
      fprintf(fileID, '                               Field = %#+.14e (eVm^3s^-1)\n', power.field);
      fprintf(fileID, '           Elastic collisions (gain) = %#+.14e (eVm^3s^-1)\n', power.elasticGain);
      fprintf(fileID, '           Elastic collisions (loss) = %#+.14e (eVm^3s^-1)\n', power.elasticLoss);
      fprintf(fileID, '                          CAR (gain) = %#+.14e (eVm^3s^-1)\n', power.carGain);
      fprintf(fileID, '                          CAR (loss) = %#+.14e (eVm^3s^-1)\n', power.carLoss);
      fprintf(fileID, '     Excitation inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.excitationIne);
      fprintf(fileID, '  Excitation superelastic collisions = %#+.14e (eVm^3s^-1)\n', power.excitationSup);
      fprintf(fileID, '    Vibrational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.vibrationalIne);
      fprintf(fileID, ' Vibrational superelastic collisions = %#+.14e (eVm^3s^-1)\n', power.vibrationalSup);
      fprintf(fileID, '     Rotational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.rotationalIne);
      fprintf(fileID, '  Rotational superelastic collisions = %#+.14e (eVm^3s^-1)\n', power.rotationalSup);
      fprintf(fileID, '               Ionization collisions = %#+.14e (eVm^3s^-1)\n', power.ionizationIne);
      fprintf(fileID, '               Attachment collisions = %#+.14e (eVm^3s^-1)\n', power.attachmentIne);
      fprintf(fileID, '             Electron density growth = %#+.14e (eVm^3s^-1) +\n', power.eDensGrowth);
      fprintf(fileID, ' %s\n', repmat('-', 1, 73));
      fprintf(fileID, '                       Power Balance = %#+.14e (eVm^3s^-1)\n', power.balance);
      fprintf(fileID, '              Relative Power Balance = % #.14e%%\n\n', power.relativeBalance*100);
      fprintf(fileID, '           Elastic collisions (gain) = %#+.14e (eVm^3s^-1)\n', power.elasticGain);
      fprintf(fileID, '           Elastic collisions (loss) = %#+.14e (eVm^3s^-1) +\n', power.elasticLoss);
      fprintf(fileID, ' %s\n', repmat('-', 1, 73));
      fprintf(fileID, '            Elastic collisions (net) = %#+.14e (eVm^3s^-1)\n\n', power.elasticNet);
      fprintf(fileID, '                          CAR (gain) = %#+.14e (eVm^3s^-1)\n', power.carGain);
      fprintf(fileID, '                          CAR (gain) = %#+.14e (eVm^3s^-1) +\n', power.carLoss);
      fprintf(fileID, ' %s\n', repmat('-', 1, 73));
      fprintf(fileID, '                           CAR (net) = %#+.14e (eVm^3s^-1)\n\n', power.carNet);
      fprintf(fileID, '     Excitation inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.excitationIne);
      fprintf(fileID, '  Excitation superelastic collisions = %#+.14e (eVm^3s^-1) +\n', power.excitationSup);
      fprintf(fileID, ' %s\n', repmat('-', 1, 73));
      fprintf(fileID, '         Excitation collisions (net) = %#+.14e (eVm^3s^-1)\n\n', power.excitationNet);
      fprintf(fileID, '    Vibrational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.vibrationalIne);
      fprintf(fileID, ' Vibrational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', power.vibrationalSup);
      fprintf(fileID, ' %s\n', repmat('-', 1, 73));
      fprintf(fileID, '        Vibrational collisions (net) = %#+.14e (eVm^3s^-1)\n\n', power.vibrationalNet);
      fprintf(fileID, '     Rotational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.rotationalIne);
      fprintf(fileID, '  Rotational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', power.rotationalSup);
      fprintf(fileID, ' %s\n', repmat('-', 1, 73));
      fprintf(fileID, '         Rotational collisions (net) = %#+.14e (eVm^3s^-1)\n', power.rotationalNet);
      
      % power balance by gases
      gases = fields(power.gases);
      powerByGas = power.gases;
      for i = 1:length(gases)
        gas = gases{i};
        fprintf(fileID, '\n%s\n\n', [repmat('*', 1, 37) ' ' gas ' ' repmat('*', 1, 39-length(gas))]);
        fprintf(fileID, '     Excitation inelastic collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).excitationIne);
        fprintf(fileID, '  Excitation superelastic collisions = %#+.14e (eVm^3s^-1) +\n', powerByGas.(gas).excitationSup);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '         Excitation collisions (net) = %#+.14e (eVm^3s^-1)\n\n', powerByGas.(gas).excitationNet);
        fprintf(fileID, '    Vibrational inelastic collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).vibrationalIne);
        fprintf(fileID, ' Vibrational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', powerByGas.(gas).vibrationalSup);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '        Vibrational collisions (net) = %#+.14e (eVm^3s^-1)\n\n', powerByGas.(gas).vibrationalNet);
        fprintf(fileID, '     Rotational inelastic collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).rotationalIne);
        fprintf(fileID, '  Rotational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', powerByGas.(gas).rotationalSup);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '         Rotational collisions (net) = %#+.14e (eVm^3s^-1)\n\n', powerByGas.(gas).rotationalNet);
        fprintf(fileID, '               Ionization collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).ionizationIne);
        fprintf(fileID, '               Attachment collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).attachmentIne);
      end
      % close file
      fclose(fileID);
      
    end
    
    function saveLookUpTable(output, electronKinetics)
      
      % name of the files containing the different lookup tables
      persistent fileName1;
      persistent fileName2;
      persistent fileName3;
      persistent fileName4;
      persistent fileName5;
      
      % local copies of different variables (for performance reasons)
      workCond = electronKinetics.workCond;
      power = electronKinetics.power;
      swarmParams = electronKinetics.swarmParam;
      rateCoeffAll = electronKinetics.rateCoeffAll;
      rateCoeffExtra = electronKinetics.rateCoeffExtra;
      eedf = electronKinetics.eedf;
      
      % initialize the files in case it is needed
      if isempty(fileName1)
        % create file names
        fileName1 = [output.folder filesep 'lookUpTableSwarm.txt'];
        fileName2 = [output.folder filesep 'lookUpTablePower.txt'];
        fileName3 = [output.folder filesep 'lookUpTableRateCoeff.txt'];
        % open files
        fileID1 = fopen(fileName1, 'wt');
        fileID2 = fopen(fileName2, 'wt');
        fileID3 = fopen(fileName3, 'wt');
        % write file headers
        fprintf(fileID3, [repmat('#', 1, 80) '\n# %-76s #\n'], 'ID   Description');
        strFile3 = '';
        for i = 1:length(rateCoeffAll)
          fprintf(fileID3, '# %-4d %-71s #\n', rateCoeffAll(i).collID, rateCoeffAll(i).collDescription);
          strAux = sprintf('R%d_ine(m^3s^-1)', rateCoeffAll(i).collID);
          strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          if 2 == length(rateCoeffAll(i).value)
            strAux = sprintf('R%d_sup(m^3s^-1)', rateCoeffAll(i).collID);
            strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          end
        end
        fprintf(fileID3, '#%s#\n# %-76s #\n#%s#\n# %-76s #\n', repmat(' ', 1, 78), ...
          '*** Extra rate coefficients ***', repmat(' ', 1, 78), 'ID   Description');
        for i = 1:length(rateCoeffExtra)
          fprintf(fileID3, '# %-4d %-71s #\n', rateCoeffExtra(i).collID, rateCoeffExtra(i).collDescription);
          strAux = sprintf('R%d_ine(m^3s^-1)', rateCoeffExtra(i).collID);
          strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          if 2 == length(rateCoeffExtra(i).value)
            strAux = sprintf('R%d_sup(m^3s^-1)', rateCoeffExtra(i).collID);
            strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          end
        end
        fprintf(fileID3, [repmat('#', 1, 80) '\n\n']);
        if isa(electronKinetics, 'Boltzmann')
          if electronKinetics.isTimeDependent
            fprintf(fileID1, '%-21s ', 'Time(s)');
            fprintf(fileID2, '%-21s ', 'Time(s)');
            fprintf(fileID3, '%-21s ', 'Time(s)');
            % create lookup table for the eedf
            fileName4 = [output.folder filesep 'lookUpTableEedf.txt'];
            fileID4 = fopen(fileName4, 'wt');
            % add first line with energies to eedf lookup table (eedfs will be saved as rows)
            fprintf(fileID4, '%-21.14e ', [0 electronKinetics.energyGrid.cell]);
            fprintf(fileID4, '\n');
            fclose(fileID4);
            % create lookup table for the electron density (if needed)
            if electronKinetics.eDensIsTimeDependent
              fileName5 = [output.folder filesep 'lookUpTableElectronDensity.txt'];
              fileID5 = fopen(fileName5, 'wt');
              fprintf(fileID5, '%-21s %-21s\n', 'time(s)', 'ne(m^-3)\n');
              fclose(fileID5);
            end
          end
          if output.isSimulationHF
            fprintf(fileID1, [repmat('%-21s ', 1, 10) '\n'], 'RedField(Td)', 'RedDiff((ms)^-1)', 'RedMob((msV)^-1)', ...
              'R[RedMobHF]((msV)^-1)', 'I[RedMobHF]((msV)^-1)', 'RedDiffE(eV(ms)^-1)', 'RedMobE(eV(msV)^-1)', ...
              'MeanE(eV)', 'CharE(eV)', 'EleTemp(eV)');
          else
            fprintf(fileID1, [repmat('%-21s ', 1, 11) '\n'], 'RedField(Td)', 'RedDiff((ms)^-1)', 'RedMob((msV)^-1)', ...
              'DriftVelocity(ms^-1)', 'RedTow(m^2)', 'RedAtt(m^2)', 'RedDiffE(eV(ms)^-1)', 'RedMobE(eV(msV)^-1)', ...
              'MeanE(eV)', 'CharE(eV)', 'EleTemp(eV)');
          end
          fprintf(fileID2, '%-21s ', 'RedField(Td)');
          fprintf(fileID3, '%-21s ', 'RedField(Td)');
        else
          if output.isSimulationHF
            fprintf(fileID1, [repmat('%-21s ', 1, 10) '\n'], 'EleTemp(eV)', 'RedField(Td)', 'RedDiff(1/(ms))', ...
              'RedMob(1/(msV))', 'R[RedMobHF](1/(msV))', 'I[RedMobHF](1/(msV))', 'RedDiffE(eV/(ms))', ...
              'RedMobE(eV/(msV))', 'MeanE(eV)', 'CharE(eV)');
          else
            fprintf(fileID1, [repmat('%-21s ', 1, 11) '\n'], 'EleTemp(eV)', 'RedField(Td)', 'RedDiff(1/(ms))', ...
            'RedMob(1/(msV))', 'DriftVelocity(m/s)', 'RedTow(m2)', 'RedAtt(m2)', 'RedDiffE(eV/(ms))', 'RedMobE(eV/(msV))', 'MeanE(eV)', ...
            'CharE(eV)');            
          end
          fprintf(fileID2, '%-21s ', 'EleTemp(eV)');
          fprintf(fileID3, '%-21s ', 'EleTemp(eV)');
        end
        fprintf(fileID2, [repmat('%-21s ', 1, 21) '\n'], 'PowerField(eVm^3s^-1)', ...
          'PwrElaGain(eVm^3s^-1)', 'PwrElaLoss(eVm^3s^-1)', 'PwrElaNet(eVm^3s^-1)', 'PwrCARGain(eVm^3s^-1)', ...
          'PwrCARLoss(eVm^3s^-1)', 'PwrCARNet(eVm^3s^-1)', 'PwrEleGain(eVm^3s^-1)', 'PwrEleLoss(eVm^3s^-1)', ...
          'PwrEleNet(eVm^3s^-1)', 'PwrVibGain(eVm^3s^-1)', 'PwrVibLoss(eVm^3s^-1)', 'PwrVibNet(eVm^3s^-1)', ...
          'PwrRotGain(eVm^3s^-1)', 'PwrRotLoss(eVm^3s^-1)', 'PwrRotNet(eVm^3s^-1)', 'PwrIon(eVm^3s^-1)', ...
          'PwrAtt(eVm^3s^-1)', 'PwrGroth(eVm^3s^-1)', 'PwrBalance(eVm^3s^-1)', 'RelPwrBalance');
        fprintf(fileID3, '%s\n', strFile3);
        % close files
        fclose(fileID1);
        fclose(fileID2);
        fclose(fileID3);
      end
      
      % check if eedf lookup table needs to be saved (and append new line with data)
      if ~isempty(fileName4)
        fileID4 = fopen(fileName4, 'at');
        fprintf(fileID4, '%-21.14e ', workCond.currentTime);
        fprintf(fileID4, '%-21.14e ', eedf);
        fprintf(fileID4, '\n');
        fclose(fileID4);
      end
      % check if electron density data needs to be saved (and append new line with data)
      if ~isempty(fileName5)
        fileID5 = fopen(fileName5, 'at');
        fprintf(fileID5, '%#.14e %#.14e\n',workCond.currentTime, workCond.electronDensity);
        fclose(fileID5);
      end

      % open files
      fileID1 = fopen(fileName1, 'at');
      fileID2 = fopen(fileName2, 'at');
      fileID3 = fopen(fileName3, 'at');
      % append new lines with data
      if isa(electronKinetics, 'Boltzmann')
        if electronKinetics.isTimeDependent
          fprintf(fileID1, '%-+21.14e ', workCond.currentTime);
          fprintf(fileID2, '%-+21.14e ', workCond.currentTime);
          fprintf(fileID3, '%-+21.14e ', workCond.currentTime);
        end
        if output.isSimulationHF
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 10) '\n'], ...
            workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, ...
            real(swarmParams.redMobilityHF), imag(swarmParams.redMobilityHF), swarmParams.redDiffCoeffEnergy, ...
            swarmParams.redMobilityEnergy, swarmParams.meanEnergy, swarmParams.characEnergy, swarmParams.Te);
        else
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 11) '\n'], ...
            workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, swarmParams.driftVelocity, ...
            swarmParams.redTownsendCoeff, swarmParams.redAttCoeff, swarmParams.redDiffCoeffEnergy, ...
            swarmParams.redMobilityEnergy, swarmParams.meanEnergy, swarmParams.characEnergy, swarmParams.Te);
        end
        fprintf(fileID2, '%-+21.14e ', workCond.reducedField);
        fprintf(fileID3, '%-+21.14e ', workCond.reducedField);
      else
        if output.isSimulationHF
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 10) '\n'], ...
            swarmParams.Te, workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, ...
            real(swarmParams.redMobilityHF), imag(swarmParams.redMobilityHF), swarmParams.redDiffCoeffEnergy, ...
            swarmParams.redMobilityEnergy, swarmParams.meanEnergy, swarmParams.characEnergy);
        else
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 11) '\n'], ...
            swarmParams.Te, workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, ...
            swarmParams.driftVelocity, swarmParams.redTownsendCoeff, swarmParams.redAttCoeff, ...
            swarmParams.redDiffCoeffEnergy, swarmParams.redMobilityEnergy, swarmParams.meanEnergy, ...
            swarmParams.characEnergy);
        end
        fprintf(fileID2, '%-+21.14e ', workCond.electronTemperature);
        fprintf(fileID3, '%-+21.14e ', workCond.electronTemperature);
      end
      fprintf(fileID2, [repmat('%-+21.14e ', 1, 20) '%19.14e%%\n'], power.field, ...
        power.elasticGain, power.elasticLoss, power.elasticNet, power.carGain, power.carLoss, power.carNet, ...
        power.excitationSup, power.excitationIne, power.excitationNet, power.vibrationalSup, power.vibrationalIne, ...
        power.vibrationalNet, power.rotationalSup, power.rotationalIne, power.rotationalNet, power.ionizationIne, ...
        power.attachmentIne, power.eDensGrowth, power.balance, power.relativeBalance*100);
      for i = 1:length(rateCoeffAll)
        fprintf(fileID3, '%-21.14e ', rateCoeffAll(i).value(1));
        if 2 == length(rateCoeffAll(i).value)
          fprintf(fileID3, '%-21.14e ', rateCoeffAll(i).value(2));
        end
      end
      for i = 1:length(rateCoeffExtra)
        fprintf(fileID3, '%-21.14e ', rateCoeffExtra(i).value(1));
        if 2 == length(rateCoeffExtra(i).value)
          fprintf(fileID3, '%-21.14e ', rateCoeffExtra(i).value(2));
        end
      end
      fprintf(fileID3, '\n');
      % close files
      fclose(fileID1);
      fclose(fileID2);
      fclose(fileID3);
      
    end
    
  end

end
