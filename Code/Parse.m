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

classdef Parse
  
  properties (Constant = true)
    
    inputFolder = 'Input';
    commentChar = '%';
    wildCardChar = '*';
    stateRegExp = ['\s*(?<quantity>\d+\s*)?(?<gasName>\w+)\((?<ionCharg>[+-])?'...
      '(?:,)?(?<eleLevel>[\w''.\[\]/+*-]+){1}(?:,v=)?(?<vibLevel>(?<=,v=)[\w+*-]+)?'...
      '(?:,J=)?(?<rotLevel>(?<=,J=)[\d+*-]+)?\)\s*'];
    electronRegExp = '\s*(?<quantity>\d+)?((?:[eE]\s*)(?=[\s]|[+][^\(])|(?:[eE]$))';
    
  end
  
  methods (Static)
    
    function [setupStruct, unparsed] = setupFile(fileName)
    % setupFile Reads the configuration file for the simulation and returns
    % a setup structure with all the information.
      
      [structArray, unparsed] = file2structArray([Parse.inputFolder filesep fileName]);
      setupStruct = structArray2struct(structArray);
      
    end
    
    function LXCatEntryArray = LXCatFiles(fileName)
    % LXCatFiles Reads LXCat files, parse their content and returns an
    % structure array 'LXCatEntryArray' with all the information.
      
      % definition of regular expressions to parse LXCat file
      LXCatRegExp1 = 'PROCESS: (?<reactants>.+?)(?<direction>->|<->) (?<products>.+?), (?<type>\w+)';
      LXCatRegExp2 = 'E = (?<threshold>[\d.]+) eV';
      LXCatRegExp3 = '\[(?<reactants>.+?)(?<direction>->|<->)(?<products>.+?), (?<type>\w+)\]';
      % create a cell array with filenames in case only one file is
      % received as input
      if ischar(fileName)
        fileName = {fileName};
      end
      % create an empty struct array of LXCat entries 
      LXCatEntryArray = struct.empty;
      % loop over different LXCat files that have to be read
      for i = 1:length(fileName)
        % open LXCat file
        fileID = fopen([Parse.inputFolder filesep fileName{i}], 'r');
        if fileID<0
          error(' Unable to open LXCat file: %s\n', fileName{i});
        end
        % parse LXCat file
        while ~feof(fileID)
          description = regexp(fgetl(fileID), LXCatRegExp1, 'names', 'once');
          if ~isempty(description)
            parameter = regexp(fgetl(fileID), LXCatRegExp2, 'names', 'once');
            description = regexp(fgetl(fileID),LXCatRegExp3, 'names', 'once');
            while ~strncmp(fgetl(fileID), '-----', 5)
            end
            rawCrossSection = (fscanf(fileID, '%f', [2 inf]));
            % add LXCat entry information into LXCatEntry struct array
            LXCatEntryArray = addLXCatEntry(description, parameter, rawCrossSection, LXCatEntryArray);
          end
        end
        % close LXCat file
        fclose(fileID);
      end

    end

    function gasAndValueArray = gasPropertyFile(fileName)
    % gasPropertyFile Reads a file with gas properties, parse its content
    % and returns an structure array 'gasAndValueArray' with all the
    % information.
    
      regExp = '(?<gasName>\w+)\s*(?<valueStr>.+)';
      fileID = fopen([Parse.inputFolder filesep fileName], 'r');
      if fileID == -1
        error('\t Unable to open file: %s\n', [Parse.inputFolder filesep fileName]);
      end
      
      gasAndValueArray = struct.empty;
      while ~feof(fileID)
        line = removeComments(fgetl(fileID));
        if isempty(line)
          continue
        end
        gasProperty = regexp(line, regExp, 'names', 'once');
        if ~isempty(gasProperty)
          gasAndValueArray(end+1).gasName = gasProperty.gasName;
          gasAndValueArray(end).value = str2num(gasProperty.valueStr);
        end
      end
      
      fclose(fileID);
      
    end
    
    function parsedEntry = gasPropertyEntry(entry)
    % gasPropertyEntry parses an entry (string) of the input file that
    % contains information related to a certain gas property. It returns an
    % structure with the parsed information. 
    
      regExp = '(?<gasName>\w+)\s*=\s*(?<value>.+)';
      parsedEntry = regexp(entry, regExp, 'names');
      if isempty(parsedEntry)
        parsedEntry(end+1).fileName = entry;
        parsedEntry = rmfield(parsedEntry, 'gasName');
        parsedEntry = rmfield(parsedEntry, 'value');
      end
      
    end
    
    function stateAndValueArray = statePropertyFile(fileName)
    % statePropertyFile Reads a file with state properties, parse its
    % content and returns an structure array 'stateAndValueArray' with all
    % the information.
    
      regExp = [Parse.stateRegExp '\s*(?<valueStr>.+)'];
      fileID = fopen([Parse.inputFolder filesep fileName], 'r');
      if fileID == -1
        error('\t Unable to open file: %s\n', [Parse.inputFolder filesep fileName]);
      end
      
      stateAndValueArray = struct.empty;
      while ~feof(fileID)
        line = removeComments(fgetl(fileID));
        if isempty(line)
          continue
        end
        stateProperty = regexp(line, regExp, 'names', 'once');
        if ~isempty(stateProperty)
          stateAndValueArray(end+1).gasName = stateProperty.gasName;
          stateAndValueArray(end).ionCharg = stateProperty.ionCharg;
          stateAndValueArray(end).eleLevel = stateProperty.eleLevel;
          stateAndValueArray(end).vibLevel = stateProperty.vibLevel;
          stateAndValueArray(end).rotLevel = stateProperty.rotLevel;
          stateAndValueArray(end).value = str2num(stateProperty.valueStr);
        end
      end
      
    end
    
    function parsedEntry = statePropertyEntry(entry)
    % statePropertyEntry parses an entry (string) of the input file that
    % contains information related to a certain state property. It returns
    % an structure with the parsed information. 
    
      regExp = [Parse.stateRegExp '\s*=\s*(?<constant>[\d.eE\s()*/+-]+)?(?<functionHandle>@[\w]+)?'...
        '(?<function>\w+)?(?(function)@?)(?<argument>.+)?\s*'];
      parsedEntry = regexp(entry, regExp, 'names');
      if isempty(parsedEntry)
        parsedEntry(end+1).fileName = entry;
        parsedEntry = rmfield(parsedEntry, 'gasName');
        parsedEntry = rmfield(parsedEntry, 'ionCharg');
        parsedEntry = rmfield(parsedEntry, 'eleLevel');
        parsedEntry = rmfield(parsedEntry, 'vibLevel');
        parsedEntry = rmfield(parsedEntry, 'rotLevel');
        parsedEntry = rmfield(parsedEntry, 'function');
        parsedEntry = rmfield(parsedEntry, 'argument');
      elseif isempty(parsedEntry.constant)
        if isempty(parsedEntry.functionHandle)
          if isempty(parsedEntry.argument)
            parsedEntry.argument = cell(0);
          else
            parsedEntry.argument = strsplit(parsedEntry.argument, ',');
            for i=1:length(parsedEntry.argument)
              numericArgument = str2num(parsedEntry.argument{i});
              if ~isnan(numericArgument)
                parsedEntry.argument{i} = numericArgument;
              else
                parsedEntry.argument{i} = parsedEntry.argument{i};
              end
            end
          end
        else
          parsedEntry.function = 'functionHandle';
          parsedEntry.argument = {str2func(parsedEntry.functionHandle)};
        end
      else
        parsedEntry.function = 'constantValue';
        parsedEntry.argument = {str2num(parsedEntry.constant)};
      end
      parsedEntry = rmfield(parsedEntry, 'quantity');
      parsedEntry = rmfield(parsedEntry, 'constant');
      parsedEntry = rmfield(parsedEntry, 'functionHandle');
      
    end
    
  end
  
end

function [structArray, rawLine] = file2structArray(file)
% file2structArray Reads an input file and and creates an array of input
% structs with all the information.
%
% See also structArray2struct
      
  fileID = fopen(file, 'r');
  if fileID == -1
    error('\t Unable to open file: %s\n', file);
  end
  structArray = struct.empty;
  rawLine = cell.empty;
  while ~feof(fileID)
    rawLine{end+1} = removeComments(fgetl(fileID));
    line = strtrim(rawLine{end});
    if isempty(line)
      rawLine = rawLine(1:end-1);
      continue;
    elseif strncmp(line, '-', 1)
      structArray(end).value{end+1} = str2value(line(2:end));
    else
      nameAndValue=strtrim(strsplit(line,':'));
      structArray(end+1).name = nameAndValue{1};
      structArray(end).value = str2value(nameAndValue{2});
      structArray(end).level = regexp(rawLine{end}, '\S', 'once');
    end
  end
  fclose(fileID);
      
end

function structure = structArray2struct(structArray)
% structArray2struct Convert an array of structures into a single structure.
%
% See also file2structArray
    
  i = 1;
  iMax = length(structArray);
  while i<=iMax
    if (i<iMax && structArray(i).level<structArray(i+1).level)
      counter = 1;
      for j = i+2:iMax
        if (structArray(j).level>structArray(i).level)
          counter = counter+1;
        else
          break
        end
      end
      structure.(structArray(i).name) = ...
        structArray2struct(structArray(i+1:i+counter));
      i = i+counter+1;
    else
      structure.(structArray(i).name) = structArray(i).value;
      i = i+1;
    end
  end
      
end

function textLineClean = removeComments(textLine)
% removeComments Return a string without comments.
%
% Comment delimiter is defined in the property commentChar of the Parse class.
% Note: Comment symbol is NOT ignored even if it's inside a string.
  
  idx = regexp(textLine, Parse.commentChar, 'once'); 
  if isempty(idx)
    textLineClean = textLine;
  else
    textLineClean = textLine(1:idx-1);
  end
end

function value = str2value(str)
% str2value Converts a string to a value. Posible values are: numeric,
% logical or string

  [value, is_number_or_logical] = str2num(strtrim(str));
  if ~is_number_or_logical
    value = strtrim(str);
  end
      
end

function LXCatEntryArray = addLXCatEntry(description, parameter, rawCrossSection, LXCatEntryArray)
% addLXCatEntry analyses the information of a particular LXCat entry and adds it to the structure array
% LXCatEntryArray.
      
  LXCatEntryArray(end+1).type = description.type;
  if strcmp(description.direction, '->')
    LXCatEntryArray(end).isReverse = false;
  elseif strcmp(description.direction, '<->')
    LXCatEntryArray(end).isReverse = true;
  end
  LXCatEntryArray(end).target = regexp(description.reactants, Parse.stateRegExp, 'names', 'once');
  if isempty(LXCatEntryArray(end).target)
    error(['I can not find a target in the collision:\n%s\nthat match the regular expression for a state.\n'...
      'Please check your LXCat files'], [description.reactants description.direction description.products]);
  end
  electronArray = regexp(description.reactants, Parse.electronRegExp, 'names');
  numElectrons = 0;
  for i = 1:length(electronArray)
    if isempty(electronArray(i).quantity)
      numElectrons = numElectrons+1;
    else
      numElectrons = numElectrons+str2double(electronArray(i).quantity);
    end
  end
  LXCatEntryArray(end).reactantElectrons = numElectrons;
  productArray = regexp(description.products, Parse.stateRegExp, 'names');
  LXCatEntryArray(end).productArray = removeDuplicatedStates(productArray);
  if isempty(LXCatEntryArray(end).productArray)
    error(['I can not find a product in the collision:\n%s\nthat match the regular expression for a state.\n'...
      'Please check your LXCat files'], [description.reactants description.direction description.products]);
  end
  electronArray = regexp(description.products, Parse.electronRegExp, 'names');
  numElectrons = 0;
  for i = 1:length(electronArray)
    if isempty(electronArray(i).quantity)
      numElectrons = numElectrons+1;
    else
      numElectrons = numElectrons+str2double(electronArray(i).quantity);
    end
  end
  LXCatEntryArray(end).productElectrons = numElectrons;
  if isempty(parameter)
    LXCatEntryArray(end).threshold = 0;
  else
    LXCatEntryArray(end).threshold = str2double(parameter.threshold);
  end
  LXCatEntryArray(end).rawCrossSection = rawCrossSection;
  
end

function newStateArray = removeDuplicatedStates(stateArray)
  
  if isempty(stateArray)
    newStateArray = struct.empty;
    return
  end
  
  numStates = length(stateArray);
  newStateArray = stateArray(1);
  if isempty(newStateArray.quantity)
    newStateArray.quantity = 1;
  else
    newStateArray.quantity = str2double(newStateArray.quantity);
  end
  for i = 2:numStates
    for j = 1:length(newStateArray)
      if strcmp(stateArray(i).gasName, newStateArray(j).gasName) && ...
          strcmp(stateArray(i).ionCharg, newStateArray(j).ionCharg) && ...
          strcmp(stateArray(i).eleLevel, newStateArray(j).eleLevel) && ...
          strcmp(stateArray(i).vibLevel, newStateArray(j).vibLevel) && ...
          strcmp(stateArray(i).rotLevel, newStateArray(j).rotLevel)
        if isempty(stateArray(i).quantity)
          newStateArray(j).quantity = newStateArray(j).quantity+1;
        else
          newStateArray(j).quantity = newStateArray(j).quantity+str2double(stateArray(i).quantity);
        end
        break;
      end
      if j == length(newStateArray)
        newStateArray(end+1) = stateArray(i);
        if isempty(stateArray(i).quantity)
          newStateArray(end).quantity = 1;
        else
          newStateArray(end).quantity = str2double(stateArray(i).quantity);
        end
      end
    end
  end
  
end
