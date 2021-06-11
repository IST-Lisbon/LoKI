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

% Script to compare two EEDFs

% obtain first EEDF
uiwait(msgbox('Select first EEDF file.'));
[file1, path1] = uigetfile('*.txt');
while ~isequal(file1, 'eedf.txt')
  uiwait(warndlg('Please select an EEDF file (''eedf.txt'')','Warning'));
  [file1, path1] = uigetfile('*.txt');
end
fid1 = fopen([path1 file1], 'r');
fgetl(fid1);
data1 = (fscanf(fid1, '%f', [3 inf]));
fclose(fid1);
legend1 = inputdlg('Define legend for first EEDF:');

% obtain second EEDF
uiwait(msgbox('Select second EEDF file.'));
[file2, path2] = uigetfile('*.txt');
while ~isequal(file2, 'eedf.txt')
  uiwait(warndlg('Please select an EEDF file (''eedf.txt'')','Warning'));
  [file2, path2] = uigetfile('*.txt');
end
fid2 = fopen([path2 file2], 'r');
fgetl(fid2);
data2 = (fscanf(fid2, '%f', [3 inf]));
fclose(fid2);
legend2 = inputdlg('Define legend for second EEDF:');

% plot both EEDFs
figure;
semilogy(data1(1,:), data1(2,:), 'r-', data2(1,:), data2(2,:), 'b-');
xlabel('Energy (eV)');
ylabel('Distribution Function (eV^{-3/2})');
legend(gca, {legend1{1} legend2{1}});