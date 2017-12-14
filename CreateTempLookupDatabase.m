%% Create temperature lookup database using RADCAL_Add_Paths

% be sure to enter the correct number of path lengths and input wavelengths
% in the text files in the radcal directory

clear; clc; close all; format compact; format shortg;

%% Inputs

% Enter database Optical Thickness (path) and temperature ranges:
PathRange = .01:.01:.05;
TempRange = 300:50:1000;
wavelengths = linspace(1,5,100)';

% Enter database .mat file output name:
SaveName = 'Database_2017_11_20.mat';

%% Air Composition:
f_CH4 = 0;
f_H2O = 0.01165;
f_CO2 = 0.00041856;
f_O2 = 0.2095;
f_Ar = 0.0093;
f_N2 = 1 - f_CH4 - f_H2O - f_CO2 - f_O2 - f_Ar;
%radcal inputs should probably get modified to partial pressures eventually...

%% Load files:

% Select response profile files:
defaultDirectory = 'Z:\Windows.Documents\My Documents\CIRE Lab\Caplan\Codes\IR Response Profiles';
[profName, profPath] = uigetfile('*.txt', 'Select response profiles (hold ctrl to select multiple)',...
    defaultDirectory, 'MultiSelect', 'on');

% RADCAL working directory:
path_RADCAL = uigetdir('C:\Users\caplanda\RADCAL', 'Select RADCAL working directory for iteration');

path_MATLAB = cd(path_RADCAL); %change directories to the radcal working directory (must be a mapped folder, not UNC)
addpath(path_MATLAB); % Add the matlab working folder to the current workspace

%% Build 2D lookup table:

[PathLength,Temp] = meshgrid(PathRange,TempRange);
dlmwrite([path_RADCAL, '\input_wavelengths_dan.txt'], wavelengths);
h = waitbar(0,'Building Database...');
tic
for i = 1:length(TempRange)
    for j = 1:length(PathRange)
        %Create and overwrite new RADCAL input file:
        inputRADCAL = [PathLength(i,j), f_CH4, f_H2O, f_CO2, f_N2, f_O2, Temp(i,j)];
        dlmwrite([path_RADCAL, '\input.txt'], inputRADCAL);
        
        system('RADCAL.exe'); %run RADCAL
        
        forward = dlmread([path_RADCAL, '\forward.txt']); %Load radcal output wavelenths file
        backward = FilterApplyIter(forward, profName, profPath, 0, 0); %Apply response profiles
        intensity(i,j) = trapz(backward(:,1),backward(:,2)); %Integrate to get total intensity
    end
    waitbar(i/length(TempRange))
end
toc
close(h)

% Plot results:
surf(Temp,PathLength,intensity)
xlabel('Temperature (K)');
ylabel('Path Length (m)');
zlabel('Integrated Intensity (W/m^2-sr)');

Composition = {'CH4', f_CH4; 'H2O', f_H2O; 'CO2', f_CO2; 'N2', f_N2; 'O2', f_O2}

% Save file:
save(SaveName,'intensity','Temp','PathLength','Composition');
% use Tempvalue = griddata(intensity,PathLength,Temp,[intensityvalue],[pathlengthvalue])

