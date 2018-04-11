%% Hot air temperature deconvolution
clear; clc; close all; format compact; format shortg;

%% User Inputs

testDir = 'J:\Kernel IR Data\1_26_2018 HeatGun';
datadir = [testDir,'\Datasets'];

camcalFile = 'TEST CAL FILE - NOT ACCURATE.txt';
BGfile = '';
TestSheetFile = 'TestMatrix.txt';

PressureCol = 3; %Test Sheet Pressure Column
VoltageCol = 4; %Test Sheet Height Column

%% Choose what to process:

% %Choose what to process using GUI:
[dataName, dataPath] = uigetfile('*.mat',...
    'Select the .mat file containing the video',...
    [datadir, '\MatData'], 'MultiSelect', 'on');
dataName = cellstr(dataName);

%%%%%% OR %%%%%%%

% % % Process all files in the folder
% dataName = cellstr(ls(datadir));
% dataName = dataName(3:end); %get rid of stupid dots added by ls


%% Load files
CamCal = dlmread([testDir, '\', camcalFile]);
TestSheet = dlmread([testDir, '\', TestSheetFile]);
% BGint = load([testDir, '\', BGfile]);

%Curve fit calibration data, using linear fit:
fit = fit(CamCal(:,1), CamCal(:,2)*10000, 'poly1'); % (*10000 to convert from W/cm^2-sr to W/m^2-sr)
calfit = @(x) fit.p1*x + fit.p2;

%% Folder Prep:
IntDataDir = [testDir,'\Intensity_Data'];
if exist(IntDataDir) ~= 7
    mkdir(IntDataDir);
end


%% Iterate over folder contents:
RedBaloons = 1;

dataMatFile = matfile([datadir,'\' dataName{RedBaloons}]); %make data vars accessible
Pressure = TestSheet(RedBaloons,PressureCol); %Number of atmospheres
Voltage = TestSheet(RedBaloons,VoltageCol);
DP = dataName{RedBaloons}(1:regexp(dataName{RedBaloons},'.mat')-1);
fprintf(['Oooh, ',DP,'!  Om nom nom... \n']);

%% Convert to intensity and Background Subtract
eventcount = 1;
eventframes = dataMatFile.nFrames;
KernCounts = dataMatFile.Data;
BGint = zeros(size(KernCounts(:,:,1)));

%Convert from raw photon counts to intensity, and background subtract:
for i = 1:eventcount
%     BGint{i,1} = calfit(BGmeanCounts{i,1});
    for j = 1:eventframes(i)
        KernInt{i,j} = calfit(double(KernCounts(:,:,j))) - BGint;
    end
end


%% Save in same format as Kernel_Isolator_v2 (with only one event with many frames)

KernData.InputIsolator.CamCal = camcalFile;
KernData.Pressure = Pressure;
KernData.Voltage = Voltage;
KernData.Intensity = KernInt;
KernData.EventCount = eventcount;
KernData.EventFrames = eventframes;
KernData.DP = DP;
KernData.fps = dataMatFile.fps;
% KernData.time = time;
% KernData.framenum = framenum;
% KernData.Notes = Notes;
% KernData.SparkRate = SparkRate;

save([IntDataDir,'\',DP,'_Int.mat'], 'KernData');

fprintf('\nThat was Tasty!\n\n\n');



