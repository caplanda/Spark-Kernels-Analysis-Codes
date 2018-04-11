%% Dataset Cropper
clear; clc; close all;
try close(h); catch; end


% Note: this needs to be processed on all dataset files and the layout
% snapshot file prior to running the LayoutDefiner code.

testDir = 'J:\Kernel IR Data\1_26_2018 HeatGun';

[fileName, filePath] = uigetfile('*.mat', ...
    'Select the .mat file containing the crop indicies', [testDir]);

%% Folder Prep:
DataDir = [testDir,'\Datasets'];
if exist(DataDir) ~= 7
    mkdir(DataDir);
end


tic
fprintf('Loading file... \n\n');
load([filePath, fileName]);

%% Select and load file

% %Choose what to process using GUI:
[dataName, dataPath] = uigetfile('*.mat',...
    'Select the .mat file containing the video',...
    [testDir, '\MatData'], 'MultiSelect', 'on');
dataName = cellstr(dataName);

h = waitbar(0,'Processing...');
RedBaloons = 99; keepVars = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
for RedBaloons = 1:length(dataName)
    try %overall try-catch
        
        dataMatFile = load([dataPath, '\', dataName{RedBaloons}]); %make data vars accessible
        
        Data = dataMatFile.Data;
        Data = rot90(Data); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rotate %%%%
        
        % cropind = uint8(cropind);
        Data = Data(cropind(1):cropind(2), cropind(3):cropind(4),:);
        fName = dataMatFile.fName;
        fPath = dataMatFile.fPath;
        fps = dataMatFile.fps;
        h = dataMatFile.h;
        maximum = dataMatFile.maximum;
        nFrames = dataMatFile.nFrames;
        w = dataMatFile.w;
        
        save( strcat(testDir,'\Datasets\',dataName{RedBaloons}), 'Data', 'fName',...
            'fPath', 'fps', 'h', 'maximum', 'nFrames', 'w');
        
    catch Error_Overall %overall try-catch
        fprintf(strcat('Analysis failed on:', dataName{RedBaloons}, ' with the following error:\n'));
        Error_Overall
        fprintf('Comment out the overall try,catch,end structure to get error line numbers...');
    end %overall try-catch
    clearvars('-except',keepVars{:}); %Clear all loop variables
    waitbar(RedBaloons/length(dataName))
end
toc
try
    close(h)
catch
end