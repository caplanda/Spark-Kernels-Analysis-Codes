%% Determines the noise levels for the test day

clear; clc; close all;

testDir = 'E:\Kernel IR Data\2017_12_19';
BoxFile = 'NoiseBox.mat';

%% Load files

IntDataDir = [testDir,'\Intensity_Data'];
load([testDir,'\',BoxFile]);

%% Choose what to process:

% % %Choose what to process using GUI:
% [dataName, dataPath] = uigetfile('*.mat',...
%     'Select the .mat file containing the video',...
%     [IntDataDir, '\MatData'], 'MultiSelect', 'on');
% dataName = cellstr(dataName);

%%%%%% OR %%%%%%%

% % Process all files in the folder
dataName = cellstr(ls(IntDataDir));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Chew through all files in folder:

h = waitbar(0,'Chewing through folder contents...');
RedBaloons = 99; keepVars = 0; Box = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
clear Box
tic
for RedBaloons = 1:length(dataName)
%     try %Overall try-catch
        load([IntDataDir,'\', dataName{RedBaloons}]);
        k = 0;
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                k = k+1;
                Box(:,:,k) = KernData.Intensity{i,j}(yBox(1):yBox(2),xBox(1):xBox(2));
            end
        end
        
%     catch
%     end
end
close(h)
toc

noiseAvg = mean(Box(:));
noiseStd = std(Box(:));

save([testDir,'\Noise.mat'],'noiseAvg', 'noiseStd');


fprintf('All Done!\n\n');





