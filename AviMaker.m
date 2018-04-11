%%  Creates avi movies of multiple selected kernels
clear; clc; close all;

%% Choose files

testDir = 'C:\Users\Danny\Desktop\SciTech Download\DataPoints\2017_12_19';

vidName = 'Int06'; %Save Name

% DataDir = [testDir,'\Temperature_Data'];
DataDir = testDir;

colorange = [0, 15];
map = parula;
% colorange = [0, 3000];
% map = hot;

FramesBetween = 3; %number of frames between each event in video
aviFramerate = 2;
Delay = 0.2; %Delay time between frames in gif

crop = [10, 50, 56, 79]; %[upperleftx, upperlefty, lowerrightx, lowerrighty]
%Change crop to be relative to igniter tip.....


%% Load file
[fileName, filePath] = uigetfile('*.mat', ...
    'Select the .mat file containing data to be imaged', DataDir,...
    'MultiSelect', 'on');
fileName = cellstr(fileName);
%% create movie matrix

k = 0; Data = struct;
RedBaloons = 99; LastError = 0; keepVars = 0; VidMat = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
clear VidMat;
for RedBaloons = 1:length(fileName)
    load([filePath, fileName{RedBaloons}]);
    
    for i = 1:KernData.EventCount
        for j = 1:KernData.EventFrames
            if ~isempty(KernData.Intensity{i,j})
                k = k+1;
                VidMat(:,:,k) = KernData.Intensity{i,j};
            end
            
        end
        for m = 1:FramesBetween
            k = k+1;
            VidMat(:,:,k) = 0;
        end
    end
    
    
    
    
    clearvars('-except',keepVars{:}); %Clear all loop variables
end


%% Write AVI Video File

% %Make a stupid matlab movie object...
% for i = 1:size(VidMat,3)
%     imshow(VidMat(:,:,i),'Colormap', parula, 'DisplayRange', colorange, 'InitialMagnification','fit');
%     Vid(i) = getframe;
% end
%
% %Write it to a damn file
% v = VideoWriter( strcat(testDir, vidName, '.avi'), 'Uncompressed AVI');
% v.FrameRate = aviFramerate;
% open(v)
% writeVideo(v,Vid)
% close(v)


%% Write animated gif to file

for idx = 1:size(VidMat,3)
    if idx == 1
        imwrite(VidMat(:,:,idx), parula,[testDir, vidName,'.gif'],'gif','LoopCount',Inf,'DelayTime',Delay);
    else
        imwrite(VidMat(:,:,idx), parula,[testDir, vidName,'.gif'],'gif','WriteMode','append','DelayTime',Delay);
    end
end


