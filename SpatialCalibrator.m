%% SpatialCalibrator
% This script will load a .mat image file (converted using sfimgReadLM.m)
% of a calibration object (probably a ruler) and calculate the distance and 
% area that each pixel represents. This will be saved a txt file.

% Inputs: .mat file of a ruler or other calibration object
% Outputs: .txt file with the pixel size, in mm

clear; clc; close all;

%% User Inputs
testDir = 'J:\Kernel IR Data\1_26_2018 HeatGun';


%% Graphical file selector

snapdir = '\Snapshots';

[fileName, filePath] = uigetfile('*.mat', 'Select IR spatial calibration image', [testDir, snapdir])
switch fileName
    case 0
        msgbox('No vegan diet, no vegan powers!','Oops...');
        clear
        return
    otherwise
        eval(['addpath ''' filePath '''']);
        fileName=fileName(1:regexp(fileName,'.mat')-1);
end

%% Load the data file
matFile = matfile([fileName, '.mat']); %make the variables in the mat file accessible
theData = matFile.Data; %Load the Data matrix
colormapRange = [min(theData(:)), max(theData(:))];
imshow(theData, colormapRange); %show the data as an image

%% Calibration: Select two points on the ruler
% Try to pick the same spot on two different tick marks. 
% Pixelation error can be minimizing by choosing points that are far apart from each other.

[x, y] = ginputc(2,'Color', 'g'); % get x and y coordinates of two different tick marks
close Figure 1
distance = input('Enter the distance between the selected points in mm: ');
pixSize = distance/sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 ); % mm per pixel
pixSize = pixSize / 1000; % meters per pixel
% pixArea = pixSize^2 %Pixel area (m^2)

fprintf('Pixel size is %0.3e m/pixel\n',pixSize);

%Save name:
REF = 'REF'; LOS = 'LOS';
calname = input('Enter name of view being calibrated (LOS or REF): ');

fig1 = figure;
imshow(theData, colormapRange); hold on
plot(x,y,'-+g');

%% Save
dlmwrite([testDir,'\SpatialCal_',calname,'.txt'], pixSize); %Save the pixel size as a text file
saveas(fig1, strcat(filePath,fileName,'_',calname, '.tif'));
