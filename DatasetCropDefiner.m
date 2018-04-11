%% Video Crop Selector
clear; clc; close all;

testDir = 'J:\Kernel IR Data\1_26_2018 HeatGun';

%% Select and load file

snapdir = '\Snapshots';

[fileName, filePath] = uigetfile('*.mat', ...
    'Select the .mat file containing the layout image', [testDir,snapdir]);

tic
fprintf('Loading file... \n\n');
load([filePath, fileName]);

colormapRange = [min(Data(:)), maximum];

%% show image
fig1 = figure;
Data = rot90(Data);  %%%%%%%%%%%%%%%%%%%%%%% DELETE THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imshow(Data, colormapRange)
a = text(2,-15, '', 'FontSize',12, 'Color','r');

%% Define Crop Edges
a.String = 'Select Crop Edges: [Top, Right, Bottom, Left]';
[xWin, yWin] = ginputc(4, 'Color', 'g');
a.String = '';

cropind = uint16([yWin(1),yWin(3), xWin(4),xWin(2)]);

%% Plot the selected area for visual confirmation:
hold on
windowX = [xWin(4), xWin(2), xWin(2), xWin(4), xWin(4)];
windowY = [yWin(1), yWin(1), yWin(3), yWin(3), yWin(1)];
plot(windowX, windowY,'c');

figure
imshow(Data(cropind(1):cropind(2), cropind(3):cropind(4)), colormapRange);

%% Save crop indicies:
save( strcat(testDir,'\','CropIndicies.mat'), 'cropind');
saveas(fig1, strcat(filePath, 'Crop.tif'));