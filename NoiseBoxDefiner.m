%% Defines the image box used for determining the noise levels
clear; clc; close all;


%% User Inputs

testDir = 'C:\Users\caplanda\2017_12_22b';

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
imshow(Data, colormapRange)
a = text(2,-15, '', 'FontSize',12, 'Color','r');

%% Define Window Edges
a.String = 'Select empty space area corners: NW then SE]';
[xBox, yBox] = ginputc(2, 'Color', 'g');
a.String = '';

%% Round to nearest pixel
xBox = round(xBox); yBox = round(yBox);


%% Plot the selected areas for visual confirmation:
hold on
boxX = [xBox(1), xBox(2), xBox(2), xBox(1), xBox(1)];
boxY = [yBox(1), yBox(1), yBox(2), yBox(2), yBox(1)];

plot(boxX, boxY, 'r');

%% Save the offsets

layoutfileused = fileName;
save( strcat(testDir,'\','NoiseBox.mat'), 'xBox', 'yBox', 'layoutfileused');
saveas(fig1, strcat(filePath, 'NoiseBox.png'));

fprintf('\nScript complete.');
toc