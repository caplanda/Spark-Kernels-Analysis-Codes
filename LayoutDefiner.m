%% Defines the layout of the camera field of view
clear; clc; close all;
% Note: you can hit backspace to delete a point that was just clicked


%% User Inputs

%Default directory (change this for convenience if needed):
% testDir = 'C:\Users\caplanda\Kernel IR Data\20170420 Igniter Test 3\rotated';
% testDir = '\\depot.engr.oregonstate.edu\users\caplanda\Windows.Documents\My Documents\CIRE Lab\Caplan\Data\Open Air Data\20170413 Igniter Test 1';
% testDir = 'J:\Caplan\Kernel IR Data\2017_11_10';
% testDir = 'C:\Users\Danny\Documents\ResearchIR Data\2017_11_10';
testDir = 'J:\Kernel IR Data\2017_12_19';

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
% Data = rot90(Data);  %%%%%%%%%%%%%%%%%%%%%%% DELETE THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imshow(Data, colormapRange)
a = text(2,-15, '', 'FontSize',12, 'Color','r');

%% Define Window Edges
a.String = 'Select Window Edges: [Top, Right, Bottom, Left]';
[xWin, yWin] = ginputc(4, 'Color', 'g');
a.String = '';

%% Define Mirror Edge
a.String = 'Select Mirror Edge (just above then just below)'; 
[xMir, yMir] = ginputc(2, 'Color', 'g');
a.String = '';

%% Define Igniter Edges
a.String = ...
    'Igniter Edges: los(Bottom, Tip, Top), Reflected(Bottom, Tip, Top)';
[xIgn, yIgn] = ginputc(6, 'Color', 'g');
a.String = '';

%% Round to nearest pixel
xWin = round(xWin); yWin = round(yWin);
xMir = round(xMir); yMir = round(yMir);
xIgn = round(xIgn); yIgn = round(yIgn);

%% Plot the selected areas for visual confirmation:
hold on
windowX = [xWin(4), xWin(2), xWin(2), xWin(4), xWin(4)];
windowY = [yWin(1), yWin(1), yWin(3), yWin(3), yWin(1)];
plot(windowX, windowY,'c');

mirrorX = [xWin(4), xWin(2), xWin(2), xWin(4)];
mirrorY = [yMir(1), yMir(1), yMir(2), yMir(2)];
plot(mirrorX, mirrorY, 'm');

losIgniterX = [xWin(2), xIgn(2), xIgn(2), xWin(2)];
losIgniterY = [yIgn(1), yIgn(1), yIgn(3), yIgn(3)];
plot(losIgniterX,losIgniterY, 'g');

refIgniterX = [xWin(2), xIgn(5), xIgn(5), xWin(2)];
refIgniterY = [yIgn(4), yIgn(4), yIgn(6), yIgn(6)];
plot(refIgniterX, refIgniterY, 'g');

%% Define active areas and igniter tip location indices
%note, indices are y then x since they are matrix indices not coordinates
refind = [yWin(1), yMir(1), xWin(4), xWin(2)]; %mirror view
losind = [yMir(2), yWin(3), xWin(4), xWin(2)]; %line of sight view
lostip = [mean([yIgn(1),yIgn(3)]), xIgn(2)]; %los igniter tip
reftip = [mean([yIgn(4),yIgn(6)]), xIgn(5)];  %reflected igniter tip
fprintf('Tip x-shift is %0.0f pixels (should be zero)\n',xIgn(2)-xIgn(5));
fprintf('los igniter width is %0.0f pixels \n', yIgn(1)-yIgn(3));
fprintf('Reflected igniter width is %0.0f pixels (should match los)\n',...
    yIgn(4)-yIgn(6));

plot(reftip(2), reftip(1),'rx'); %igniter tip locations on the layout image
plot(lostip(2),lostip(1),'rx');

%% Adjust igniter tip locations to account for cropped image indices
lostip(1) = lostip(1) - yMir(2);
lostip(2) = lostip(2) - xWin(4);
reftip(1) = reftip(1) - yWin(1);
reftip(2) = reftip(2) - xWin(4);

%% Plot los and mir images for visual confirmation
% figure
% imshow(Data(refind(1):refind(2), refind(3):refind(4)), colormapRange);
% hold on
% plot(reftip(2),reftip(1), 'rx');
% figure
% imshow(Data(losind(1):losind(2), losind(3):losind(4)), colormapRange);
% hold on
% plot(lostip(2),lostip(1), 'rx');

%% Save the offsets

layoutfileused = fileName;
save( strcat(testDir,'\','Layout.mat'), 'refind', 'losind', 'lostip', 'reftip', 'layoutfileused');
saveas(fig1, strcat(filePath, 'Layout.tif'));

fprintf('\nScript complete.');
toc