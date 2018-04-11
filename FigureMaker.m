%% Spark Kernels Image Figure Maker
% Creates a vertical stack of cropped images with a colorbar.
clear; clc; close all; format compact; format shortg;

%% User Defined Values

testDir = 'E:\Kernel IR Data\2017_11_20';

DataDir = [testDir,'\Temperature_Data'];
colorange = [0, 15];
map = parula;
% colorange = [0, 3000];
% map = hot;

crop = [10, 50, 56, 79]; %[upperleftx, upperlefty, lowerrightx, lowerrighty]
EventNumber = 2;
spacing = 2; %pixel spacing between images

NumImgs = 4; %number of stacked images
StrStrt = -40; %Starting x location for time labels



%% Load file
[fileName, filePath] = uigetfile('*.mat', ...
    'Select the .mat file containing data to be imaged', DataDir);

load([filePath, fileName]);

%% create figure

ImgH = crop(4) - crop(2);
ImgW = crop(3) - crop(1);

ImageH = NumImgs*ImgH + (NumImgs-1)*spacing;
ImageW = ImgW+1;

Image = zeros(ImageH, ImageW);

y = 1;
for j = 1:NumImgs
    try %%%%%%%% Enter KernData.[what you want to plot] on the next line %%%%%%%%%%%%%%%%%%%
        Image(y:y+ImgH,1:ImageW) = KernData.intLOS{EventNumber,j}(crop(2):crop(4),crop(1):crop(3));
        y = y + ImgH + spacing;
    catch ErrImgBuilder
        ErrImgBuilder
    end
end

fig1 = figure;
imshow(Image,'Colormap', map, 'DisplayRange', colorange); hold on
h = colorbar;

y = 1 + ImgH;
for j = 1:NumImgs-1
    rectangle('Position', [0,y,ImageW+1,spacing], 'FaceColor', 'w');
    y = y + ImgH + spacing;
end

y = 1;
for j = 1:NumImgs
    try
        a = text(StrStrt,y+ImgH/2, '', 'FontSize',12, 'Color','k');
        %         a.String = strcat('t = ',KernData.time{j},'ms');
        a.String = sprintf('t = %0.2f ms',KernData.time{EventNumber,j});
        y = y + ImgH + spacing;
    catch
    end
end

title('1.0 atm');
ylabel(h, 'W/m^2-sr');
% set(get(h,'title'),'string','W/m^2-sr');




