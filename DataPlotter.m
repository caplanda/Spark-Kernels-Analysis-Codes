%% Data Plotter

clear; clc; close all; format compact; format shortg;
tic

%% User Defined Values

DataName = 'J:\Kernel IR Data\CompiledDataSciTech\DataTable_Combined.mat';

load(DataName);
tablerows = size(DataStr,1);

%% Align igniter tips and match dimensions

% Find the distances from the igniter tip to all edges for each frame:
distances = zeros(tablerows,4);
for i = 1:tablerows
    imgsize(i,:) = size(DataStr(i).Temperature); %native image size (for both temperature and los intensity)
    distances(i,1) = DataStr(i).IgTip(2); % left edge to tip
    distances(i,2) = imgsize(i,2) - DataStr(i).IgTip(2); % tip to right edge
    distances(i,3) = DataStr(i).IgTip(1); % top edge to tip
    distances(i,4) = imgsize(i,1) - DataStr(i).IgTip(1); % tip to bottom edge
end
distances = round(distances);

% Determine the tip to edge distances that will work for all frames:
Dist(1) = max(unique(distances(:,1)));
Dist(2) = max(unique(distances(:,2)));
Dist(3) = max(unique(distances(:,3)));
Dist(4) = max(unique(distances(:,4)));

% Find locations for the upper left corner of each image within the new image:
NWcorner = zeros(tablerows,2);
for i = 1:tablerows
    NWcorner(i,1) = Dist(3) - distances(i,3);
    NWcorner(i,2) = Dist(1) - distances(i,1);
end

imgSizeNew = [Dist(4)+Dist(3), Dist(1)+Dist(2)];
for i = 1:tablerows
    DataStr(i).ResizeTemp = zeros(imgSizeNew);
    DataStr(i).ResizeTemp(NWcorner(i,1)+1:NWcorner(i,1)+imgsize(i,1),...
        NWcorner(i,2)+1:NWcorner(i,2)+imgsize(i,2)) = DataStr(i).Temperature;
    
    DataStr(i).ResizeInt = zeros(imgSizeNew);
    DataStr(i).ResizeInt(NWcorner(i,1)+1:NWcorner(i,1)+imgsize(i,1),...
        NWcorner(i,2)+1:NWcorner(i,2)+imgsize(i,2)) = DataStr(i).intLOS;
end

%%

DataTable = struct2table(DataStr);

pressures = unique(DataTable.Pressure);
voltages = unique(DataTable.Voltage);
times = unique(DataTable.Time);

PixThresh = 10; %Number of times a pixel needs to be active to be counted in average maps

%% Average Temperature Map
for p = 1:length(pressures)
    rows = DataTable.Pressure == pressures(p);
    Temps = DataTable.ResizeTemp(rows);
    for i = 1:size(Temps,1)
        TempsMat(:,:,i) = Temps{i};
    end
    TempPix = TempsMat > 0;
    PixCounts(:,:,p) = sum(TempPix,3);
    PixUse = PixCounts(:,:,p) > PixThresh;
    AvgTemp(:,:,p) = mean(TempsMat,3,'omitnan');
    AvgTemp = AvgTemp .* PixUse;
    clear TempsMat TempPix PixUse
end

figs = 0;
colorange = [400, 1800];
for i = 1:length(pressures)
    figs = figs + 1;
    figure(figs)
    imshow(AvgTemp(:,:,i),'Colormap', hot, 'DisplayRange', colorange)
    title(['Pressure = ', num2str(pressures(i)), ' atm']);
    c = colorbar;
    ylabel(c, 'Kelvin');
end

colorange = [0, 350];
for i = 1:length(pressures)
    figs = figs+1;
    figure(figs)
    imshow(PixCounts(:,:,i), 'Colormap', jet, 'DisplayRange', colorange);
    title(['Pressure = ', num2str(pressures(i)), ' atm']);
    c = colorbar;
    ylabel(c, 'Pixel Counts');
end

%% Phase Average Temperatures
for p = 1:length(pressures)
    figs = figs + 1;
    figure(figs)
    for t = 1:length(times)
        rows = DataTable.Pressure == pressures(p) & DataTable.Time == times(t);
        Temps = DataTable.ResizeTemp(rows);
        for i = 1:size(Temps,1)
            TempsMat(:,:,i) = Temps{i};
        end
        if exist('TempsMat')
            TempPix = TempsMat > 0;
            PixCounts(:,:,t,p) = sum(TempPix,3);
            PixUse = PixCounts(:,:,t,p) > PixThresh;
            AvgTemp(:,:,t,p) = mean(TempsMat,3,'omitnan');
            AvgTemp = AvgTemp .* PixUse;
            
            subplot(1,length(times),t)
            imshow(AvgTemp(:,:,t,p),'Colormap', hot, 'DisplayRange', colorange)
            title([num2str(pressures(p)), ' atm, ', num2str(times(t)), ' ms']);
            c = colorbar;
            ylabel(c, 'Kelvin');
            clear rows Temps TempsMat
        end
    end
end



%%
fprintf('All Done!\n\n');
toc