%% Data Plotter

clear; clc; close all; format compact; format shortg;
tic

%% User Defined Values

DataName = '\\depot.engr.oregonstate.edu\users\caplanda\Windows.Documents\My Documents\CIRE Lab\Caplan\Data\Okhovat Discrepancy\2std_CO2_000625_H2O_01\DataTable_Combined.mat';

timetol = 0.25; %time tollerance for determining if similar times are the 'same'
presstol = 0.03; %pressure tollerance for determining if similar pressures are the 'same'
minKern = 50; %minimum number of kernels at a particular time or pressure for that time/pressure to be included

PixThresh = 10; %Number of times a pixel needs to be active to be counted in average maps

percentile = 10; %percentage to use for percentile maps (not a decimal)

load(DataName);
tablerows = size(DataStr,1);

set(groot,'defaultAxesLineStyleOrder',{'x:','+:','o:','*:', 's:', 'd:', '^:', 'v:', '<:', '>:' });
set(groot,'defaultAxesColorOrder', [0    0.4470    0.7410]);

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
Dist(1) = max(unique(distances(:,1))); %Left edge
Dist(2) = max(unique(distances(:,2))); %Right edge
Dist(3) = max(unique(distances(:,3))); %Top edge
Dist(4) = max(unique(distances(:,4))); %Bottom edge

% Find locations for the upper left corner of each image within the new image:
NWcorner = zeros(tablerows,2);
for i = 1:tablerows
    NWcorner(i,1) = Dist(3) - distances(i,3);
    NWcorner(i,2) = Dist(1) - distances(i,1);
    DataStr(i).IgTipNew(1) = DataStr(i).IgTip(1) + NWcorner(i,1);
    DataStr(i).IgTipNew(2) = DataStr(i).IgTip(2) + NWcorner(i,2);
end

% Overlay all images on the new image size, with tips aligned:
imgSizeNew = [Dist(4)+Dist(3), Dist(1)+Dist(2)];
for i = 1:tablerows
    DataStr(i).ResizeTemp = zeros(imgSizeNew);
    DataStr(i).ResizeTemp(NWcorner(i,1)+1:NWcorner(i,1)+imgsize(i,1),...
        NWcorner(i,2)+1:NWcorner(i,2)+imgsize(i,2)) = DataStr(i).Temperature;
    
    DataStr(i).ResizeInt = zeros(imgSizeNew);
    DataStr(i).ResizeInt(NWcorner(i,1)+1:NWcorner(i,1)+imgsize(i,1),...
        NWcorner(i,2)+1:NWcorner(i,2)+imgsize(i,2)) = DataStr(i).intLOS;
end

% Modify the origin location, so that location axes display properly
PixSize = mean(extractfield(DataStr,'PixSizeLOS'));
xWorldLimits = [-Dist(1), Dist(2)] * PixSize * 100;
yWorldLimits = [-Dist(3), Dist(4)] * PixSize * 100;


%% Play intensity data as movie

for i = 1:size(DataStr,1)
    IntVid(:,:,i) = DataStr(i).ResizeInt;
end

figure
handle = implay(IntVid);

% Set the colormap:
handle.Visual.ColorMap.UserRangeMin = 0;
handle.Visual.ColorMap.UserRangeMax = 15;
handle.Visual.ColorMap.UserRange = 1;


%%

DataTable = struct2table(DataStr);


% pressures = unique(DataTable.Pressure);
% voltages = unique(DataTable.Voltage);





%% Find the unique times within a tolerance, and occuring a minimum number of times

timesRaw = unique(DataTable.Time);
n = 0;
timesTol = uniquetol(timesRaw, timetol/max(timesRaw(:)));
for i = 1:length(timesTol)
    TimeCount = sum(DataTable.Time==timesTol(i));
    if TimeCount > minKern
        n = n+1;
        rows = abs(DataTable.Time - timesTol(i)) <= timetol;
        times(n) = mean(DataTable.Time(rows));
    end
    clear TimeCount
end

%% Find the unique pressures within a tolerance, and occuring a minimum number of times

pressRaw = unique(DataTable.Pressure);
n = 0;
pressTol = uniquetol(pressRaw, presstol/max(pressRaw(:)));
for i = 1:length(pressTol)
    PressCount = sum(DataTable.Pressure==pressTol(i));
    PC(i) = PressCount;
    if PressCount > minKern
        n = n+1;
        rows = abs(DataTable.Pressure - pressTol(i)) <= presstol;
        pressures(n) = mean(DataTable.Pressure(rows));
    end
    clear PressCount
end

%% Average Temperature Map

%Generate average and Percentile Temperature matricies:
for p = 1:length(pressures)
    rows = abs(DataTable.Pressure - pressures(p)) <= presstol; %rows for each pressure
    Temps = DataTable.ResizeTemp(rows); %Temperatures for each pressure
    for i = 1:size(Temps,1)
        TempsMat(:,:,i) = Temps{i}; %Put temperatures together into a matrix
    end
    TempPix = TempsMat > 0; %Pixels that have temperature info
    PixCounts(:,:,p) = sum(TempPix,3); %Map of number of times each pixel is active
    PixUse = PixCounts(:,:,p) > PixThresh; %Map of active pixels that occur more times than the threshold
    AvgTemp(:,:,p) = mean(TempsMat,3,'omitnan'); %Average temperature map for each pressure
    PercentTemp(:,:,p) = prctile(TempsMat,percentile,3);
    AvgTemp = AvgTemp .* PixUse; %Average temperature map displaying only pixels that occur often enough
    PercentTemp = PercentTemp .* PixUse;
    clear TempsMat TempPix PixUse
end

% clear PixCounts PixUse %PixCounts is used in plots later, and should be from intensity not temperature

%Generate average and percentile Intensity matricies:
for p = 1:length(pressures)
    rows = abs(DataTable.Pressure - pressures(p)) <= presstol; %rows for each pressure
    NumKern(p) = sum(rows); %number of kernels for each pressure
    Ints = DataTable.ResizeInt(rows); %Temperatures for each pressure
    for i = 1:size(Ints,1)
        IntMat(:,:,i) = Ints{i}; %Put temperatures together into a matrix
    end
    %     IntPix = IntMat > 0; %Pixels that have temperature info
    %     PixCounts(:,:,p) = sum(IntPix,3); %Map of number of times each pixel is active
    PixUse = PixCounts(:,:,p) > PixThresh; %Map of active pixels that occur more times than the threshold
    AvgInt(:,:,p) = mean(IntMat,3,'omitnan'); %Average temperature map for each pressure
    PercentInt(:,:,p) = prctile(IntMat,percentile,3);
    AvgInt = AvgInt .* PixUse; %Average temperature map displaying only pixels that occur often enough
    PercentInt = PercentInt .* PixUse;
    AvgIntWorld(p) = imref2d(size(AvgInt(:,:,p)),xWorldLimits,yWorldLimits);
    World(p) = imref2d(size(PercentInt(:,:,p)),xWorldLimits,yWorldLimits);
    clear IntMat IntPix PixUse
end

figs = 0;

%Plot the average temperature map
colorange = [400, 1400];
for i = 1:length(pressures)
    figs = figs + 1;
    figure(figs);
    imshow(AvgTemp(:,:,i), World(1), 'Colormap', hot, 'DisplayRange', colorange, 'InitialMagnification','fit')
    axis on; hold on
    plot(DataStr(i).IgTipNew(2), DataStr(i).IgTipNew(1), 'xr');
    title(['Pressure = ', num2str(pressures(i),1), ' atm']);
    xlabel('(cm)');
    ylabel('(cm)');
    c = colorbar;
    ylabel(c, 'Kelvin');
    
    %For sensitivity Analysis:
    if abs(pressures(i)-1) < 0.03
        Avgpeak = prctile(reshape(AvgTemp(:,:,i),[],1),90)
    end
    
    
end

%Plot the active pixel map (based on temperature):
% colorange = [0, 80];
% for i = 1:length(pressures)
%     figs = figs + 1;
%     figure(figs);
%     imshow(PixCounts(:,:,i)/NumKern(i)*100, World(i), 'Colormap', pink, 'DisplayRange', colorange, 'InitialMagnification','fit');
%     title(['Pressure = ', num2str(pressures(i),1), ' atm']);
%     c = colorbar;
%     ylabel(c, 'Normalized Location Frequency %');
% end

% %Plot the percentile temperature map
% colorange = [400, 1200];
% for i = 1:length(pressures)
%     figs = figs + 1;
%     figure(figs);
%     imshow(PercentTemp(:,:,i),World(i),'Colormap', hot, 'DisplayRange', colorange, 'InitialMagnification','fit')
%     axis on; hold on
%     plot(DataStr(i).IgTipNew(2), DataStr(i).IgTipNew(1), 'xr');
%     title(['Pressure = ', num2str(pressures(i),1), ' atm']);
%     xlabel('(cm)');
%     ylabel('(cm)');
%     c = colorbar;
%     ylabel(c, 'Kelvin');
% end

%Plot the percentile intensity map
% colorange = [0, 1.8];
% for i = 1:length(pressures)
%     figs = figs + 1;
%     figure(figs);
%     imshow(PercentInt(:,:,i),World(i), 'Colormap', parula, 'DisplayRange', colorange, 'InitialMagnification','fit')
%     axis on; hold on
%     plot(0,0, 'xr');
%     title(['Pressure = ', num2str(pressures(i),1), ' atm']);
%     xlabel('(cm)');
%     ylabel('(cm)');
%     c = colorbar;
%     ylabel(c, 'W/m^2-sr');
% end

% figure
% for i = 1:size(DataStr,1)
%     plot(DataStr(i).IgTipNew(2), DataStr(i).IgTipNew(1), 'x'); hold on
% end



%% Phase Average Temperatures
colorange = [400, 2000];
ncols = 5;
for p = 1:length(pressures)
    nrows = ceil(length(times)/ncols);
    fig = figure('Position', [25, 50, 1800, 500*nrows]);
    for t = 1:length(times)
        rows = abs(DataTable.Pressure - pressures(p)) <=presstol & abs(DataTable.Time - times(t)) <= timetol;
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
            
            subplot(nrows,ncols,t);
            imshow(AvgTemp(:,:,t,p), World(1),'Colormap', hot, 'DisplayRange', colorange, 'InitialMagnification','fit')
            title([num2str(pressures(p),1), ' atm, ', num2str(times(t),2), ' ms']);
            xlabel('(cm)');
            ylabel('(cm)');
            %             c = colorbar;
            %             ylabel(c, 'Kelvin');
            clear rows Temps TempsMat
            
            %For sensitivity Analysis
            if abs(pressures(p)-1) < 0.03
                if abs(times(t)-0.68) < 0.03
                    phasepeak = prctile(reshape(AvgTemp(:,:,t,p),[],1),90)
                end
            end
            
            
        end
    end
end

% 
% %% Phase Average Intensities
% clear IntMat Ints
% colorange = [0, 10];
% ncols = 5;
% for p = 1:length(pressures)
%     nrows = ceil(length(times)/ncols);
%     fig = figure('Position', [25, 50, 1800, 500*nrows]);
%     for t = 1:length(times)
%         rows = abs(DataTable.Pressure - pressures(p)) <=presstol & abs(DataTable.Time - times(t)) <= timetol;
%         Ints = DataTable.ResizeInt(rows);
%         for i = 1:size(Ints,1)
%             IntMat(:,:,i) = Ints{i};
%         end
%         if exist('IntMat')
%             %             IntPix = IntMat > 0;
%             %             PixCounts(:,:,t,p) = sum(IntPix,3);
%             PixUse = PixCounts(:,:,t,p) > PixThresh;
%             AvgInt(:,:,t,p) = mean(IntMat,3,'omitnan');
%             AvgInt = AvgInt .* PixUse;
%             
%             subplot(nrows,ncols,t);
%             imshow(AvgInt(:,:,t,p), World(1),'Colormap', parula, 'DisplayRange', colorange, 'InitialMagnification','fit')
%             
%             hold on
%             plot(DataStr(i).IgTipNew(2), DataStr(i).IgTipNew(1), 'xr');
%             
%             title([num2str(pressures(p),1), ' atm, ', num2str(times(t),2), ' ms']);
%             xlabel('(cm)');
%             ylabel('(cm)');
%             %             c = colorbar;
%             %             ylabel(c, 'W/m^2-sr');
%         end
%         clear rows Ints IntMat
%     end
% end
% 
% %% Colorbar maker
% 
figure
caxis([400, 2000])
colormap hot
c = colorbar;
c.Label.String = 'Temperature, Kelvin';
% 
% %% Isobar plots
% marker = {'rx:','g+:','bo:','c*:', 'ms:', 'kd:', 'r^:', 'gv:', 'b<:', 'c>:' };
% Font = 12;
% 
% % Avg Temp vs Time:
% clear AvgTemp AvgInt
% figure
% for p = 1:length(pressures)
%     for t = 1:length(times)
%         rows = abs(DataTable.Pressure - pressures(p)) <=presstol & abs(DataTable.Time - times(t)) <= timetol;
%         TempCell = DataTable.Temperature(rows);
%         KernTemp = 0;
%         for i = 1:size(TempCell,1)
%             TempArray = reshape(TempCell{i},[1,numel(TempCell{i})]);
%             KernTemp(i) = mean(TempArray,'omitnan');
%         end
%         AvgTemp(t) = mean(KernTemp,'omitnan');
%         clear TempArray KernTemp TempCell rows
%         if AvgTemp(t) == 0
%             AvgTemp(t) = NaN;
%         end
%     end
%     plot(times, AvgTemp, marker{p}); hold on
%     legendInfo{p} = ['P = ' num2str(pressures(p),1) 'atm'];
%     clear AvgTemp
% end
% legend(legendInfo, 'FontSize', Font);
% xlabel('Time (ms)', 'FontSize', Font);
% ylabel('Average Temperature (K)', 'FontSize', Font);
% set(gca,'fontsize',Font)
% 
% % Area vs Time:
% figure
% for p = 1:length(pressures)
%     for t = 1:length(times)
%         rows = abs(DataTable.Pressure - pressures(p)) <=presstol & abs(DataTable.Time - times(t)) <= timetol;
%         Vol(t) = mean(DataTable.Volume(rows),'omitnan')*(100^3);
%     end
%     plot(times, Vol, marker{p}); hold on
%     legendInfo{p} = ['P = ' num2str(pressures(p),1) 'atm'];
%     clear Vol
% end
% legend(legendInfo, 'FontSize', Font);
% xlabel('Time (ms)', 'FontSize', Font);
% ylabel('Volume (cm^3)', 'FontSize', Font);
% set(gca,'fontsize',Font)
% 
% % Avg Intensity vs Time:
% figure
% for p = 1:length(pressures)
%     for t = 1:length(times)
%         rows = abs(DataTable.Pressure - pressures(p)) <=presstol & abs(DataTable.Time - times(t)) <= timetol;
%         IntCell = DataTable.Intensity(rows);
%         KernInt = 0;
%         for i = 1:size(IntCell,1)
%             IntArray = reshape(IntCell{i},[1,numel(IntCell{i})]);
%             KernInt(i) = mean(IntArray,'omitnan');
%         end
%         AvgInt(t) = mean(KernInt,'omitnan');
%         clear KernInt IntArray IntCell rows
%         if AvgInt(t) == 0
%             AvgInt(t) = NaN;
%         end
%     end
%     plot(times, AvgInt, marker{p}); hold on
%     legendInfo{p} = ['P = ' num2str(pressures(p),1) 'atm'];
%     clear AvgInt
% end
% legend(legendInfo, 'FontSize', Font);
% xlabel('Time (ms)', 'FontSize', Font);
% ylabel('Average Intensity (W/m^2-sr)', 'FontSize', Font);
% set(gca,'fontsize',Font)


%%
fprintf('All Done!\n\n');
toc