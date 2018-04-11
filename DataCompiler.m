%% Data Plotter

clear; clc; close all; format compact; format shortg;

%% User Defined Values

testDir = 'J:\Kernel IR Data\1_26_2018 HeatGun';

TempDataDir = [testDir,'\Temperature_Data'];

%% Choose what to process:

%Enter names of files to be processed
% dataName = {'DP-000006_Int.mat'};

%%%%%% OR %%%%%%%

%Process all files in the folder
dataName = cellstr(ls(TempDataDir));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Chew through all files in folder:

h = waitbar(0,'Chewing through folder contents...');
k = 0; Data.DP(1) = 0;
RedBaloons = 99; LastError = 0; keepVars = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
for RedBaloons = 1:length(dataName)
    try %overall try-catch
        load([TempDataDir,'\', dataName{RedBaloons}]);
        fprintf(['Oooh, ',KernData.DP,'!  Om nom nom... \n']);
        fprintf('Looks like there are %d events with up to %d frames.\n',KernData.EventCount,max(KernData.EventFrames));
        
        % Build data structure, spreadsheet style:
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames
                if ~isempty(KernData.Intensity{i,j})
                    k = k + 1;
                    Data(k).DP = KernData.DP;
                    Data(k).Pressure = KernData.Pressure;
                    Data(k).Voltage = KernData.Voltage;
                    Data(k).Intensity = KernData.Intensity{i,j};
                    Data(k).Time = KernData.time{i,j};
                    Data(k).Temperature = KernData.temperature{i,j};
                    Data(k).Tmean = KernData.Tmean{i,j};
                    Data(k).Tpeak = max(KernData.temperature{i,j}(:));
                    Data(k).fps = KernData.fps;
                    Data(k).binLOS = KernData.binLOS{i,j};
                    Data(k).binREF = KernData.binREF{i,j};
                    Data(k).intLOS = KernData.intLOS{i,j};
                    Data(k).intREF = KernData.intREF{i,j};
                    Data(k).PixSizeLOS = KernData.PixSizeLOS;
                    Data(k).PixSizeREF = KernData.PixSizeREF;
                    
                    %Calculate and store mean kernel intensities:
                    intLOSkern = KernData.binLOS{i,j} .* KernData.intLOS{i,j};
                    Data(k).intLOSavg = sum(intLOSkern(:)) / length(find(intLOSkern));
                    intREFkern = KernData.binREF{i,j} .* KernData.intREF{i,j};
                    Data(k).intREFavg = sum(intREFkern(:)) / length(find(intREFkern));
                    
                    %                     Data(k).framenum = KernData.time{i,j} / 1000 * KernData.fps;
                end
            end
        end
        
        fprintf('\nThat was Tasty!\n\n\n');
        
    catch Error_Overall %overall try-catch
        fprintf(strcat('Analysis failed on:', dataName{RedBaloons}, ' with the following error:\n'));
        LastError = Error_Overall
        fprintf('Comment out the overall try,catch,end structure to get error line numbers...');
    end %overall try-catch
    clearvars('-except',keepVars{:}); %Clear all loop variables
    waitbar(RedBaloons/length(dataName))
end
close all
close(h)


DataTable = struct2table(Data);

pressures = unique(extractfield(Data,'Pressure'));
voltages = unique(extractfield(Data,'Voltage'));
times = unique(extractfield(Data,'Time'));

%% Average Plots over time

for p = 1:length(pressures)
    i = 0;
    isobar(p).Pressure = pressures(p);
    isobar(p).LegendEntry = ['P=',num2str(pressures(p)),' bar'];
    for k = 1:size(Data,2)
        if Data(k).Pressure == pressures(p) && ~isnan(Data(k).Tmean)
            i = i+1;
            isobar(p).Tmean(i) = Data(k).Tmean;
            isobar(p).Time(i) = Data(k).Time;
            isobar(p).Tpeak(i) = Data(k).Tpeak;
            isobar(p).intLOSavg(i) = Data(k).intLOSavg;
            isobar(p).area(i) = sum(Data(k).binLOS(:));
            isobar(p).PixSizeLOS(i) = Data(k).PixSizeLOS;
        end
    end
end



figure(1); figure(2); figure(3);
lnspc = {'xr:', '+r:', 'or:', 'xg:', '+g:', 'og:', 'xb:', '+b:', 'ob:'};
for p = 1:length(pressures)
    isobar(p).TimeAvg = times;
    for t = 1:length(times)
        i = 0;
        for j = 1:length(isobar(p).Time)
            if isobar(p).Time(j) == times(t)
                i = i+1;
                TmeanVar(i) = isobar(p).Tmean(j);
                %                 TpeakVar(i) = isobar(p).Tpeak(j);
                intLOSavgVar(i) = isobar(p).intLOSavg(j);
                areaVar(i) = isobar(p).area(j) * (isobar(p).PixSizeLOS(j)*100)^2;
            end
        end
        if exist('TmeanVar')
            isobar(p).TmeanAvg(t) = mean(TmeanVar);
        else
            isobar(p).TmeanAvg(t) = NaN;
        end
        
        if exist('intLOSavgVar')
            isobar(p).intLOSavgAvg(t) = mean(intLOSavgVar);
        else
            isobar(p).intLOSavgAvg(t) = NaN;
        end
        
        if exist('areaVar')
            isobar(p).areaAvg(t) = mean(areaVar);
        else
            isobar(p).areaAvg(t) = NaN;
        end
        clear TmeanVar intLOSavgVar areaVar
    end
    
    figure(1)
    plot(isobar(p).TimeAvg, isobar(p).TmeanAvg, lnspc{p}); hold on
%     title('Average Mean Temperature');
    xlabel('Time (ms)');
    ylabel('Temperature (K)');
    grid on
    axis([0 4 700 1500]);
%     xticks(0.68*[0 1 2 3 4 5 6]);
    
    figure(2)
    plot(isobar(p).TimeAvg, isobar(p).intLOSavgAvg, lnspc{p}); hold on
%     title('Average Mean Intensity');
    xlabel('Time (ms)');
    ylabel('Intensity (W/m^2-sr)');
    grid on
    axis([0 4 0 3.5]);
%     xticks(0.68*[0 1 2 3 4 5 6]);
    
    figure(3)
    plot(isobar(p).TimeAvg, isobar(p).areaAvg, lnspc{p}); hold on
%     title('Average Kernel Area');
    xlabel('Time (ms)');
    ylabel('Area (cm^2)');
    grid on
    axis([0 4 0 3.5]);
%     xticks(0.68*[0 1 2 3 4 5 6]);
end
figure(1)
legend(isobar.LegendEntry)
figure(2)
legend(isobar.LegendEntry)
figure(3)
legend(isobar.LegendEntry)
figs = 3;


%% Average Temperature Map
for p = 1:length(pressures)
    rows = DataTable.Pressure == pressures(p);
    Temps = DataTable.Temperature(rows);
    for i = 1:size(Temps,1)
        TempsMat(:,:,i) = Temps{i};
    end
    AvgTemp(:,:,p) = mean(TempsMat,3,'omitnan');
end

colorange = [300, 3000];
for i = 1:length(pressures)
    figs = figs + 1;
    figure(figs)
    imshow(AvgTemp(:,:,i),'Colormap', hot, 'DisplayRange', colorange)
    title(['Pressure = ', num2str(pressures(i)), ' atm']);
    colorbar
end



%% Phase Average Temperatures
for p = 1:length(pressures)
    figs = figs + 1;
    figure(figs)
    for t = 1:length(times)
        rows = DataTable.Pressure == pressures(p) & DataTable.Time == times(t);
        Temps = DataTable.Temperature(rows);
        for i = 1:size(Temps,1)
            TempsMat(:,:,i) = Temps{i};
        end
        if exist('TempsMat')
        AvgTemp(:,:,t,p) = mean(TempsMat,3,'omitnan');
        subplot(1,length(times),t)
        imshow(AvgTemp(:,:,t,p),'Colormap', hot, 'DisplayRange', colorange)
        title([num2str(pressures(p)), ' atm, ', num2str(times(t)), ' ms']);
        colorbar
        clear rows Temps TempsMat
        end
    end
end

    


%%
fprintf('All Done!\n\n');
toc