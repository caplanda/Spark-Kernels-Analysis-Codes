%% Data Plotter

clear; clc; close all; format compact; format shortg;

%% User Defined Values

testDir = 'J:\Kernel IR Data\2017_11_20';

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
toc

%% Make Plots

pressures = unique(extractfield(Data,'Pressure'));
voltages = unique(extractfield(Data,'Voltage'));
times = unique(extractfield(Data,'Time'));


% 
% for v = 1:length(voltages)
%     figure('Name',['Voltage = ', num2str(voltages(v)),' Volts'],'NumberTitle','off')
%     
%     lnspc = {'xr', '+r', 'or', 'xg', '+g', 'og', 'xb', '+b', 'ob'};
%     for p = 1:length(pressures)
%         i = 0;
%         isovolt(v).isobar(p).Pressure = pressures(p);
%         isovolt(v).isobar(p).LegendEntry = ['P=',num2str(pressures(p)),' bar'];
%         for k = 1:size(Data,2)
%             if Data(k).Pressure == pressures(p) && Data(k).Voltage == voltages(v)
%                 i = i+1;
%                 isovolt(v).isobar(p).Tmean(i) = Data(k).Tmean;
%                 isovolt(v).isobar(p).Time(i) = Data(k).Time;
%                 isovolt(v).isobar(p).Tpeak(i) = Data(k).Tpeak;
%             end
%         end
%         subplot(2,2,1)
%         plot(isovolt(v).isobar(p).Time,isovolt(v).isobar(p).Tmean,lnspc{p}); hold on
%         title('Mean Temp over Time');
%         xlabel('Time (ms)');
%         ylabel('Temperature (K)');
%         grid on
%         subplot(2,2,2)
%         plot(isovolt(v).isobar(p).Time,isovolt(v).isobar(p).Tpeak,lnspc{p}); hold on
%         title('Peak Temp over Time');
%         xlabel('Time (ms)');
%         ylabel('Temperature (K)');
%         grid on
%     end
%     
%     lnspc = {'xr:', '+r:', 'or:', 'xg:', '+g:', 'og:', 'xb:', '+b:', 'ob:'};
%     for p = 1:length(pressures)
%         isovolt(v).isobar(p).TimeAvg = times;
%         for t = 1:length(times)
%             i = 0;
%             for j = 1:length(isovolt(v).isobar(p).Time)
%                 if isovolt(v).isobar(p).Time(j) == times(t) && ~isnan(isovolt(v).isobar(p).Tmean(j))
%                     i = i+1;
%                     TmeanVar(i) = isovolt(v).isobar(p).Tmean(j);
%                     TpeakVar(i) = isovolt(v).isobar(p).Tpeak(j);
%                 end
%             end
%             isovolt(v).isobar(p).TmeanAvg(t) = mean(TmeanVar);
%             isovolt(v).isobar(p).TpeakAvg(t) = mean(TpeakVar);
%             isovolt(v).isobar(p).NumKerns = length(TmeanVar);
%          end
%         subplot(2,2,3)
%         plot(isovolt(v).isobar(p).TimeAvg, isovolt(v).isobar(p).TmeanAvg, lnspc{p}); hold on
%         title('Average Mean Temp over Time');
%         xlabel('Time (ms)');
%         ylabel('Temperature (K)');
%         grid on
%         subplot(2,2,4)
%         plot(isovolt(v).isobar(p).TimeAvg, isovolt(v).isobar(p).TpeakAvg, lnspc{p}); hold on
%         title('Average Peak Temp over Time');
%         xlabel('Time (ms)');
%         ylabel('Temperature (K)');
%         grid on
%     end
%     legend(isovolt(v).isobar.LegendEntry)
% end
% 
% 
% figure('Name',['Mean Temperature'],'NumberTitle','off')
% 
% for v = 1:length(voltages)
%     for p = 1:length(pressures)
%         subplot(1,3,v)
%         plot(isovolt(v).isobar(p).TimeAvg, isovolt(v).isobar(p).TmeanAvg, lnspc{p}); hold on
%         title(['Supply = ',num2str(voltages(v)),' V']);
%         xlabel('Time (ms)');
%         ylabel('Temperature (K)');
%         grid on
%         axis([0 6*0.68 700 1300]);
%         xticks(0.68*[0 1 2 3 4 5 6]);
%     end
% end
% legend(isovolt(v).isobar.LegendEntry)
% 
% 
% 
% clear v

% for v = 1:length(voltages)
%     figure('Name',['Voltage = ', num2str(voltages(v)),' Volts'],'NumberTitle','off')

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
        isobar(p).TmeanAvg(t) = mean(TmeanVar);
%         isobar(p).TpeakAvg(t) = mean(TpeakVar);
        isobar(p).intLOSavgAvg(t) = mean(intLOSavgVar);
        isobar(p).areaAvg(t) = mean(areaVar);
    end
    
    figure(1)
    plot(isobar(p).TimeAvg, isobar(p).TmeanAvg, lnspc{p}); hold on
%     title('Average Mean Temperature');
    xlabel('Time (ms)');
    ylabel('Temperature (K)');
    grid on
    axis([0 6*0.68 750 1200]);
    xticks(0.68*[0 1 2 3 4 5 6]);
    
    figure(2)
    plot(isobar(p).TimeAvg, isobar(p).intLOSavgAvg, lnspc{p}); hold on
%     title('Average Mean Intensity');
    xlabel('Time (ms)');
    ylabel('Intensity (W/m^2-sr)');
    grid on
    axis([0 6*0.68 0 3]);
    xticks(0.68*[0 1 2 3 4 5 6]);
    
    figure(3)
    plot(isobar(p).TimeAvg, isobar(p).areaAvg, lnspc{p}); hold on
%     title('Average Kernel Area');
    xlabel('Time (ms)');
    ylabel('Area (cm^2)');
    grid on
    axis([0 6*0.68 0 3.5]);
    xticks(0.68*[0 1 2 3 4 5 6]);
end
figure(1)
legend(isobar.LegendEntry)
figure(2)
legend(isobar.LegendEntry)
figure(3)
legend(isobar.LegendEntry)


%     legend(isovolt(v).isobar.LegendEntry)
% end

%%
fprintf('All Done!\n\n');
