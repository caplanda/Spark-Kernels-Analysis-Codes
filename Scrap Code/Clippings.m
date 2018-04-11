%% Code Clippings


%% TempAnalysis Clippings:


% Separates kernel indices into spark events, and fills in bg images in the
% beginning of each event.

eventcount = 1;
eventframecount = nBG;
for i = 1:length(kerndex)-1
    eventframecount = eventframecount + 1;
    eventdex(eventcount,eventframecount) = {kerndex(i)};
    if kerndex(i+1)-kerndex(i) > 1
        eventframecount = nBG;
        eventcount = eventcount + 1;
    end
end
eventdex(eventcount,eventframecount+1) = {kerndex(end)}; %fix loop indexing problem

%Add in the prior background frame indicies:
for i = 1:eventcount
    for j = 1:nBG
    eventdex{i,j} = eventdex{i,nBG+1} - j;
    end
end
%NOTE: The preBG frames are in reverse order. This shouldn't matter since
%they are going to be averaged anyway. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non-binarized edge detection algorithms:

% Sobel
[edgeSobel, threshSobel] = edge(lostest, 'Sobel', 'both');

% Prewitt
[edgePrewitt, threshPrewitt] = edge(lostest, 'Prewitt', 'both');

% Roberts
[edgeRoberts, threshRoberts] = edge(lostest, 'Roberts', 'both');

% LoG
sigmaLoG = 1;
[edgeLoG, threshLoG] = edge(lostest, 'log', [], sigmaLoG);

% Canny
sigmaCanny = 1;
[edgeCanny, threshCanny] = edge(lostest, 'Canny', [], sigmaCanny);

% Plot:
figure

subplot(2,3,1)
imshow(lostest, [0 0.005])
title('original');

subplot(2,3,2)
imshow(edgeSobel);
title(strcat('Sobel, thresh = ', num2str(threshSobel)));

subplot(2,3,3)
imshow(edgePrewitt);
title(strcat('Prewitt, thresh = ', num2str(threshPrewitt)));

subplot(2,3,4)
imshow(edgeRoberts);
title(strcat('Roberts, thresh = ', num2str(threshRoberts)));

subplot(2,3,5)
imshow(edgeLoG);
title(strcat('LoG, thresh = ', num2str(threshLoG)));

subplot(2,3,6)
imshow(edgeCanny);
title(strcat('Canny, thresh = ', num2str(threshCanny)));

%%

% %Line of sight %%%%%%%%%%%%%%%
% lostestbin = lostest > binThresh;
% lostestbin(:,realtip(2)-1:end) = 0;
% 
% [losedgeSobel, losthreshSobel] = edge(lostestbin, 'Sobel', 'both');
% losoverlaySobel = losedgeSobel + lostest;
% 
% sigmaCanny = 1.2;
% [losedgeCanny, losthreshCanny] = edge(lostestbin, 'Canny', [], sigmaCanny);
% losoverlayCanny = losedgeCanny + lostest;
% 
% sigmaLoG = 1.2;
% [losedgeLoG, losthreshLoG] = edge(lostestbin, 'log', [], sigmaLoG);
% losoverlayLoG = losedgeLoG + lostest;
% 
% % Reflected %%%%%%%%%%%%%%%%%%%%%%%%
% reftestbin = reftest > binThresh;
% reftestbin(:,reftip(2)-1:end) = 0;

% [refedgeSobel, refthreshSobel] = edge(reftestbin, 'Sobel', 'both');
% refoverlaySobel = refedgeSobel + reftest;
% 
% sigmaCanny = 1.2;
% [refedgeCanny, refthreshCanny] = edge(reftestbin, 'Canny', [], sigmaCanny);
% refoverlayCanny = refedgeCanny + reftest;
% 
% sigmaLoG = 1.2;
% [refedgeLoG, refthreshLoG] = edge(reftestbin, 'log', [], sigmaLoG);
% refoverlayLoG = refedgeLoG + reftest;

% % Plot:
% colormaprange = [0, max(lostest(:))];
% colorscheme = 'jet'; %doesn't work....
% 
% % line of sight %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% subplot(2,3,1)
% imshow(lostest, colormaprange); hold on
% title('BG subtracted IR intensity');
% colormap(colorscheme);
% 
% subplot(2,3,2)
% imshow(lostestbin)
% title(strcat('Binarized with threshold = ', num2str(binThresh)));
% colormap(colorscheme);
% 
% subplot(2,3,4)
% imshow(losoverlaySobel, colormaprange);
% title('Sobel');
% colormap(colorscheme);
% 
% subplot(2,3,5)
% imshow(losoverlayCanny, colormaprange);
% title(strcat('Canny with sigma = ',num2str(sigmaCanny)));
% colormap(colorscheme);
% 
% subplot(2,3,6)
% imshow(losoverlayLoG, colormaprange);
% title(strcat('LoG with sigma = ', num2str(sigmaLoG)));
% colormap(colorscheme);
% 
% % reflected %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% colormaprange = [0, max(reftest(:))];
% 
% figure 
% subplot(2,3,1)
% imshow(reftest, colormaprange); hold on
% title('BG subtracted IR intensity');
% colormap(colorscheme);
% 
% subplot(2,3,2)
% imshow(reftestbin)
% title(strcat('Binarized with threshold = ', num2str(binThresh)));
% colormap(colorscheme);
% 
% subplot(2,3,4)
% imshow(refoverlaySobel, colormaprange);
% title('Sobel');
% colormap(colorscheme);
% 
% subplot(2,3,5)
% imshow(refoverlayCanny, colormaprange);
% title(strcat('Canny with sigma = ',num2str(sigmaCanny)));
% colormap(colorscheme);
% 
% subplot(2,3,6)
% imshow(refoverlayLoG, colormaprange);
% title(strcat('LoG with sigma = ', num2str(sigmaLoG)));
% colormap(colorscheme);


%% Determine line of sight depth of each pixel



% losedges = losedgeLoG;
% refedges = refedgeLoG;
% 
% %maximum width of the kernel in line of sight view:
% [losrow, loscol] = find(losedges);
% losleftedge = min(loscol);
% losrightedge = max(loscol);
% losmaxwidth = losrightedge - losleftedge
% 
% %maximum depth of each kernel (max height in the reflected view):
% [refrow, refcol] = find(refedges);
% reftopedge = min(refrow);
% refbottomedge = max(refrow);
% refmaxheight = refbottomedge - reftopedge
% 
% %Determine width of each row:
% for i = losleftedge:losrightedge
%     temp = find(losedges(:,i)); %temporary variable
%     loswidth(i) = max(temp) - min(temp);
% end
% clear temp
% 
% for i = reftopedge:refbottomedge
%     temp = find(refedges(i,:));
%     refheight(i) = max(temp) - min(temp);
% end

%% Gaussian Filter

% sigma = .1; %standard deviation parameter for Gaussian filter
% test = losInt{1,1};
% filttest = imgaussfilt(test,sigma);
% figure
% subplot(1,2,1)
% imshow(test, [0 0.005]);
% title('original');
% subplot(1,2,2)
% imshow(filttest, [0 0.005]);
% title('Gaussian filter');

%% Edge Detection

% close all

% Kernel event and frame to test the code on:
eventnumber = 3;
framenumber = 1;
lostest = losInt{eventnumber,framenumber};
reftest = refInt{eventnumber,framenumber};

% Simple binary theshold method
% binThresh = 0.0001;
binThresh = 0.0015;

lostestbin = lostest > binThresh;
lostestbin(:,realtip(2)-1:end) = 0;

reftestbin = reftest > binThresh/0.85; %0.85 because of reflection losses
reftestbin(:,reftip(2)-1:end) = 0;

losbound = bwboundaries(lostestbin);
refbound = bwboundaries(reftestbin);


losedges = zeros(size(lostest));
for j = 1:size(losbound,1)
    for i = 1:size(losbound{j},1)
        losedges(losbound{j}(i,1), losbound{j}(i,2)) = 1;
    end
end

refedges = zeros(size(reftest));
for j = 1:size(refbound,1)
    for i = 1:size(refbound{j},1)
        refedges(refbound{j}(i,1), refbound{j}(i,2)) = 1;
    end
end

loscolormaprange = [0, max(lostest(:))];
refcolormaprange = [0, max(reftest(:))];
figure
subplot(2,2,1)
imshow(reftest, refcolormaprange);
title('Reflected');
subplot(2,2,2)
imshow(refedges+reftest, refcolormaprange);
subplot(2,2,3)
imshow(lostest, loscolormaprange);
title('Line of sight');
subplot(2,2,4)
imshow(losedges+lostest, loscolormaprange);

%% Make Plots old

        %         for i = 1:KernData.EventCount
        %             for j = 1:KernData.EventFrames
        %                 if ~isempty(KernData.Intensity{i,j})
        %                     k = k + 1;
        %                     Data.DP{k} = KernData.DP;
        %                     Data.Pressure{k} = KernData.Pressure;
        %                     Data.Voltage{k} = KernData.Voltage;
        %                     Data.Intensity{k} = KernData.Intensity{i,j};
        %                     Data.Time{k} = KernData.time{i,j};
        %                     Data.Temperature{k} = KernData.temperature{i,j};
        %                     Data.Tmean{k} = KernData.Tmean{i,j};
        %                     Data.Tpeak{k} = max(KernData.temperature{i,j}(:));
        %                     Data.fps{k} = KernData.fps;
        %                         end
        %             end
        %         end

% for k = 1:length(Data.DP)
%     if Data.Time{k} == 1 * (1000/Data.fps{k}) && Data.Voltage{k} == 20
%         One.DP{k} = Data.DP{k};
%         One.Pressure{k} = Data.Pressure{k};
%         One.Voltage{k} = Data.Voltage{k};
%         One.Intensity{k} = Data.Intensity{k};
%         One.Time = Data.Time{k};
%         One.Temperature{k} = Data.Temperature{k};
%         One.Tmean{k} = Data.Tmean{k};
%         One.Tpeak{k} = Data.Tpeak{k};
%         One.fps{k} = Data.fps{k};
%     elseif Data.Time{k} == 2 * (1000/Data.fps{k}) && Data.Voltage{k} == 20
%         Two.DP{k} = Data.DP{k};
%         Two.Pressure{k} = Data.Pressure{k};
%         Two.Voltage{k} = Data.Voltage{k};
%         Two.Intensity{k} = Data.Intensity{k};
%         Two.Time = Data.Time{k};
%         Two.Temperature{k} = Data.Temperature{k};
%         Two.Tmean{k} = Data.Tmean{k};
%         Two.Tpeak{k} = Data.Tpeak{k};
%         Two.fps{k} = Data.fps{k};
%     elseif Data.Time{k} == 3 * (1000/Data.fps{k}) && Data.Voltage{k} == 20
%         Three.DP{k} = Data.DP{k};
%         Three.Pressure{k} = Data.Pressure{k};
%         Three.Voltage{k} = Data.Voltage{k};
%         Three.Intensity{k} = Data.Intensity{k};
%         Three.Time = Data.Time{k};
%         Three.Temperature{k} = Data.Temperature{k};
%         Three.Tmean{k} = Data.Tmean{k};
%         Three.Tpeak{k} = Data.Tpeak{k};
%         Three.fps{k} = Data.fps{k};
%     elseif Data.Time{k} == 4 * (1000/Data.fps{k}) && Data.Voltage{k} == 20
%         Four.DP{k} = Data.DP{k};
%         Four.Pressure{k} = Data.Pressure{k};
%         Four.Voltage{k} = Data.Voltage{k};
%         Four.Intensity{k} = Data.Intensity{k};
%         Four.Time = Data.Time{k};
%         Four.Temperature{k} = Data.Temperature{k};
%         Four.Tmean{k} = Data.Tmean{k};
%         Four.Tpeak{k} = Data.Tpeak{k};
%         Four.fps{k} = Data.fps{k};
%     elseif Data.Time{k} == 5 * (1000/Data.fps{k}) && Data.Voltage{k} == 20
%         Five.DP{k} = Data.DP{k};
%         Five.Pressure{k} = Data.Pressure{k};
%         Five.Voltage{k} = Data.Voltage{k};
%         Five.Intensity{k} = Data.Intensity{k};
%         Five.Time = Data.Time{k};
%         Five.Temperature{k} = Data.Temperature{k};
%         Five.Tmean{k} = Data.Tmean{k};
%         Five.Tpeak{k} = Data.Tpeak{k};
%         Five.fps{k} = Data.fps{k};
%     elseif Data.Time{k} == 6 * (1000/Data.fps{k}) && Data.Voltage{k} == 20
%         Six.DP{k} = Data.DP{k};
%         Six.Pressure{k} = Data.Pressure{k};
%         Six.Voltage{k} = Data.Voltage{k};
%         Six.Intensity{k} = Data.Intensity{k};
%         Six.Time = Data.Time{k};
%         Six.Temperature{k} = Data.Temperature{k};
%         Six.Tmean{k} = Data.Tmean{k};
%         Six.Tpeak{k} = Data.Tpeak{k};
%         Six.fps{k} = Data.fps{k};
%     elseif Data.Voltage{k} == 15 || Data.Voltage{k} == 10
%     else
%         k
%     end
% end
%
% fig1 = figure;
% plot(cell2mat(One.Pressure),cell2mat(One.Tmean), 'xb'); hold on
% plot(cell2mat(Two.Pressure),cell2mat(Two.Tmean), 'or');
% plot(cell2mat(Three.Pressure),cell2mat(Three.Tmean), '+g');
% plot(cell2mat(Four.Pressure),cell2mat(Four.Tmean), '^m');
% plot(cell2mat(Five.Pressure),cell2mat(Five.Tmean), 'vk');
% % plot(cell2mat(Six.Pressure),cell2mat(Six.Tmean), '*c');
%
% title('Mean Temperature vs Pressure at 20V supply');
% xlabel('Pressure (atm)');
% ylabel('Mean Path Averaged Temperature (K)');
% legend([num2str(round(One.Time,1)),' \mus'], [num2str(round(Two.Time,1)),' \mus'],...
%     [num2str(round(Three.Time,1)),' \mus'], [num2str(round(Four.Time,1)),' \mus'],...
%     [num2str(round(Five.Time,1)),' \mus']);
%
% One.TmeanAvg = mean(cell2mat(One.Tmean));
% One.PressureAvg = mean(cell2mat(One.Pressure));
% Two.TmeanAvg = mean(cell2mat(Two.Tmean));
% Two.PressureAvg = mean(cell2mat(Two.Pressure));
% Three.TmeanAvg = mean(cell2mat(Three.Tmean));
% Three.PressureAvg = mean(cell2mat(Three.Pressure));
% Four.TmeanAvg = mean(cell2mat(Four.Tmean));
% Four.PressureAvg = mean(cell2mat(Four.Pressure));
% Five.TmeanAvg = mean(cell2mat(Five.Tmean));
% Five.PressureAvg = mean(cell2mat(Five.Pressure));
%
% fig2 = figure;
% plot(


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

