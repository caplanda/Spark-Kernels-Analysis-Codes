%% Spark Kernels Temperature Analysis Code
clear; clc; close all; format compact; format shortg;

%% User Defined Values

testDir = 'E:\Kernel IR Data\2017_12_22a';
layoutFile = 'Layout.mat';
REFspatialFile = 'SpatialCal_REF.txt';
LOSspatialFile = 'SpatialCal_LOS.txt';
tempdatabaseFile = 'Database_2017_12_22_MultiPressure.mat';
NoiseFile = 'Noise.mat';

% noiseAvg = 0.034; %Average noise level
% noiseStd = 0.058; %Average noise standard deviation
num_noiseStd = 4; %Threshold for kernel detection (num standard deviations above noise)
artSize = 5; %Maximum artifact size (pixels) for noise reduction
peakNum = 10; %Top ____ temps to be averaged to get the peak temp

Notes = '';

%% Load files

IntDataDir = [testDir,'\Intensity_Data'];
TempDataDir = [testDir,'\Temperature_Data'];

try
    if exist(TempDataDir) ~= 7
        mkdir(TempDataDir);
    end
catch
    fprintf('Check the test directory name...');
end
load([testDir, '\', tempdatabaseFile]);
load([testDir, '\', layoutFile]);
load([testDir,'\', NoiseFile]);
PixSizeREF = dlmread([testDir, '\', REFspatialFile]);
PixSizeLOS = dlmread([testDir, '\', LOSspatialFile]);

%% Choose what to process:

% %Choose what to process using GUI:
% [dataName, dataPath] = uigetfile('*.mat',...
%     'Select the .mat file containing the video',...
%     [IntDataDir, '\MatData'], 'MultiSelect', 'on');
% dataName = cellstr(dataName);

%%%%%% OR %%%%%%%

% %Process all files in the folder
dataName = cellstr(ls(IntDataDir));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Chew through all files in folder:

h = waitbar(0,'Chewing through folder contents...');
RedBaloons = 99; LastError = 0; keepVars = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
for RedBaloons = 1:length(dataName)
    try %Overall try-catch
        load([IntDataDir,'\', dataName{RedBaloons}]);
        fprintf(['Oooh, ',KernData.DP,'!  Om nom nom... \n']);
        fprintf('Looks like there are %d events with up to %d frames.\n',KernData.EventCount,max(KernData.EventFrames));
        
        %% Separate out LOS and REF views:
        
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                intREF{i,j} = KernData.Intensity{i,j}(refind(1):refind(2), refind(3):refind(4));
                intLOS{i,j} = KernData.Intensity{i,j}(losind(1):losind(2), losind(3):losind(4));
            end
        end
        sizeREF = size(intREF{1,1});
        sizeLOS = size(intLOS{1,1});
        
        %% Edge Detection
        % Simple binary theshold method
        
        
        binThresh = noiseAvg + num_noiseStd*noiseStd;
        
        % Kernel event and frame to test the code on:
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                
                %Binarize:
                binLOS{i,j} = intLOS{i,j} > binThresh;
                binLOS{i,j}(:,lostip(2)-1:end) = 0; %Exclude pixels from the igniter tip and over
                binREF{i,j} = intREF{i,j} > binThresh/0.85; %0.85 because of reflection losses
                binREF{i,j}(:,reftip(2)-1:end) = 0;
                
                %Remove artifacts:
                ccLOS = bwconncomp(binLOS{i,j},4);
                ccpixLOS = cellfun(@numel,ccLOS.PixelIdxList);
                artLOS = find(ccpixLOS <= artSize);
                for k = 1:length(artLOS)
                    binLOS{i,j}(ccLOS.PixelIdxList{artLOS(k)}) = 0;
                end
                
                ccREF = bwconncomp(binREF{i,j},4);
                ccpixREF = cellfun(@numel,ccREF.PixelIdxList);
                artREF = find(ccpixREF <= artSize);
                for k = 1:length(artREF)
                    binREF{i,j}(ccREF.PixelIdxList{artREF(k)}) = 0;
                end
                
                %Remove edge reflections
                binLOS{i,j} = imclearborder(binLOS{i,j},4);
                binREF{i,j} = imclearborder(binREF{i,j},4);
                
                %kernel edge locations
                boundLOS = bwboundaries(binLOS{i,j});
                boundREF = bwboundaries(binREF{i,j});
                
                %Organize edges into a matrix
                edgeLOS = zeros(sizeLOS);
                for m = 1:size(boundLOS,1)
                    for n = 1:size(boundLOS{m},1)
                        edgeLOS(boundLOS{m}(n,1), boundLOS{m}(n,2)) = 1;
                    end
                end
                edgesLOS{i,j} = edgeLOS; %store in a cell structure
                
                edgeREF = zeros(sizeREF);
                for m = 1:size(boundREF,1)
                    for n = 1:size(boundREF{m},1)
                        edgeREF(boundREF{m}(n,1), boundREF{m}(n,2)) = 1;
                    end
                end
                edgesREF{i,j} = edgeREF;
                
            end
        end
        
        %% Determine kernel height and depth
        
        % Align the igniter tip in the los and ref images (x matters, y doesn't)
        tipoffset  = abs(lostip(2) - reftip(2)); %igniter tip x-dir offset betweein the los and ref views
        filler = zeros(sizeREF(1),tipoffset);
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                if lostip(2) > reftip(2)
                    temp = cat(2, filler, binREF{i,j}); %add NaN columns on left
                    binREF{i,j} = temp(:,1:end-tipoffset); %crop to original size (right is cropped off)
                elseif lostip(2) < reftip(2)
                    temp = cat(2, binREF{i,j}, filler); %add NaN columns to right
                    binREF{i,j} = temp(:,tipoffset+1:end); %crop to original size (left is cropped off)
                end
            end
        end
        
        % Find left and right edges of the kernel in both views
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                loskerncols = max(binLOS{i,j}, [], 1);
                losleft{i,j} = find(loskerncols,1,'first');
                losright{i,j} = find(loskerncols,1,'last');
                
                refkerncols = max(binREF{i,j}, [], 1);
                refleft{i,j} = find(refkerncols,1,'first');
                refright{i,j} = find(refkerncols,1,'last');
            end
        end
        
        % Eliminate any columns that do not correspond to a column in the other view
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                for k = 1:sizeLOS(2)
                    if sum(binLOS{i,j}(:,k)) == 0
                        binREF{i,j}(:,k) = 0;
                    elseif sum(binREF{i,j}(:,k)) == 0
                        binLOS{i,j}(:,k) = 0;
                    end
                end
            end
        end
        
        % Find left and right edges of the los kernel (now same as ref)
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                kerncols = max(binLOS{i,j}, [], 1);
                leftedge{i,j} = find(kerncols,1,'first');
                rightedge{i,j} = find(kerncols,1,'last');
                %         width{i,j} = rightedge{i,j} - leftedge{i,j};
            end
        end
        
        % Determine the major and minor axes of each elipse slice
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                for k = leftedge{i,j}:rightedge{i,j}
                    if sum(binLOS{i,j}(:,k)) ~= 0
                        lostop{i,j}(k) = find(binLOS{i,j}(:,k),1,'first');
                        losbot{i,j}(k) = find(binLOS{i,j}(:,k), 1, 'last');
                        reftop{i,j}(k) = find(binREF{i,j}(:,k),1,'first');
                        refbot{i,j}(k) = find(binREF{i,j}(:,k), 1, 'last');
                        losaxis{i,j}(k) = losbot{i,j}(k) - lostop{i,j}(k);
                        refaxis{i,j}(k) = refbot{i,j}(k) - reftop{i,j}(k);
                        loscenter{i,j}(k) = (lostop{i,j}(k) + losbot{i,j}(k)) / 2;
                    end
                end
            end
        end
        
        d = @(a,b,x) 2*sqrt((1-x^2/a^2)*b^2); % chord length of an ellipse
        
        % Determine the depth for each los pixel using ellipse approximation
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                depth{i,j} = NaN(sizeLOS);
                for k = leftedge{i,j}:rightedge{i,j}
                    if sum(binLOS{i,j}(:,k)) ~= 0
                        if losaxis{i,j}(k) >= refaxis{i,j}(k)
                            for m = lostop{i,j}(k):losbot{i,j}(k)
                                depth{i,j}(m,k) = d(losaxis{i,j}(k)/2, refaxis{i,j}(k)/2, (m-loscenter{i,j}(k)));
                            end
                        elseif losaxis{i,j}(k) < refaxis{i,j}(k)
                            for m = lostop{i,j}(k):losbot{i,j}(k)
                                depth{i,j}(m,k) = d(refaxis{i,j}(k)/2, losaxis{i,j}(k)/2, (m-loscenter{i,j}(k)));
                            end
                        else
                            fprintf('Error in ellipse loop at: i=%d, j=%d, k=%d, m=%d \n', i,j,k,m);
                        end
                    end
                end
                depth{i,j} = depth{i,j} * PixSizeREF;
            end
        end
        
        % Determine the ellipsoid stack volume
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                volumeVar = 0;
                for k = leftedge{i,j}:rightedge{i,j}
                    if sum(binLOS{i,j}(:,k)) ~= 0
                        a = refaxis{i,j}(k)/2*PixSizeREF;
                        b = losaxis{i,j}(k)/2*PixSizeLOS;
                        volumeVar = volumeVar + (pi*a*b)*PixSizeLOS;
                    end
                end
                Volume{i,j} = volumeVar;
            end
        end
        
        %% Temperature Lookup
        
%         %     h = waitbar(0,'Looking up temperatures...');
%         %     tic
        PressDex = find(Database.Pressures==KernData.Pressure); %Index for which interpolation surface to use
        for i = 1:KernData.EventCount
            for j = 1:KernData.EventFrames(i)
                Tsum{i,j} = 0;
                Tcount{i,j} = 0;
                temperature{i,j} = NaN(sizeLOS);
                for k = leftedge{i,j}:rightedge{i,j}
                    if sum(binLOS{i,j}(:,k)) ~= 0
                        for m = lostop{i,j}(k):losbot{i,j}(k)
                            if binLOS{i,j}(m,k) == 1
                                temperature{i,j}(m,k) = ...
                                    griddata(Database.Int(PressDex).IntGrid, Database.Path, Database.Temp,...
                                    intLOS{i,j}(m,k),depth{i,j}(m,k));
                                if ~isnan(temperature{i,j}(m,k))
                                    Tsum{i,j} = Tsum{i,j} + temperature{i,j}(m,k);
                                end
                                Tcount{i,j} = Tcount{i,j} + 1;
                            end
                        end
                    end
                end
                Tmean{i,j} = Tsum{i,j} / Tcount{i,j};
            end
            %         waitbar(i/KernData.EventCount)
        end
%         %     toc
%         %     close(h)
        
        
        
        %% Save temperature data to the existing kernel data file
        
        KernData.InputTemp.Database = tempdatabaseFile;
        KernData.InputTemp.noiseAvg = noiseAvg;
        KernData.InputTemp.noiseStd = noiseStd;
        KernData.InputTemp.num_noiseStd = num_noiseStd;
        KernData.InputTemp.artSize = artSize;
        KernData.temperature = temperature;
        KernData.Tmean = Tmean;
        KernData.PixSizeREF = PixSizeREF;
        KernData.PixSizeLOS = PixSizeLOS;
        KernData.intREF = intREF;
        KernData.intLOS = intLOS;
        KernData.binLOS = binLOS;
        KernData.binREF = binREF;
        KernData.Volume = Volume;
        KernData.IgTip = lostip;
        KernData.Notes = Notes;
        
        save([TempDataDir,'\', dataName{RedBaloons}],'KernData');
        
        
        %% plot a particular kernel
%         i = 1;
%         j = 1;
%         
%         close all;
%         figure
%         subplot(2,1,1)
%         imshow(binREF{i,j});
%         subplot(2,1,2)
%         imshow(binLOS{i,j});
%         
%         
% %         fprintf('\nAfter %0.2f ms: \n',KernData.time{i,j});
% %         temps = temperature{i,j}(~isnan(temperature{i,j}(:)));
% %         sortedtemps = sort(temps, 'descend');
% %         PeakTemp = mean(sortedtemps(1:peakNum))
% %         MeanTemp = Tmean{i,j}
% 
%         axesCol = leftedge{i,j} + 16;
%         lineWidth = 6;
%                 
%         %Plot kernel and edge detection:
%         intensitycolormaprange = [0, 25];
%         figure
%         
%         subplot(2,4,1)
%         imshow(intREF{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange);
%         ylabel('Reflected');
%         subplot(2,4,2)
%         imshow(edgesREF{i,j}*1e6+intREF{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange);
%         
%         subplot(2,4,5)
%         imshow(intLOS{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange);
%         ylabel('Line of sight');
%         subplot(2,4,6)
%         imshow(edgesLOS{i,j}*1e6+intLOS{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange);
%         
%         %kernel edge locations
%         boundLOS = bwboundaries(binLOS{i,j});
%         boundREF = bwboundaries(binREF{i,j});
%         
%         %Organize edges into a matrix
%         edgeLOS = zeros(sizeLOS);
%         for m = 1:size(boundLOS,1)
%             for n = 1:size(boundLOS{m},1)
%                 edgeLOS(boundLOS{m}(n,1), boundLOS{m}(n,2)) = 1;
%             end
%         end
%         edgesLOS{i,j} = edgeLOS; %store in a cell structure
%         
%         edgeREF = zeros(sizeREF);
%         for m = 1:size(boundREF,1)
%             for n = 1:size(boundREF{m},1)
%                 edgeREF(boundREF{m}(n,1), boundREF{m}(n,2)) = 1;
%             end
%         end
%         edgesREF{i,j} = edgeREF;
%         
%         subplot(2,4,3)
%         imshow(edgesREF{i,j}*1e6+intREF{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange);
%         subplot(2,4,4)
%         imshow(edgesREF{i,j}*1e6+intREF{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange); hold on
%         plot([axesCol, axesCol], [refbot{i,j}(axesCol)+0.5, reftop{i,j}(axesCol)-0.5], 'r', 'LineWidth', lineWidth);
%         
%         subplot(2,4,7)
%         imshow(edgesLOS{i,j}*1e6+intLOS{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange);
%         subplot(2,4,8)
%         imshow(edgesLOS{i,j}*1e6+intLOS{i,j},'Colormap', parula, 'DisplayRange', intensitycolormaprange); hold on
%         plot([axesCol, axesCol], [losbot{i,j}(axesCol)+0.5, lostop{i,j}(axesCol)-0.5], 'r', 'LineWidth', lineWidth);
%         
%         
% %         %Plot temperature:
% %         % temperaturecolormaprange = [min(temperature{eventnum,framenum}(:)), max(temperature{eventnum,framenum}(:))];
% %         temperaturecolormaprange = [min(Temp(:)), max(Temp(:))];
% %         figure
% %         imshow(temperature{i,j}, 'Colormap', hot, 'DisplayRange', temperaturecolormaprange);
% %         title('Temperature (K)');
% %         colorbar
        
        fprintf('\nThat was Tasty!\n\n\n');
        
    catch Error_Overall %overall try-catch
        fprintf(strcat('Analysis failed on:', dataName{RedBaloons}, ' with the following error:\n'));
        LastError = Error_Overall
        fprintf('Comment out the overall try-catch-end structure to get error line numbers...');
    end %overall try-catch
    clearvars('-except',keepVars{:}); %Clear all loop variables
    waitbar(RedBaloons/length(dataName))
end

close(h)
toc
fprintf('All Done!\n\n');