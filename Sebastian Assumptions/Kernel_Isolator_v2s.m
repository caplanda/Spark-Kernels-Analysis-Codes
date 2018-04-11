%% Spark Kernels Temperature Analysis Code
clear; clc; close all; format compact; format shortg;

%% User Defined Values

testDir = 'E:\Kernel IR Data\2017_12_19\Okhovat Assumptions';

camcalFile = 'CalLookup_INT0.08208_20171130.txt';
TestSheetFile = 'TestMatrix.txt';
PressureCol = 3; %Test Sheet Pressure Column
VoltageCol = 4; %Test Sheet Voltage Column

nBG = 2; %number of background frames for each event
PlasmaThresh = 3500; %mean first frame intensity threshold for plasma existence
avenum = 500; %number of frames in the moving average
stdnum = .1; %theshold number of standard deviations from the moving average



Notes = {''};

%% Load files

datadir = [testDir,'\Datasets'];
CamCal = dlmread([testDir, '\', camcalFile]);
TestSheet = dlmread([testDir, '\', TestSheetFile]);

slashes = regexp(testDir,'\');
testdate = testDir(slashes(end)+1:end);

%Curve fit calibration data, using linear fit:
fit = fit(CamCal(:,1), CamCal(:,2)*10000, 'poly1'); % (*10000 to convert from W/cm^2-sr to W/m^2-sr)
calfit = @(x) fit.p1*x + fit.p2;


%% Choose what to process:

% % %Choose what to process using GUI:
% [dataName, dataPath] = uigetfile('*.mat',...
%     'Select the .mat file containing the video',...
%     [datadir, '\MatData'], 'MultiSelect', 'on');
% dataName = cellstr(dataName);

%%%%%% OR %%%%%%%

% % Process all files in the folder
dataName = cellstr(ls(datadir));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Chew through all files in folder:

h = waitbar(0,'Chewing through folder contents...');
RedBaloons = 99; keepVars = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
for RedBaloons = 1:length(dataName)
%     try %overall try-catch
        Pressure = TestSheet(RedBaloons,PressureCol); %Number of atmospheres
        Voltage = TestSheet(RedBaloons,VoltageCol);
        dataMatFile = matfile([datadir,'\' dataName{RedBaloons}]); %make data vars accessible
        DP = dataName{RedBaloons}(1:regexp(dataName{RedBaloons},'.mat')-1);
        fprintf(['Oooh, ',DP,'!  Om nom nom... \n']);
        
        %% Folder Prep
        
                IntImgDir = [testDir,'\Intensity_Images'];
                if exist(IntImgDir) ~= 7
                    mkdir(IntImgDir);
                end
        
                DPdir = [IntImgDir,'\',DP];
                if exist(DPdir) ~= 7
                    mkdir(DPdir);
                end
        
                figdir = [DPdir,'\figdir'];
                if exist(figdir) ~= 7
                    mkdir(figdir);
                end
        
        IntDataDir = [testDir,'\Intensity_Data'];
        if exist(IntDataDir) ~= 7
            mkdir(IntDataDir);
        end
        
        %% Find all frames with kernels
        %         meanCounts = mean(mean(dataMatFile.Data,1),2); %Time trace of the mean frame intensity
        %         meanCounts = reshape(meanCounts,[length(meanCounts),1]);
        %
        %         movingaverage = movmean(meanCounts,avenum);
        %         stddev = std(meanCounts);
        %         cutoff = movingaverage + stdnum*stddev;
        %
        %         kerndex = find(meanCounts>cutoff); %array of indicies for frames that have a kernel
        %
        %
        %
        %         %% Create spark event cell structure of indicies where kernels exist
        %
        %         %Separate out individual spark events (sort the array of kernel indicies
        %         %into a cell structure where each row contains the indicies for however
        %         %many frames exist for each event):
        %         eventcountfull = 1; %event number
        %         eventframecount = 0; %frame number for each event
        %         for i = 1:length(kerndex)-1
        %             eventframecount = eventframecount + 1;
        %             eventdexfull(eventcountfull,eventframecount) = {kerndex(i)};
        %             if kerndex(i+1)-kerndex(i) > 1
        %                 eventframecount = 0;
        %                 eventcountfull = eventcountfull + 1;
        %             end
        %         end
        %         eventdexfull(eventcountfull,eventframecount+1) = {kerndex(end)};
        %         %(the loop doesn't catch the last frame)
        %
        %         %Determine the spark frequency
        %         firstdex = cell2mat(eventdexfull(:,1));
        %         framesbetween = diff(firstdex);
        %         SparkRate = 1/ (mean(framesbetween) / dataMatFile.fps);
        %
        %
        %         %cell array of indicies for selected # of pre-event background frames):
        %         for i = 1:eventcountfull
        %             Tic = nBG+1;
        %             for j = 1:nBG
        %                 Tic = Tic - 1;
        %                 BGdexfull{i,j} = eventdexfull{i,1} - Tic;
        %             end
        %         end
        %
        %         %Remove all events that do not contain plasma on the first frame
        %         %This is because plasma in the first frame is the only way to tell timing
        %         %without an external sync source
        %         eventcount = 0;
        %         eventframes = 0;
        %         for i = 1:eventcountfull
        %             firstframemeancounts = mean2(dataMatFile.Data(:,:,eventdexfull{i,1}));
        %             eventframeswithplasma = sum(~cellfun(@isempty,eventdexfull(i,:)),2);
        %             if firstframemeancounts >= PlasmaThresh && eventframeswithplasma > 1 && eventdexfull{i,1} > nBG
        %                 eventcount = eventcount + 1;
        %                 eventframes(eventcount,1) = eventframeswithplasma - 1;
        %                 eventdex(eventcount,1:eventframeswithplasma-1) = eventdexfull(i,2:eventframeswithplasma);
        %                 BGdex(eventcount,:) = BGdexfull(i,:);
        %             end
        %         end
        
        %% Locating events using Sebastian's method:
        
        start = 1;
        finnish = size(dataMatFile.Data,3);
        tic1 = 0;
        for j = start:finnish-1
            m = max(max(dataMatFile.Data(:,:,j+1))) - max(max(dataMatFile.Data(:,:,j)));
            if m>3000 && m<=6000
                tic1 = tic1 + 1;
                plasdex(tic1) = j;
            elseif m>6000
                tic1 = tic1+1;
                plasdex(tic1) = j+1;
            end
        end
        
        tic2=0;
        tic3=0;
        for k=2:length(plasdex)                                                            %Notes the first frame of each spark?
            tic2=tic2+1;
            diff(tic2)=plasdex(k)-plasdex(k-1);
            if diff(tic2)>1
                tic3=tic3+1;
                eventfirstdex(tic3)=plasdex(k-1);
            end
            
            if k==length(plasdex)
                tic3=tic3+1;
                eventfirstdex(tic3)=plasdex(k);
            end
        end
        eventcount = tic3;
        eventframes = 5*ones(1,eventcount);
        
        % Addapt sebastian's location method to my background subtraction method:
        
        tic4 = 0;
        for i = 1:eventcount
            for j = 1:eventframes(i)
                eventdex{i,j} = eventfirstdex(i) + mod(tic4,5);
                tic4 = tic4+1;
            end
        end
        
        tic5 = 0;
        for i = 1:eventcount
            for j = 1:eventframes(i)
                BGdex{i,j} = eventfirstdex(i) -1 -mod(tic4,2);
                tic5 = tic5+1;
            end
        end
        
        %% Create cell arrays of useful data for each spark event
        
        %Make a cell structure with the raw data based on the indicies structure:
        for i = 1:eventcount
            for j = 1:eventframes(i)
                KernCounts{i,j} = dataMatFile.Data(:,:,eventdex{i,j});
            end
        end
        
        %Build a background data cell array:
        for i = 1:eventcount
            for j = 1:nBG
                BGcounts{i,j} = dataMatFile.Data(:,:,BGdex{i,j});
            end
        end
        
        %% Convert to intensity and background subtract
        
        %Mean spatial BG intensities for each event:
        for i = 1:eventcount
            %temporarily make a 3D numeric matrix of the data:
            for j = 1:nBG
                tempvar(:,:,j) = BGcounts{i,j};
            end
            %Average the BG frames into a single image for each event, and store in
            %a cell structure:
            BGmeanCounts{i,1} = mean(tempvar,3);
        end
        clear tempvar
        
        %Convert from raw photon counts to intensity, and background subtract:
        for i = 1:eventcount
            BGint{i,1} = calfit(BGmeanCounts{i,1});
            for j = 1:eventframes(i)
                KernInt{i,j} = calfit(double(KernCounts{i,j})) - BGint{i,1};
            end
        end
        
        %%Time after Plasma
        for i = 1:eventcount
            for j = 1:eventframes(i)
                time{i,j} = j * (1/dataMatFile.fps)*1000;
                framenum{i,j} = j;
            end
        end
        
        %% Plot the kernels and save images
        
%                 % Kernel existance plot:
%                 fig1 = figure;
%                 plot(meanCounts, 'b'); hold on
%                 plot(cutoff,'r--');
%                 plot([1,length(meanCounts)], [PlasmaThresh PlasmaThresh], 'g:');
%                 title('Kernel Existance');
%                 xlabel('time (Frame number)');
%                 ylabel('Mean frame intensity');
%                 legend('Mean frame intensity', 'Kernel Threshold', 'Plasma Threshold');
%                 saveas(fig1, [figdir,'\Threshold_Plot.fig']);
%                 saveas(fig1, [DPdir,'\Threshold_Plot.tif']);
        
                % Intensity images:
                ncols = 5;
                colorange = [0, 5];
                for i = 1:eventcount
                    nrows = ceil(eventframes(i)/ncols);
                    fig = figure('Position', [25, 50, 1600, 300*nrows]);
                    for j = 1:eventframes(i)
                        subplot(nrows,ncols,j)
                        imshow(KernInt{i,j},'Colormap', parula, 'DisplayRange', colorange);
                    end
                    eventname = sprintf('Event %d',i);
                    title([testdate, ', ', num2str(Pressure), ' atm, ',DP,', ' eventname]);
                    saveas(fig, [figdir,'\',eventname,'.fig']);
                    saveas(fig, [DPdir,'\',eventname,'.png']);
                end
                close all
        
        %% Save data to file
        KernData.InputIsolator.CamCal = camcalFile;
        KernData.InputIsolator.nBG = nBG;
        KernData.InputIsolator.PlasmaThresh = PlasmaThresh;
        KernData.InputIsolator.aveNum = avenum;
        KernData.InputIsolator.stdNum = stdnum;
        KernData.TestDate = testdate;
        KernData.Voltage = Voltage;
        KernData.Pressure = Pressure;
        KernData.Intensity = KernInt;
        KernData.EventCount = eventcount;
        KernData.EventFrames = eventframes;
        KernData.DP = DP;
        KernData.fps = dataMatFile.fps;
        KernData.time = time;
        KernData.framenum = framenum;
        KernData.Notes = Notes;
%         KernData.SparkRate = SparkRate;
        
        save([IntDataDir,'\',DP,'_Int.mat'], 'KernData');
        
        fprintf('\nThat was Tasty!\n\n\n');
%     catch Error_Overall %overall try-catch
%         fprintf(strcat('Analysis failed on:', dataName{RedBaloons}, ' with the following error:\n'));
%         Error_Overall
%         fprintf('Comment out the overall try,catch,end structure to get error line numbers...');
%     end %overall try-catch
    clearvars('-except',keepVars{:}); %Clear all loop variables
    waitbar(RedBaloons/length(dataName))
    
end

close(h)
toc
%%
fprintf('All Done!\n\n');
