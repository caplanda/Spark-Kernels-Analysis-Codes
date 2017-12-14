%% Spark Kernels Visualization Code
clear; clc; close all; format compact; format shortg;

%% File Info:

% testDir = 'C:\Users\caplanda\Kernel IR Data\2017_11_20';
testDir = 'F:\Kernel IR Data\2017_11_20';

layoutFile = 'Layout.mat';
camcalFile = 'CalLookup_INT0.08208_20171130.txt';

nBG = 3; %number of frames prior to the beginning of each event that will be averaged and used for background subtraction
% mirrorReflectivity = 0.85; %Percent of light reflected

testdate = '11/20/2017';


%% Folder/file prep

datadir = [testDir,'\Datasets'];
IntensityDir = [testDir,'\Intensity_Images'];
if exist(IntensityDir) ~= 7
    mkdir(IntensityDir);
end


fprintf('Loading files... \n');
load([testDir, '\', layoutFile]);
CamCal = dlmread([testDir, '\', camcalFile]);

%Curve fit calibration data, using linear fit:
fit = fit(CamCal(:,1), CamCal(:,2)*10000, 'poly1'); % (10000 to convert from W/cm^2-sr to W/m^2-sr)
calfit = @(x) fit.p1*x + fit.p2;

%% Choose what to process:

%Enter names of tiff folders to be converted
% vidName = {'DP3_Ignition'};

%%%%%% OR %%%%%%%

%Process all folders in the tiffs folder
dataName = cellstr(ls(datadir));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Load files

% [dataName, dataPath] = uigetfile('*.mat',...
%     'Select the .mat file containing the video',...
%     [datadir, '\MatData']);

for k = 1:length(dataName)
    
    %Load reflected and line of sight data from the .mat file:
    dataMatFile = matfile([datadir,'\', dataName{k}]); %make data vars accessible
    refdata = dataMatFile.Data(refind(1):refind(2), refind(3):refind(4),:);
    losdata = dataMatFile.Data(losind(1):losind(2), losind(3):losind(4),:);
    
    %% Make folders for saving (if they don't exist already)
    DP = dataName{k}(1:regexp(dataName{k},'.mat')-1);
    
    DPdir = [IntensityDir,'\',DP];
    if exist(DPdir) ~= 7
        mkdir(DPdir);
    end
    
    figdir = [DPdir,'\figdir'];
    if exist(figdir) ~= 7
        mkdir(figdir);
    end
    
    %% Find all frames with kernels
    losmean = mean(mean(losdata,1),2); %Time trace of the mean frame intensity
    losmean = reshape(losmean,[length(losmean),1]);
    
    avenum = 500; %number of frames in the moving average
    stdnum = .15; %theshold number of standard deviations from the moving average
    movingaverage = movmean(losmean,avenum);
    stddev = std(losmean);
    cutoff = movingaverage + stdnum*stddev;
    
    kerndex = find(losmean>cutoff); %array of indicies for frames that have a kernel
    
    fig1 = figure;
    plot(losmean, 'b'); hold on
    plot(cutoff,'r');
    title('Kernel Existance');
    xlabel('time (Frame number)');
    ylabel('Mean frame intensity');
    legend('Mean frame intensity', 'Cutoff');
    
    clear movingaverage stddev cutoff losmean
    
    %% Create spark event cell structure of indicies where kernels exist
    
    %Separate out individual spark events (sort the array of kernel indicies
    %into a cell structure where each row contains the indicies for however
    %many frames exist for each event):
    eventcount = 1; %event number
    eventframecount = 0; %frame number for each event
    for i = 1:length(kerndex)-1
        eventframecount = eventframecount + 1;
        eventdex(eventcount,eventframecount) = {kerndex(i)};
        if kerndex(i+1)-kerndex(i) > 1
            eventframecount = 0;
            eventcount = eventcount + 1;
        end
    end
    eventdex(eventcount,eventframecount+1) = {kerndex(end)};
    %(the loop doesn't catch the last frame)
    
    %cell array of indicies for selected # of pre-event background frames):
    for i = 1:eventcount
        tic = nBG+1;
        for j = 1:nBG
            tic = tic - 1;
            BGdex{i,j} = eventdex{i,1} - tic;
        end
    end
    
    clear eventframecount tic i j
    
    %% Create cell arrays of useful data for each spark event
    
    %Make a cell structure with the raw data based on the indicies structure:
    for i = 1:eventcount
        for j = 1:sum(~cellfun(@isempty,eventdex(i,:)),2) %(# kernel frames in each event)
            losRaw{i,j} = losdata(:,:,eventdex{i,j});
            refRaw{i,j} = refdata(:,:,eventdex{i,j});
        end
    end
    
    %Build a background data cell array:
    for i = 1:eventcount
        for j = 1:sum(~cellfun(@isempty,BGdex(i,:)),2) %(# kernel frames in each event)
            losBGRaw{i,j} = losdata(:,:,BGdex{i,j});
            refBGRaw{i,j} = refdata(:,:,BGdex{i,j});
        end
    end
    %Bug: if first event is too close to the beginning of the data, then
    %background frames prior to event will be out of bounds.
    
    %Find frames with plasma: (maybe)
    for i = 1:eventcount
        for j = 1:sum(~cellfun(@isempty,BGdex(i,:)),2) %(# kernel frames in each event)
            maxval{i,j} = max(losRaw{i,j}(:));
        end
    end
    %try eliminating frames where the max value is above some threshold. Can
    %also say only if above threshold in both the los and ref view. Or only if
    %above threshold for max value as well as average value. Or both.
    
    clear refdata losdata i j
    
    %% Convert to intensity and background subtract
    
    %Mean spatial BG intensities for each event:
    for i = 1:eventcount
        %temporarily make a 3D numeric matrix of the data:
        for j = 1:sum(~cellfun(@isempty,losBGRaw(i,:)),2) %(# kernel frames in each event)
            lostemporary(:,:,j) = losBGRaw{i,j};
            reftemporary(:,:,j) = refBGRaw{i,j};
        end
        %Average the BG frames into a single image for each event, and store in
        %a cell structure:
        losBGRawMean{i,1} = mean(lostemporary,3);
        refBGRawMean{i,1} = mean(reftemporary,3);
    end
    
    %Convert from raw photon counts to intensity, and background subtract:
    for i = 1:eventcount
        %     losBG{i,1} =  interp1( CamCal(:,1), CamCal(:,2), double(losBGRawMean{i,1}) );
        %     refBG{i,1} =  interp1( CamCal(:,1), CamCal(:,2), double(refBGRawMean{i,1}) );
        losBG{i,1} = calfit(losBGRawMean{i,1});
        refBG{i,1} = calfit(refBGRawMean{i,1});
        for j = 1:sum(~cellfun(@isempty,losRaw(i,:)),2) %(# kernel frames in each event)
            %         losInt{i,j} = interp1( CamCal(:,1), CamCal(:,2), double(losRaw{i,j}) );
            losInt{i,j} = calfit(double(losRaw{i,j}));
            losInt{i,j} = losInt{i,j} - losBG{i,1};
            
            %         refInt{i,j} = interp1( CamCal(:,1), CamCal(:,2), double(refRaw{i,j}) );
            refInt{i,j} = calfit(double(refRaw{i,j}));
            refInt{i,j} = refInt{i,j} - refBG{i,1};
        end
    end
    
    % clear i j lostemporary reftemporary losBGRawMean refBGRawMean losBG refBG CamCal losRaw refRaw losBGRaw refBGRaw
    
    
    %% Plot the kernels
    ncols = 5;
    colorange = [0, 100];
    for i = 1:eventcount
        eventframes = sum(~cellfun(@isempty,losInt(i,:)),2); %(# kernel frames in each event)
        nrows = ceil(eventframes/ncols);
        fig = figure('Position', [25, 50, 1600, 300*nrows]);
        for j = 1:eventframes
            subplot(nrows,ncols,j)
            imshow(losInt{i,j},'Colormap', parula, 'DisplayRange', colorange);
        end
        title([testdate, ', ', dataName{k}, ', Event ', num2str(i)]);
        eventname = sprintf('Event_%d',i);
        saveas(fig, [figdir,'\',eventname,'.fig']);
        saveas(fig, [DPdir,'\',eventname,'.tif']);
    end
    close all
    clear kerndex
end
