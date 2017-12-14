%% Spark Kernels Temperature Images
clear; clc; close all; format compact; format shortg;

%% User Defined Values

testDir = 'J:\Kernel IR Data\2017_11_20';

TempDataDir = [testDir,'\Temperature_Data'];

colorange = [0, 3000];



%% Choose what to process:

%Enter names of files to be processed
% dataName = {'DP-000006_Int.mat'};

%%%%%% OR %%%%%%%

%Process all files in the folder
dataName = cellstr(ls(TempDataDir));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Chew through all files in folder:

h = waitbar(0,'Chewing through folder contents...');
RedBaloons = 99; ErrorLog = 0; keepVars = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
for RedBaloons = 1:length(dataName)
    try %overall try-catch
        load([TempDataDir,'\', dataName{RedBaloons}]);
        fprintf(['Oooh, ',KernData.DP,'!  Om nom nom... \n']);
        fprintf('Looks like there are %d events with up to %d frames.\n',KernData.EventCount,max(KernData.EventFrames));
        
        %% Folder Prep
        
        TempImgDir = [testDir,'\Temperature_Images'];
        if exist(TempImgDir) ~= 7
            mkdir(TempImgDir);
        end
        
        DPdir = [TempImgDir,'\',KernData.DP];
        if exist(DPdir) ~= 7
            mkdir(DPdir);
        end
        
        figdir = [DPdir,'\figdir'];
        if exist(figdir) ~= 7
            mkdir(figdir);
        end
        
        %% Plot the kernels and save images
        ncols = 5;
        for i = 1:KernData.EventCount
            nrows = ceil(KernData.EventFrames(i)/ncols);
            fig = figure('Position', [25, 50, 1600, 300*nrows]);
            for j = 1:KernData.EventFrames(i)
                subplot(nrows,ncols,j)
                imshow(KernData.temperature{i,j},'Colormap', hot, 'DisplayRange', colorange);
                colorbar
            end
            eventname = sprintf('Event %d',i);
            title([KernData.TestDate, ', ', num2str(KernData.Pressure), ' atm, ',KernData.DP,', ' eventname]);
            saveas(fig, [figdir,'\',eventname,'.fig']);
            saveas(fig, [DPdir,'\',eventname,'.png']);
        end
        close all
        
        fprintf('\nThat was Tasty!\n\n\n');
    catch Error_Overall %overall try-catch
        fprintf(strcat('Analysis failed on:', dataName{RedBaloons}, ' with the following error:\n'));
        ErrorLog(RedBaloons) = Error_Overall
        fprintf('Comment out the overall try,catch,end structure to get error line numbers...');
    end %overall try-catch
    clearvars('-except',keepVars{:}); %Clear all loop variables
    waitbar(RedBaloons/length(dataName))
end
close all
close(h)
toc
fprintf('All Done!\n\n');