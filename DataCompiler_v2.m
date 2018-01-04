%% Data Compiler
%Compiles all the data for a particular test directory into a single table
%The data is first compiled into a spreadsheet style structure, and then
%converted to a table and saved. 

clear; clc; close all; format compact; format shortg;

%% User Defined Values

testDir = 'E:\Kernel IR Data\2017_11_20';

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
k = 0; Data = struct;
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
                    k = k + 1; %Kernel frame counter - each frame of each event gets its own row in the struct spreadsheet
                    Data(k).Date = KernData.TestDate;
                    Data(k).DP = KernData.DP;
                    Data(k).fps = KernData.fps;
                    Data(k).Pressure = KernData.Pressure;
                    Data(k).Voltage = KernData.Voltage;
                    Data(k).SparkRate = KernData.SparkRate;
                    Data(k).Time = KernData.time{i,j};
                    Data(k).Intensity = KernData.Intensity{i,j};
                    Data(k).Temperature = KernData.temperature{i,j};
                    Data(k).Tmean = KernData.Tmean{i,j};
                    Data(k).Tpeak = max(KernData.temperature{i,j}(:));
                    Data(k).Volume = KernData.Volume{i,j};
                    Data(k).binLOS = KernData.binLOS{i,j};
                    Data(k).binREF = KernData.binREF{i,j};
                    Data(k).intLOS = KernData.intLOS{i,j};
                    Data(k).intREF = KernData.intREF{i,j};
                    Data(k).PixSizeLOS = KernData.PixSizeLOS;
                    Data(k).PixSizeREF = KernData.PixSizeREF;
                    Data(k).IgTip = KernData.IgTip;
                    
                    %Calculate and store mean kernel intensities:
                    intLOSkern = KernData.binLOS{i,j} .* KernData.intLOS{i,j};
                    Data(k).intLOSavg = sum(intLOSkern(:)) / length(find(intLOSkern));
                    intREFkern = KernData.binREF{i,j} .* KernData.intREF{i,j};
                    Data(k).intREFavg = sum(intREFkern(:)) / length(find(intREFkern));
                    
                    Data(k).InputIsolator = KernData.InputIsolator;
                    Data(k).InputTemp = KernData.InputTemp;
                    Data(k).Notes = KernData.Notes;
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

DataTable = struct2table(Data); %Convert the struct to a table to make organizing and plotting easier

slashes = regexp(testDir,'\');
testDate = testDir(slashes(end)+1:end);
save([testDir,'\',testDate,'_DataTable'],'DataTable');



%%
fprintf('All Done!\n\n');
toc


