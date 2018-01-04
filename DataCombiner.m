%% Data Plotter

clear; clc; close all; format compact; format shortg;

%% User Defined Values

dataDir = 'E:\Kernel IR Data\CompiledDataSciTech';

%Process all files in the folder
dataName = cellstr(ls(dataDir));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Chew through all files in folder:

h = waitbar(0,'Chewing through folder contents...');
k = 0; Data.Date(1) = 0;
RedBaloons = 1; LastError = 0; keepVars = 0; TempTable = table;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
load([dataDir,'\', dataName{RedBaloons}]); %Load first table
TempTable = DataTable; %Store in temporary table
for RedBaloons = 2:length(dataName)
    load([dataDir,'\', dataName{RedBaloons}]); %Load next table
    TempTable = [TempTable; DataTable]; %Append to temporary table
    
    clearvars('-except',keepVars{:}); %Clear all loop variables
    waitbar(RedBaloons/length(dataName))
end
close all
close(h)

DataStr = table2struct(TempTable); %Rename temporary table, and then save:
save([dataDir,'\DataTable_Combined'],'DataStr');

%%
fprintf('All Done!\n\n');
toc


