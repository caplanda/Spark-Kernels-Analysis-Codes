%% Adds information to existing mat files

clear; clc; close all;

%% User Defined Values

Directory = 'E:\Kernel IR Data\2017_11_20\Temperature_Data';

%Load the info to append:
load('E:\Kernel IR Data\2017_11_20\Layout.mat');

%% Choose what to process:

%Enter names of files to be processed
% dataName = {'DP-000006_Int.mat'};

%%%%%% OR %%%%%%%

%Process all files in the folder
dataName = cellstr(ls(Directory));
dataName = dataName(3:end); %get rid of stupid dots added by ls

%% Chew through all files in folder:

h = waitbar(0,'Chewing through folder contents...');
k = 0; Data = struct;
RedBaloons = 99; LastError = 0; keepVars = 0;
keepVars = who; %Non-loop variables protected from being cleared on each iteration
tic
for RedBaloons = 1:length(dataName)
    load([Directory,'\', dataName{RedBaloons}]);
    
    %Append what you want:
    KernData.IgTip = lostip;
    
    save([Directory,'\', dataName{RedBaloons}], 'KernData');
    clearvars('-except',keepVars{:}); %Clear all loop variables
    waitbar(RedBaloons/length(dataName))
end
close all
close(h)

%%
fprintf('All Done!\n\n');
toc


