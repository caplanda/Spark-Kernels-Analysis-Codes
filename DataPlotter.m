clear; clc; close all;

%% Inputs:
load('J:\Kernel IR Data\2017_11_20\AllData.mat');
isolineString = 'Pressure';
xaxisString = 'Time'; %Independent
yaxisString = 'Tmean'; %Dependent


DataTable = struct2table(Data);


%% The function:

isoName = struct;

isolines = unique(extractfield(Data,isolineString));
% voltages = unique(extractfield(Data,'Voltage'));
xaxis = unique(extractfield(Data,xaxisString));

for i = 1:length(isolines)
    for x = 1:length(xaxis)
        rows = DataTable.(isolineString) == isolines(i) &...
            DataTable.(xaxisString) == xaxis(x) &...
            ~isnan(DataTable.(yaxisString)) &...
            DataTable.(yaxisString) > 0;
        isoName.(yaxisString)(x) = mean(DataTable.(yaxisString)(rows)); %Bug: overwrites the previous line
    end
end

isoName.(isolineString) = isolines;
isoName.(xaxisString) = xaxis;