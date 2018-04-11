%% Spark Kernals Power Trace Analysis %%

%TO RUN 
%SELECT CORRECT FILE FOR RESISTENCE
%correct directory to your computer and the desired voltage data
%folder and enter desired filename in csvwrite on line 106

clear; clc; clf; close all;



%% Load Files %%

%change directory on your computer
dataDir = 'C:\Users\Danny\Desktop\SciTech Download\Spark-Kernels-Analysis-Codes\Spark kernel power\2017_12_22 Voltage Data';
testDir = 'C:\Users\Danny\Desktop\SciTech Download\Spark-Kernels-Analysis-Codes\Spark kernel power';
TestSheetFile = '2017_12_22 test sheet.txt';

dataName = cellstr(ls(dataDir));
dataName = dataName(3:end); %remove ls dots
%% Process files in the folder %%

%FIND TRUE RESISTENCE VALUES
R = 2; %ohm
V_supply = 15; %V
keepVars_1 = who; %protect non-loop variables from being cleared

%call pressures
TestSheet = dlmread([testDir, '\', TestSheetFile]);
PressureCol = 2; %test sheet pressure column

for count_1 = 1:length(dataName)
    %load file
    file = dataName{count_1};
    Trace = csvread([dataDir, '\', file], 40, 0); %row offset to remove garbage leftover from previous putty recording
    
    Trace(length(Trace),:) = []; %to remove the last line that putty might not have finished when test was ended
    
    %pull timestamp and voltage through thte resistor from the trace
    t = Trace(:,1);
    V = Trace(:,2);
    
    %filter the voltage
    V_filtered = medfilt1(V);
    
    %find current through the resistor
    I = V_filtered./R;
       
    
    %% Find multiple period in the trace 
    k = 1;   %row variable for the current matrix
    k_old = 1; %initial value same as k for logic
    j = 1;   %column variable for current matrix
    
    I_new = 0;
    I_old = 0;
    I_older = 0;
    %I_oldest = 0; %all four cannot begin as 0 in case the first point chosen is 0 or first row of the matrix would be empty
    
    I_period = zeros(40,30); %periods appear to be around 30 data points and there are around 25 periods per sample. This ensures sufficient room to store all the data.
    t_period = zeros(40,30);
    
    for count_2 = 1:length(Trace)
        
        if k ~= k_old 
            j = 1;     %if a new period is found the column variable is reset to the left side of the matrix
        end
        
        I_new = I(count_2);
        I_period(k,j) = I_new; 
        t_period(k,j) = t(count_2);
        
        j = j+1; %advance column
        k_old = k;  %keep track for new row detection
    
        %advance row if the end of the period is reached each period from 
        %the current trace is stored in a new row for per spark calculations
        if (I_new<=0.05) && (I_old<=0.05) && (I_older<=0.05) %with new resistence values do not equal zero
            k = k+1;    %if all 3 are less than 0.00 the period is over and progresses to the new row in the matrix
        end
        
        %update I values for period detection
        I_oldest = I_older;
        I_older = I_old;
        I_old = I_new;
        
    end
    
    %reomve zero rows from matrix in case of gaps longer than 3 0s
    t_period( all(~I_period,2), : ) = [];
    I_period( all(~I_period,2), : ) = [];
    
    %% CHECK FOR CUTOFF
    I_sum = sum(I_period,2);    %to check for cut periods
    indices = find(abs(I_sum)<1.55); %CHECK WITH DAN. If there are 3 0s in a row within a period it will split it. This removes the partial periods to avoid affecting the averages
    I_period(indices,:) = [];
    t_period(indices,:) = [];

    %% find instantaneous power
    
    P_inst = I_period.*V_supply; % Instantaneous power through ignitor
    Energy = zeros(10,1); %reset for each file
    
    %% Integrate for Total energy per period
    for count_3 = 1:size(P_inst,1)
        Energy(count_3,1) = trapz(t_period(count_3,:), P_inst(count_3,:));
    end
    l = length(Energy);
    Energy(l) = [];
    %% Mean of periods per trial linked with file name
    Pressure = TestSheet(count_1,PressureCol);
    
    E_mean(count_1,1) = dataName(count_1);
    E_mean(count_1,2) = {Pressure};
    E_mean(count_1,3) = {mean(Energy)}; %stored as string to prevent cell error. 

end
%% write to file
[n_rows,n_col] = size(E_mean);
fid = fopen('2017_12_22 E_mean Data.csv','w');
for count_4 = 1:n_rows
    fprintf(fid,'%s, %s\n',E_mean{count_4,:});
  %  fprintf(fid,'%f, %f, %f\n',c{2:end,:})
end
fclose(fid);

figure
scatter(cell2mat(E_mean(:,2)),cell2mat(E_mean(:,3))/1000, 20);
axis([0.2 1.1, 1.9, 2.6]);
xlabel('Pressure');
ylabel('Average Energy per Spark (J)');
box on