% -------------------------------------------------------------------
%                          FilterApply.m                               
% -------------------------------------------------------------------
% Written by Daniel Caplan                                         
% For use by the Blunck Research Group at Oregon State University  
% Contact caplanda@oregonstate.edu to report issues.               
% Created 3/29/2017                                           
% -------------------------------------------------------------------
% Input:                                                           
% (specIn, plotSwitch), gui response profile selection                       
%      specIn:     wavelengths in col 1, intensities in col 2      
%      plotSwitch: enter 1 to the plot input, filters, and output
%      Profiles:    wavelengths in col 1, transmisivity in col 2    
% -------------------------------------------------------------------
% Output:                                                          
% [specOut], printed list of Profiles applied, and plot if requested                                     
%      specOut:    wavelengths in col 1, modified intensity in col 2
% -------------------------------------------------------------------
% To make finding the profiles more convenient, the default search 
% directory can be set changed. Simply modify the defaultDirectory
% variable, which is the first variable in this function. 
% -------------------------------------------------------------------
% Updates:
% -------------------------------------------------------------------

function [ specOut ] = FilterApply( specIn, plotSwitch )

%% Select response profile files:

defaultDirectory = '\\depot.engr.oregonstate.edu\students\c\caplanda\Windows.Documents\My Documents\CIRE Lab\Caplan\Matlab Codes\IR Response Profiles';

[fName, fPath] = uigetfile('*.txt', 'Select response profiles (hold ctrl to select multiple)',...
    defaultDirectory, 'MultiSelect', 'on');
fName = cellstr(fName); %make sure that fName is a cell, allowing a single profile to be used.

%% Apply filters:

fprintf('Profiles Applied: \n');
specOut = specIn;
for j = 1:length(fName)
    profile = dlmread(strcat(fPath,fName{j}));
    for i = 1:length(specIn(:,1))
        if specIn(i,1) >= min(profile(:,1)) && specIn(i,1) <= max(profile(:,1))
            trans(i) = fixpt_interp1(profile(:,1), profile(:,2), specIn(i,1), ...
                float('double'), 1, float('double'), 1, 'Nearest');
        else
            trans(i) = 0;
        end
    end
    specOut(:,2) = specOut(:,2) .* trans';
    fprintf([fName{j},'\n']);
end

%% Plot (if requested):

if plotSwitch == 1
    figure
    
    subplot(3,1,1)
    plot(specIn(:,1),specIn(:,2));
    xlim([0 6]);
    title('Input Spectrum');
    grid on
    
    subplot(3,1,2)
    for j = 1:length(fName)
        profile = dlmread(strcat(fPath,fName{j}));
        plot(profile(:,1), profile(:,2)); hold on
    end
    xlim([0 6]);
    ylim([0 1]);
    title('Response profiles');
    legend(fName,'Location','west');
    grid on
    
    subplot(3,1,3)
    plot(specOut(:,1),specOut(:,2));
    xlim([0 6]);
    title('Output Spectrum');
    grid on
end

end