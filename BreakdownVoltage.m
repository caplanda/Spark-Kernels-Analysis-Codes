%% Breakdown Voltage
%Based on Paschen's law, found in "The Incluence of Spark Discharge
%Characteristics on Minimum Ignition Energy in Flowing Gases" by Ballal and
%Lefebvre

clear; clc; close all;

%Constants for air:
A = 2.8e5;
B = 5.5;

P = linspace(0.01,1,100); %atm
d = 0.2; %cm

V = A*(P*d) ./ (B - log(P*d));

plot(P,V);

% The units have not been checked! Looks like the equation might not be
% compete either.