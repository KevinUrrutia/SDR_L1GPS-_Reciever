%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Top Sheet: Entry point for NAV TEK
%Author: Kevin Urrutia
%Release Status: Pre-Release
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% Call GUI

%% Load user File (Replace with GUI)
fname_IQ = 'L1_IF20KHz_FS18MHz.bin';

%TODO: Pull sampling rate, center frequency, time to skip from file name

%Temporary file information
num_skip = 0; %[s]
samp_f = 18e6; %[Hz]
IF = 20e3; %[Hz]

%% Initialize constraints and settings
settings = initSettings(fname_IQ, num_skip, samp_f, IF);

%% Extract IQ from file
[I, Q] = IQ_parsing(settings.IQ);
%Create complex signal from data
X = I + 1i * Q;

%% create initial visualization
probeData();




