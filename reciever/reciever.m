clear; clc; close all;

%---Setup
fc = 1575.42e6; %Center Frequency 
smpl_factor = 7;
fIF = 1610476.19047612; %intermediate frequency 
Tc = 1 / (1.023e6); %code period [s]
Nc = 1023; % number of chip
fsamp = 40e6 / smpl_factor; %sampling rate
fchip = 1.023e6; %chipping rate [Hz]
fdmin = -6000; %minimum Doppler Frequency to search
fdmax = 6000; %maximum Doppler Frequency to search
deltafd = 100; %doppler resolution to search
T_sub = 1e-3; %sub-accumulation period

SVn = [1 2 3 5 6 7 10 12 13 29 34]; %satellites to search

N = floor(fsamp * (T_sub * 100)); %determine how much data to load in total, in this case 60s

%---load in data
load('mystery_data_file.mat');
X = Y(1:N);

%--Plot power spectrum
% [pxx, fd] = periodogram(X, [], [], fsamp, 'centered');
% figure;
% plot(fd, 10*log10(pxx));


Ns = floor(fsamp * T_sub); %determine how much data to load for chipping period

%---Generate PRN code Bank
PRN = createPRNBank();

%---Oversample PRN code
C = cell([1, 34]);
for i = 1:34
    C{i} = resample(PRN{i}, Ns, Nc);
end

%---plot setup
fdVec = fdmin:deltafd:fdmax;
startTimeVec = 0:length(C{1}) - 1;
[Dopplerm, timem] = meshgrid(fdVec, startTimeVec/fsamp);

%---Acquisition
while(length(X)  - Ns >= 0) 
    Xp = X(1 : Ns);
    acquisition_f(fdVec, Xp, C, Dopplerm, timem, Ns, fsamp, fIF, SVn, Tc, startTimeVec);
    X = X(Ns : end);
end

