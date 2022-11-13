clear; clc; close all;

%----Setup
dur = 0.5; %time interval of data to load
fsamp = 5e6; %IQ sampling frequency
Tc = 1e-3; %code period [s]
fchip = 1.023e6; %chipping rate [Hz]
del_chip = 3/217; % sampling interval in chips 
delt = del_chip * Tc; % sampling interval in seconds
fdmin = -6000; %minimum Doppler Frequency to search
fdmax = 6000; %maximum Doppler Frequency to search
deltafd = 1000; %doppler resolution to search

%---Generate PRN code bank
PRN = createPRNBank();

%---Load Data
[I , Q] = IQ_parsing('niData01head.bin');
X = I + 1j*Q;

%---load over accumulation period
Ns = fsamp * Tc; %you get the number of samples in one code period
Xp = X(1: Ns); 

%--oversamples the code 
C = cell([1, 34]);
for i = 1:34
    C{i} = resample(PRN{i}, fsamp, fchip); %interpolate the PRN code to the chipping rate
end

%---plot setup
fdVec = fdmin:deltafd:fdmax;
startTimeVec = 0:length(C{1})-1;
[Dopplerm, timem] = meshgrid(fdVec, startTimeVec/fsamp);

%--calculate N0 (noise denisity), can done by doing the autocorrelation with an unusused C/A code
[N0, SNk2] = calcN0(fdmax, fdmin, deltafd, Ns, fsamp, C{34});
display("Calculated N0: " + N0);

%---PRN Search
SK2 = acq_svid(fdmax, fdmin, Ns, deltafd, fsamp, C); 

%---Plot PRN
[C_N0, fd_hat, ts_hat] = PRN_plot(SK2, Tc, fsamp, fdVec, Dopplerm, timem);
