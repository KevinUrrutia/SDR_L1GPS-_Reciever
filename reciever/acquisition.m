clear; clc; close all;

%----Setup
smpl_factor = 7;
dur = 0.5; %time interval of data to load
fsamp = 40e6 / smpl_factor;
Tc = 1e-3; %code period [s]
% fIF = 1610476.19047612; %intermediate frequency 
fIF = 0;
fchip = 1.023e6; %chipping rate [Hz]
del_chip = 3/217; % sampling interval in chips 
delt = del_chip * Tc; % sampling interval in seconds
fdmin = -6000; %minimum Doppler Frequency to search
fdmax = 6000; %maximum Doppler Frequency to search
deltafd = 100; %doppler resolution to search
Nc = 1023; % number of chips

%---Generate PRN code bank
PRN = createPRNBank();

%---Load Data
[I , Q] = IQ_parsing('niData01head.bin');
X = I + 1j*Q;
% [pxx, fd] = periodogram(X, [], [], fsamp, "centered");
% plot(fd/pi,10*log10(pxx));



% load('mystery_data_file.mat');
% X = Y;

%---load over accumulation period
Ns = floor(fsamp * Tc); %you get the number of samples in one code period
Xp = X(1: Ns); 

%--oversamples the code 
C = cell([1, 34]);
for i = 1:34
    C{i} = resample(PRN{i}, Ns, Nc); %interpolate the PRN code to the chipping rate
end

%---plot setup
fdVec = fdmin:deltafd:fdmax;
startTimeVec = 0:length(C{1})-1;
[Dopplerm, timem] = meshgrid(fdVec, startTimeVec/fsamp);

%--calculate N0 (noise denisity), can done by doing the autocorrelation with an unusused C/A code
% [N0, SNk2] = calcN0(fdmax, fdmin, deltafd, Ns, fsamp, C{34}, Xp);
N0 = 1.291271654429430e4;
display("Calculated N0: " + N0);


%--PRN search
SK2 = acq_svid(fdmax, fdmin, Ns, deltafd, fsamp, C, Xp, fIF);

%--Plot PRN
[ts_hat, fd_hat, C_N0] = plotPRN(SK2, Tc, fsamp, N0, startTimeVec, fdVec, Dopplerm, timem);

