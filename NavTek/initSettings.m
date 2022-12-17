function settings = initSettings(fname_IQ, num_skip, samp_f, IF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Initializes and saves user settings. Settings are edited inside the
%%function, these are masked from the user. The only settting they have
%%acccess to are the IQ file
%%Inputs: fnameIQ -> IQ file location
%         num_skips -> seconds to skip in the file [s]
%         samp_f -> samping frequency [Hz]
%         IF -> intermediate frequency [Hz]
%%Outputs: conf_settings - > Receiever Settings (struct)
%
%
%%Author: Kevin Urrutia
%%Release Status: Pre-Release
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Processing Settings
%IQ File location
settings.IQ = fname_IQ;

%Seconds to process 
settings.msToProcess = 37000; %[milli-seconds]

%Number fo channels to be used for signal processing, num of sats
settings.numChannels = 12; 

%save sampling rate
settings.samp_f = samp_f;

%Number of seconds to skip
settings.skipNumberBytes = num_skip * settings.samp_f * 4; %there are 4 bytes per sample

%Data Type used to store one sample
settings.dataType = 'schar';

%Intermediate, sampling and code Frequencies
settings.IF = IF;
settings.codeFreqBasis = 1.023e6; %[Hz]

%number of code chips in a period
settings.codeLength = 1023;

%% Aquisition Settings

%List of satellites to look for. Exclude satellites to speed up aqusition
settings.SVID = 1:32; %[PRN Numbers]

%Band around IF to search for satellites. (Max Doppler)
settings.maxDopp = 7000; %[Hz]

%Threshold for signal presence
settings.acqThresh = 3.5; 

%Freqeuncy search step
settings.acqSearchStep = 100;

%% Constants
settings.c = 299792458; %Speed of light, [m/s]

end