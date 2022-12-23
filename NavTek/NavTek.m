classdef NavTek

    properties
        %% Post processing settings
        %IQ file name
        fnameIQ = "";
        %Number of milliseconds to process, used 36000 + transients,
        %enough Nav subframes need to be provided
        msToProcess = 37000; %[ms]
        %Number of channels used for signal processing
        numChannels = 12;
        %Move to the starting point of the post processing.
        skipNumBytes = 0;
        %Date type to store samples
        dataType = 'schar';
        %Intermediate, sampling and code frequencies
        fIF = 20e3; %[Hz]
        fsamp = 18e6; %[Hz]
        fchip = 1.023e6; %[Hz]
        %Code Length
        codeLength = 1023;

        %% Acquisition Settings
        %satellite PRN numbers
        SVn = 1:32;   
        %minimum, maximum, resolution of doppler search
        fdmin = -7000; %[Hz]
        fdmax = 7000; %[Hz]
        delta_fd = 1000; %[Hz]
        %NonCoherent integration times after 1ms coherent integration
        NonCohTime = 20; 
        %Sampling rate for downsampling 
        resamplingThresh = 8e6; %[Hz]
        %acquisition threshold for desision rule
        acqThresh = 3.5; 
        %m sequence table
        caTable;

        %% Constants
        c = 299792458;
    end

    methods (Static)
        function CA = generateGoldSeq(fsamp, fchip, codeLength, svid, digit_flag)
            S = load('G2_taps.mat');

            n_stages = 10;

            %init registers
            G1 = ones(1, n_stages);
            G2 = ones(1, n_stages);
    
            %pre-allocate
            CA = zeros(1, 1023);

            %create sequence
            for i = 1:1023
                [g1, G1] = NavTek.shift(G1, [3, 10], 10);
                [g2, G2] = NavTek.shift(G2, [2, 3, 6, 8, 9, 10], S.G2_taps{svid});

                CA(i) = mod((g1 + g2), 2);
            end
            CA(CA == 0) = -1;

            if (digit_flag)
                %==Digitizing=================================================
                % Make index array to read CA code values
                %The length of the index depends on the sampling freqeuncy 
                %the number of samples per millisecond(b/c one CA code period is 1 millisecond)
                samplesPerCode = round(fsamp / ...
                                (fchip / codeLength));
                ts = 1 / fsamp;
                tc = 1 / fchip;
    
                codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);
                %--- Correct the last index (due to number rounding issues) -----------
                codeValueIndex(end) = 1023;
    
                %--- Make the digitized version of the C/A code -----------------------
                % The "upsampled" code is made by selecting values form the CA code
                % chip array (caCode) for the time instances of each sample.
                CA = CA(codeValueIndex);
            end
        end 


        function [out, register] = shift(register, feedback, output)
            %GPS shift register
            %input feedback -> which positions to use as feedback (1 indexed)
            %      output -> which positions are used for the output
            %output -> output of shift register
            
            %calculate the output 
            out = zeros(1, size(output, 2));
            for i = 1:length(output)
                out(i) = register(output(i));
            end
            
            if (size(out, 2) > 1)
                out = mod(sum(out), 2);
            else
                out = out(1);
            end
            
            element = zeros(1, size(feedback, 2));
            for i = 1:length(feedback)
                element(i) = register(feedback(i));
            end
            fb = mod(sum(element), 2);
            
            %shift to the right
            register = circshift(register,1);
            register(1) = fb;
            
        end
    end

    methods
        function obj = NavTek()
            prompt = "Enter full IQ file path";
            disp(prompt);
            [file, path] = uigetfile("*.bin");
            if isequal(file, 0)
                error("Enter a valid file");
            else 
                disp(['User selected ', fullfile(path,file)]);
            end
            
            obj.fnameIQ = fullfile(path, file);
        end
        
        function obj = probeData(obj)
            %% Generate and plot raw data
            fid = fopen(obj.fnameIQ, 'rb');

            if (fid > 0)
                %move to the start of the where we want to start post
                %processing
                fseek(fid, obj.skipNumBytes, 'bof');

                %find the number of samples per chipping sequence
                sampsPerCode = round(obj.fsamp / (obj.fchip / obj.codeLength));

                %read in the data 
                [data, count] = fread(fid, [1, 2 * 100 * sampsPerCode], obj.dataType);

                fclose(fid);

                if (count < 2 * 100 * sampsPerCode)
                    error("Error not enough data to load in");
                end

                %initialize plots
                timeScale = 0: 1/obj.fsamp : 5e-3; 
                f = figure;
                f.Position(3:4) = [1000, 1000];

                %% Time domain plot
                data=data(1:2:end) + 1i .* data(2:2:end);
                subplot(3, 2, 4);
                plot(1000 * timeScale(1:round(sampsPerCode/2)), ...
                real(data(1:round(sampsPerCode/2))));

                axis tight;    grid on;
                title ('Time domain plot (I)');
                xlabel('Time (ms)'); ylabel('Amplitude');

                subplot(3, 2, 3);
                plot(1000 * timeScale(1:round(sampsPerCode/2)), ...
                imag(data(1:round(sampsPerCode/2))));

                axis tight;    grid on;
                title ('Time domain plot (Q)');
                xlabel('Time (ms)'); ylabel('Amplitude');

                %% Frequency domain plot
                subplot(3,2,1:2);
                [sigspec,freqv]=pwelch(data, 32768, 2048, 32768, obj.fsamp,'twosided');
                plot(([-(freqv(length(freqv)/2:-1:1));freqv(1:length(freqv)/2)])/1e6, ...
                10*log10([sigspec(length(freqv)/2+1:end);
                sigspec(1:length(freqv)/2)]));
                axis tight;
                grid on;
                title ('Frequency domain plot');
                xlabel('Frequency (MHz)'); ylabel('Magnitude');

                %% Histogram
                subplot(3, 2, 6);
                histogram(real(data), -128:128)
                dmax = max(abs(data)) + 1;
                axis tight;     adata = axis;
                axis([-dmax dmax adata(3) adata(4)]);
                grid on;        title ('Histogram (I)');
                xlabel('Bin');  ylabel('Number in bin');
            
                subplot(3, 2, 5);
                histogram(imag(data), -128:128)
                dmax = max(abs(data)) + 1;
                axis tight;     adata = axis;
                axis([-dmax dmax adata(3) adata(4)]);
                grid on;        title ('Histogram (Q)');
                xlabel('Bin');  ylabel('Number in bin');

                uiwait(f);
            end
        end

        function acq_results = aquisition(obj, data)
            %% Initialization
            %==Variable initialization for course aquisition==============
            %number of samples per chipping sequence
            samplesPerCode = round(obj.fsamp / (obj.fchip / obj.codeLength));
            %find sampling period 
            ts = 1 / obj.fsamp;
            %Find phase points of 2ms local carrier wave (1 ms for local duplicate,
            % other 1ms for zeros padding)
            phasePoints = (0: (samplesPerCode*2-1)) * 2 * pi * ts;
            %number of frequency bins for specified search band
            numFreqBins = round(obj.fdmax * 2 / obj.delta_fd) + 1;
            %carrier frequency bins to be searched
            courseFreqBins = zeros(1, numFreqBins);

            %==Initialize aquisition results==============================
            %carrier frequencies of detected signals
            acq_results.carrFreq = zeros(1, 32);
            %C/A phases of deteceted signals
            acq_results.codePhase = zeros(1, 32);
            %Correlation peak rations of detected signals
            acq_results.peak = zeros(1, 32);

            %==Variables for fine aquisition==============================
            %fine frequency search step
            fineSearchStep = 25;
            %Number of frequency bins for fine aquisition
            numFineBins = round(obj.delta_fd / fineSearchStep) + 1;
            %Carrier frequencies of fine frequency bins
            fineFreqBins = zeros(1, numFineBins);
            %Correlation values for all find freqeucny bins
            fineResult = zeros(1, numFineBins);
            %Coherent integration of 40 codes
            sumPerCode = zeros(1, 40);
            %phase points of the local carrier wave 
            finePhasePoints = (0 : (40 * samplesPerCode-1)) * 2 * pi * ts;

            %==Input signal power for GRLT statistic calculation
            sigPower = sqrt(var(data(1:samplesPerCode)) * samplesPerCode);

            %Preform search for all listed satellite numbers....
            fprintf('(');
            for PRN = obj.SVn
                %% Course Aquisition
                %generate PRN sqequences
                caCodesTable = obj.generateGoldSeq(obj.fsamp, obj.fchip, obj.codeLength, PRN, true);
                %add zero padding samples
                caCode2ms = [caCodesTable, zeros(1, samplesPerCode)];
                %search results of all frequency bins and code shifts
                results = zeros(numFreqBins, samplesPerCode * 2);
                %perform DFT of CA codes
                caCodeFreqDom = conj(fft(caCode2ms));
                Sk = cell([1, 30]);

                %make correlation for all frequency bins
                for freqBinIndex = 1:numFreqBins
                    %generate carrier wave frequency grid
                    courseFreqBins(freqBinIndex) = obj.fIF + obj.fdmax - ...
                        obj.delta_fd * (freqBinIndex - 1);
                    
                    %generate local sine and cosine
                    sigCarr = exp(-1i * courseFreqBins(freqBinIndex) * phasePoints);

                    %Complete Correlation
                    for nonCohIndex = 1: obj.NonCohTime
                        %Take 2ms of input data to do correlation
                        signal = data((nonCohIndex - 1) * samplesPerCode + 1: (nonCohIndex + 1) * samplesPerCode);
                        
                        %remove the carrier from the signal
                        I = real(sigCarr .* signal);
                        Q = imag(sigCarr .* signal);
                        %convert baseband signal to the frequency domain
                        IQfreqDom = fft(I + 1i * Q);
                        %Multiplication in frequency domain corresponds to
                        %correlation in time domain
                        convCodeIQ = IQfreqDom .* caCodeFreqDom;
                        %perform DFT and store correlation results
                        cohResults = abs(ifft(convCodeIQ));
                        %Non Coherent Integration
                        results(freqBinIndex, :) = results(freqBinIndex, :) + cohResults;
                    end
                end
                %% Find Correlation peaks for CA
                %find correlation peaks for CA and teh carrier frequency
                [~, acqCoarseBin] = max(max(results, [], 2));
                %find the code phase for the same correlation peak
                [peakSize, codePhase] = max(max(results));
                %store the GLRT statistic
                acq_results.peakMetric(PRN) = peakSize/sigPower/obj.NonCohTime;

                %% Plot Aquistion results
                S = mesh(results);
                grid on;
                xlabel('f_d [Hz]')
                ylabel('t_s [sec]')
                title(strcat(['|S_k|^2 for PRN ' int2str(PRN)]))
                pause(1);
                clf(S); 
                
                %if result is above the threshold then there is a signal
                %present 
                %% Fine carrier frequency search
                if(acq_results.peakMetric(PRN) > obj.acqThresh)
                    %indicate PRN number of detected signal
                    fprintf('%02d ', PRN);

                    %Prepare 20ms of code, carrier and input signals
                    %CA code with 10230 chips
                    caCode = obj.generateGoldSeq(obj.fsamp, obj.fchip, obj.codeLength, PRN, false);
                    % C/A code sample index
                    codeValueIndex = floor((ts * (0 : 40*samplesPerCode -1)) / ...
                                    (1/obj.fchip));
                    % C/A code samples
                    caCode40ms = caCode(rem(codeValueIndex, obj.codeLength) + 1);
                    % Take 40cm incoming signal for fine acquisition
                    sig40cm = data(codePhase:codePhase + 40*samplesPerCode -1);

                    %search fine freqency bins
                    for fineBinIndex = 1 : numFineBins
                        %--- Correlation for each code --------------------------------
                        % Carrier frequencies of the frequency bins
                        fineFreqBins(fineBinIndex) = courseFreqBins(acqCoarseBin) + ...
                            obj.delta_fd/2 - fineSearchStep * (fineBinIndex - 1);

                        % Local carrier signal
                        sigCarr40cm = exp(-1i * fineFreqBins(fineBinIndex) * finePhasePoints);
                        % Wipe off code and carrier from incoming signals
                        basebandSig = sig40cm .* caCode40ms .* sigCarr40cm;

                        % Coherent integration for each code
                        for index = 1:40
                            sumPerCode(index) = sum(basebandSig( samplesPerCode * ...
                                (index - 1) + 1 : samplesPerCode * index ));
                        end

                        %--- Search Nav bit edge for ----------------------------------
                        % 20 cases of Nav bit edge
                        maxPower = 0;
                        for comIndex = 1:20
                            % Power for 20ms coherent integration
                            comPower = abs(sum(sumPerCode(comIndex:comIndex+19)));
                            % Maximal integration power
                            maxPower = max(maxPower,comPower);
                        end % Search different NH code combiniations
                        fineResult(fineBinIndex) = maxPower;
                    end % for numOfFineBins

                    %--- Find the fine carrier freq. ----------------------------------
                    [~, maxFinBin] = max(fineResult);
                    acqResults.carrFreq(PRN) = fineFreqBins(maxFinBin);
                    % Save code phase acquisition result
                    acqResults.codePhase(PRN) = codePhase;
                    %signal found, if IF =0 just change to 1 Hz to allow processing
                    if(acqResults.carrFreq(PRN) == 0)
                        acqResults.carrFreq(PRN) = 1;
                    end
                else
                    %--- No signal with this PRN --------------------------------------
                    fprintf('. ');
                end
            end
            fprintf(')\n');
        end
        function postProcess(obj)
            %% probe data for validity
            obj.probeData();

            %open file location
            fid = fopen(obj.fnameIQ, 'rb');
            
            if (fid > 0)
                %move to the start of the processing location
                fseek(fid, 2*obj.skipNumBytes, 'bof');


                %Find the number of samples per chipping sequence
                samples_per_code = round(obj.fsamp / (obj.fchip / obj.codeLength));

                %need to read in at least 42ms of data to get good course
                %aquisition
                codeLen = max(42, obj.NonCohTime + 2);
            
                %read the data to perform acquisition
                data = fread(fid, 2*codeLen*samples_per_code, obj.dataType)';
                data = data(1:2:end) + 1i .* data(2:2:end);

                disp("Aquiring satellites");
                %% Aquisition
                acq_results = obj.aquisition(data);
            end
        end
    end
end