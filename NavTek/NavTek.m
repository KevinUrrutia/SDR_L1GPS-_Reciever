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

        %% Tracking
        VSMinterval = 49; %accumulation interval for computing VSM CN0 in ms
        dllCorrelatorSpacing = 0.5; %[chips]
        intTime = 0.001; %integration time for the PLL and DLL
        dllNoiseBandwidth = 1.5; %Hz
        dllDampeningRatio = 0.7; 
        pllNoiseBandwidth = 20; %Hz
        pllDampeningRatio = 0.7;

        %% CN0
        %Accumulation interval in tracking
        CN0_accTime = 0.001;
        %Accumulation interval for computing VSM C/N0
        CN0_VSMinterval = 40;

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

        %calculate the loop filter coefficients
        function [t1, t2] = calcLoppCoeff(LBW, zeta, k)
            %Calculate the natural frequency
            Wn = (LBW * 8 * zeta) / (4 * zeta.^2 + 1);
            t1 = k / (Wn * Wn);
            t2 = 2 * zeta / Wn;
        end

        function plotTracking(channels, tracking_results, msToProcess)
            for channelNr = channels
                if tracking_results(channelNr).status == 'T'      
                    
                    figure(channelNr +200);
                    clf(channelNr +200);
                    set(channelNr +200, 'Name', ['Channel ', num2str(channelNr), ...
                        ' (PRN ', ...
                        num2str(tracking_results(channelNr).PRN), ...
                        ') results']);
                    
                    %% Draw axes ======================================================
                    % Row 1
                    handles(1, 1) = subplot(3, 3, 1);
                    handles(1, 2) = subplot(3, 3, [2 3]);
                    % Row 2
                    handles(2, 1) = subplot(3, 3, 4);
                    handles(2, 2) = subplot(3, 3, [5 6]);
                    % Row 3
                    handles(3, 1) = subplot(3, 3, 7);
                    handles(3, 2) = subplot(3, 3, 8);
                    handles(3, 3) = subplot(3, 3, 9);
                    
                    %% Plot all figures ===============================================
                    
                    timeAxisInSeconds = (1:msToProcess)/1000;
                
                    %----- Discrete-Time Scatter Plot ---------------------------------
                    plot(handles(1, 1), tracking_results(channelNr).I_P,...
                        tracking_results(channelNr).Q_P, ...
                        '.');
                    
                    grid  (handles(1, 1));
                    axis  (handles(1, 1), 'equal');
                    title (handles(1, 1), 'Discrete-Time Scatter Plot');
                    xlabel(handles(1, 1), 'I prompt');
                    ylabel(handles(1, 1), 'Q prompt');
                    
                    %----- Nav bits ---------------------------------------------------
                    plot  (handles(1, 2), timeAxisInSeconds, ...
                        tracking_results(channelNr).I_P);
                    
                    grid  (handles(1, 2));
                    title (handles(1, 2), 'Bits of the navigation message');
                    xlabel(handles(1, 2), 'Time (s)');
                    axis  (handles(1, 2), 'tight');
                    
                    %----- PLL discriminator unfiltered--------------------------------
                    plot  (handles(2, 1), timeAxisInSeconds, ...
                        tracking_results(channelNr).pllDiscr, 'r');
                    
                    grid  (handles(2, 1));
                    axis  (handles(2, 1), 'tight');
                    xlabel(handles(2, 1), 'Time (s)');
                    ylabel(handles(2, 1), 'Amplitude');
                    title (handles(2, 1), 'Raw PLL discriminator');
                    
                    %----- Correlation ------------------------------------------------
                    plot(handles(2, 2), timeAxisInSeconds, ...
                    [sqrt(tracking_results(channelNr).I_E.^2 + ...
                    tracking_results(channelNr).Q_E.^2)', ...
                    sqrt(tracking_results(channelNr).I_P.^2 + ...
                    tracking_results(channelNr).Q_P.^2)', ...
                    sqrt(tracking_results(channelNr).I_L.^2 + ...
                    tracking_results(channelNr).Q_L.^2)'], ...
                    '-*');
                    
                    grid  (handles(2, 2));
                    title (handles(2, 2), 'Correlation results');
                    xlabel(handles(2, 2), 'Time (s)');
                    axis  (handles(2, 2), 'tight');
                        
                    hLegend = legend(handles(2, 2), '$\sqrt{I_{E}^2 + Q_{E}^2}$', ...
                        '$\sqrt{I_{P}^2 + Q_{P}^2}$', ...
                        '$\sqrt{I_{L}^2 + Q_{L}^2}$');
                        
                    %set interpreter from tex to latex. This will draw \sqrt correctly
                    set(hLegend, 'Interpreter', 'Latex');
                        
                    %----- PLL discriminator filtered----------------------------------
                    plot  (handles(3, 1), timeAxisInSeconds, ...
                        tracking_results(channelNr).pllDiscrFilt, 'b');
                        
                    grid  (handles(3, 1));
                    axis  (handles(3, 1), 'tight');
                    xlabel(handles(3, 1), 'Time (s)');
                    ylabel(handles(3, 1), 'Amplitude');
                    title (handles(3, 1), 'Filtered PLL discriminator');
                        
                    %----- DLL discriminator unfiltered--------------------------------
                    plot  (handles(3, 2), timeAxisInSeconds, ...
                        tracking_results(channelNr).dllDiscr, 'r');
                        
                    grid  (handles(3, 2));
                    axis  (handles(3, 2), 'tight');
                    xlabel(handles(3, 2), 'Time (s)');
                    ylabel(handles(3, 2), 'Amplitude');
                    title (handles(3, 2), 'Raw DLL discriminator');
                        
                    %----- DLL discriminator filtered----------------------------------
                    plot  (handles(3, 3), timeAxisInSeconds, ...
                        tracking_results(channelNr).dllDiscrFilt, 'b');
                        
                    grid  (handles(3, 3));
                    axis  (handles(3, 3), 'tight');
                    xlabel(handles(3, 3), 'Time (s)');
                    ylabel(handles(3, 3), 'Amplitude');
                    title (handles(3, 3), 'Filtered DLL discriminator');
                        
                    %----- Plot CNo----------------------------------
                    figure(channelNr +300);
                    clf(channelNr +300);
                    set(channelNr +300, 'Name', ['Channel ', num2str(channelNr), ...
                        ' (PRN ', ...
                        num2str(tracking_results(channelNr).PRN), ...
                        ') CNo']);
                    plot(tracking_results(channelNr).CNo.VSMValue);hold on
                    plot(tracking_results(channelNr).CNo.VSMValue);hold off
                    title('CNo Estimation')
                    ylabel('dB-Hz')
                    xlabel('epoch computation')
                end %if trackResults(channelNr).status == 'T' 
    
            end % for channelNr = channelList
        end

        function CNo = CNoVSM(I, Q, T) 
            Z = I.^2 + Q.^2;
            Zm = mean(Z);
            Zv = var(Z);

            Pavg = sqrt(Zm^2 - Zv);
            Nv = 0.5*(Zm - Pavg);
            CNo = 10 * log10(abs((1/T)*Pavg/(2*Nv)));
        end

        function status = navPartyChk(ndat)
            if (ndat(2) ~= 1)
                ndat(3:26)= -1 .* ndat(3:26);  % Also could just negate
            end
            
            %--- Calculate 6 parity bits ----------------------------------------------
            % The elements of the ndat array correspond to the bits showed in the table
            % 20-XIV (ICD-200C document) in the following way:
            % The first element in the ndat is the D29* bit and the second - D30*.
            % The elements 3 - 26 are bits d1-d24 in the table.
            % The elements 27 - 32 in the ndat array are the received bits D25-D30.
            % The array "parity" contains the computed D25-D30 (parity) bits.
            
            parity(1) = ndat(1)  * ndat(3)  * ndat(4)  * ndat(5)  * ndat(7)  * ...
                        ndat(8)  * ndat(12) * ndat(13) * ndat(14) * ndat(15) * ...
                        ndat(16) * ndat(19) * ndat(20) * ndat(22) * ndat(25);
            
            parity(2) = ndat(2)  * ndat(4)  * ndat(5)  * ndat(6)  * ndat(8)  * ...
                        ndat(9)  * ndat(13) * ndat(14) * ndat(15) * ndat(16) * ...
                        ndat(17) * ndat(20) * ndat(21) * ndat(23) * ndat(26);
            
            parity(3) = ndat(1)  * ndat(3)  * ndat(5)  * ndat(6)  * ndat(7)  * ...
                        ndat(9)  * ndat(10) * ndat(14) * ndat(15) * ndat(16) * ...
                        ndat(17) * ndat(18) * ndat(21) * ndat(22) * ndat(24);
            
            parity(4) = ndat(2)  * ndat(4)  * ndat(6)  * ndat(7)  * ndat(8)  * ...
                        ndat(10) * ndat(11) * ndat(15) * ndat(16) * ndat(17) * ...
                        ndat(18) * ndat(19) * ndat(22) * ndat(23) * ndat(25);
            
            parity(5) = ndat(2)  * ndat(3)  * ndat(5)  * ndat(7)  * ndat(8)  * ...
                        ndat(9)  * ndat(11) * ndat(12) * ndat(16) * ndat(17) * ...
                        ndat(18) * ndat(19) * ndat(20) * ndat(23) * ndat(24) * ...
                        ndat(26);
            
            parity(6) = ndat(1)  * ndat(5)  * ndat(7)  * ndat(8)  * ndat(10) * ...
                        ndat(11) * ndat(12) * ndat(13) * ndat(15) * ndat(17) * ...
                        ndat(21) * ndat(24) * ndat(25) * ndat(26);
            
            %--- Compare if the received parity is equal the calculated parity --------
            if ((sum(parity == ndat(27:32))) == 6)
                
                % Parity is OK. Function output is -1 or 1 depending if the data bits
                % must be inverted or not. The "ndat(2)" is D30* bit - the last  bit of
                % previous subframe. 
                status = -1 * ndat(2);
            else
                % Parity failure
                status = 0;
            end
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

                S = mesh(results);
                grid on;
                xlabel('f_d bin [Hz]')
                ylabel('t_s bin [sec]')
                title(strcat(['C/A for PRN ' int2str(PRN)]))
                pause(1);
                clf(S);

                %% Find Correlation peaks for CA
                %find correlation peaks for CA and teh carrier frequency
                [~, acqCoarseBin] = max(max(results, [], 2));
                %find the code phase for the same correlation peak
                [peakSize, codePhase] = max(max(results));
                %store the GLRT statistic
                acq_results.peakMetric(PRN) = peakSize/sigPower/obj.NonCohTime;
                
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
                    acq_results.carrFreq(PRN) = fineFreqBins(maxFinBin);
                    % Save code phase acquisition result
                    acq_results.codePhase(PRN) = codePhase;
                    %signal found, if IF =0 just change to 1 Hz to allow processing
                    if(acq_results.carrFreq(PRN) == 0)
                        acq_results.carrFreq(PRN) = 1;
                    end
                else
                    %--- No signal with this PRN --------------------------------------
                    fprintf('. ');
                end
            end
            fprintf(')\n');
        end

        function channel = sortAcquisition(obj, acqResults)
            channel = [];
            channel.PRN = 0;
            channel.acquiredFreq = 0;
            channel.codePhase = 0;
            channel.status          = '-';

            [~, PRNindexes] = sort(acqResults.peakMetric, 2, 'descend');
            channel = repmat(channel, 1, obj.numChannels);
            for ii = 1:min([obj.numChannels, sum(acqResults.carrFreq ~= 0)])
                channel(ii).PRN          = PRNindexes(ii);
                channel(ii).acquiredFreq = acqResults.carrFreq(PRNindexes(ii));
                channel(ii).codePhase    = acqResults.codePhase(PRNindexes(ii));
                
                % Set tracking into mode (there can be more modes if needed e.g. pull-in)
                channel(ii).status       = 'T';
            end
        end

        function dispFineCAResults(obj, channel)
            fprintf("\nFINE ACQUISITION RESULTS\n");
            fprintf('  *=========*=====*===============*===========*=============*========*\n');
            fprintf(  '| Channel | PRN |   Frequency   |  Doppler  | Code Offset | Status |\n');
            fprintf(  '*=========*=====*===============*===========*=============*========*\n');

            for channelNr = 1 : obj.numChannels
                if (channel(channelNr).status ~= '-')
                    fprintf('|      %2d | %3d |  %2.5e |   %5.0f   |    %6d   |     %1s  |\n', ...
                            channelNr, ...
                            channel(channelNr).PRN, ...
                            channel(channelNr).acquiredFreq, ...
                            channel(channelNr).acquiredFreq - obj.fIF, ...
                            channel(channelNr).codePhase, ...
                            channel(channelNr).status);
                else
                    fprintf('|      %2d | --- |  ------------ |   -----   |    ------   |   Off  |\n', ...
                            channelNr);
                end
            end

            fprintf('*=========*=====*===============*===========*=============*========*\n\n');
        end

        function tracking_results = tracking(obj, fid, channel)
            % intialize tracking structure array
            tracking_results.status = '-';
            %in Record of C/A start
            tracking_results.absoluteSamples = zeros(1, obj.msToProcess);
            %Freq of PRN code
            tracking_results.codeFreq = inf(1, obj.msToProcess);
            %Freq of tracked carrier wave
            tracking_results.carrFreq = inf(1, obj.msToProcess);
            %in phase correlatores
            tracking_results.I_P  = zeros(1, obj.msToProcess);
            tracking_results.I_E = zeros(1, obj.msToProcess);
            tracking_results.I_L = zeros(1, obj.msToProcess);
            %Quadrature correlators
            tracking_results.Q_P = zeros(1, obj.msToProcess);
            tracking_results.Q_E = zeros(1, obj.msToProcess);
            tracking_results.Q_L = zeros(1, obj.msToProcess);
            %Loop Descriminators
            tracking_results.dllDescr = inf(1, obj.msToProcess);
            tracking_results.dllDescrFilt = inf(1, obj.msToProcess);
            tracking_results.pllDescr = inf(1, obj.msToProcess);
            tracking_results.pllDescrFilt =inf(1, obj.msToProcess);
            %Remainder of code and carr phase
            tracking_results.remCodePhase = inf(1, obj.msToProcess);
            tracking_results.remCarrPhase = inf(1, obj.msToProcess);
            %C/N0
            tracking_results.CN0.VSMValue = zeros(1, floor(obj.msToProcess/obj.CN0_VSMinterval));
            tracking_results.CN0.VSMidx = zeros(1, floor(obj.msToProcess/obj.CN0_VSMinterval));

            %copy struct for all tracked satellites
            tracking_results = repmat(tracking_results, 1, 12);
            
            %---------------Intitalize tracking variables
            codePeriods = obj.msToProcess;
            
            %---------------DLL variables
            earlyLateSPC = obj.dllCorrelatorSpacing;
            %summation interval
            PDIcode = obj.intTime;
            %calculate loop coeffiecient values
            %calculate loop coeffiecient values
            [t1_code, t2_code] = obj.calcLoppCoeff(obj.dllNoiseBandwidth, obj.dllDampeningRatio, 1);
            
            %--------------PLL variables
            PDIcarr = obj.intTime;
            [t1_carr, t2_carr] = obj.calcLoppCoeff(obj.pllNoiseBandwidth, obj.pllDampeningRatio, 0.25);

            for channelNr = 1:obj.numChannels
                %acquistion happened correctly
                if channel(channelNr).PRN ~= 0
                    tracking_results(channelNr).PRN = channel(channelNr).PRN;

                    %move to the appropriate start of the file, this should be after the acquisition results have been grabbed
                    fseek(fid, 2*(obj.skipNumBytes + channel(channelNr).codePhase - 1), 'bof');
            
                    %generate C/A code sampled 1x/chip
                    ca_gen = obj.generateGoldSeq(obj.fsamp, obj.fchip, obj.codeLength, channel(channelNr).PRN, 0);
                    caCode = [ca_gen(1023), ca_gen, ca_gen(1)];

                    %---Perform Initializations for NCO
                    %define initial code frequency basis for NCO
                    codeFreq = obj.fchip;
                    %define residual code phase
                    remCodePhase = 0;
                    %define carrier frequency used over entire tracking period
                    carrFreq = channel(channelNr).acquiredFreq;
                    carrFreqBasis = channel(channelNr).acquiredFreq;
                    %define residual carrier phase
                    remCarrPhase = 0;

                    %code tracking loop parameters
                    oldCodeNco   = 0.0;
                    oldCodeError = 0.0;
            
                    %carrier/Costas loop parameters
                    oldCarrNco   = 0.0;
                    oldCarrError = 0.0;
            
                    %CN0 computation
                    vsmCnt = 0; 

                    %-------Process for the specificied code periods 
                    for loopCnt = 1:codePeriods
                        %read in the following block of data
                        tracking_results(channelNr).absoluteSample(loopCnt) = (ftell(fid)) / 2;
            
                        %update the phase step based on code freq (changing) and
                        %sampling freq (fixed)
                        codePhaseStep = codeFreq / obj.fsamp;

                        %fine the size of a "block" or code period in whole samples
                        blksize = ceil((obj.codeLength - remCodePhase) / codePhaseStep);

                        %read in the appropriate number of samples based on the blksize
                        [signal, samples] = fread(fid, 2*blksize, obj.dataType);

                        signal = signal';
                        signal = signal(1:2:end) + 1i.*signal(2:2:end);
            
                        %if sufficient samples were not read in end the program
                        if(samples ~= 2*blksize)
                            fprintf("\nNot able to read the specified number of samples for tracking, exiting tracking.\n");
                            return;
                        end

                        %% Setup the code phase tracking info
                        %save remCodePhase for current correlation
                        tracking_results(channelNr).remCodePhase(loopCnt) = remCodePhase;

                        % Define index into early code vector
                        tcode = (remCodePhase-earlyLateSPC) : codePhaseStep : ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSPC);
                        tcode2      = ceil(tcode) + 1;
                        earlyCode   = caCode(tcode2);

                        % Define index into late code vector
                        tcode = (remCodePhase+earlyLateSPC) : codePhaseStep : ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSPC);
                        tcode2 = ceil(tcode) + 1;
                        lateCode = caCode(tcode2);

                        %Define index to prompt code vector
                        tcode       = remCodePhase : codePhaseStep : ((blksize-1)*codePhaseStep+remCodePhase);
                        tcode2      = ceil(tcode) + 1;
                        promptCode  = caCode(tcode2);

                        %remaining code phase for each tracking update
                        remCodePhase = (tcode(blksize) + codePhaseStep) - obj.codeLength;

                        %% Generate the carrier frequency to mix to baseband
                        tracking_results(channelNr).remCarrPhase(loopCnt) = remCarrPhase;

                        %Get the argument to the sin and cos functions
                        time = (0:blksize) ./ obj.fsamp;
                        trigarg = ((carrFreq * 2 * pi) .* time) + remCarrPhase;
                        %Remaining carrier phase for each tracking update
                        remCarrPhase = rem(trigarg(blksize + 1), (2*pi));

                        %compute the signal to mixe the collected data to baseband
                        carrsig = exp(-1i .* trigarg(1:blksize));

                        %% Correlate the signal and get the size accumulated values
                        iBaseBand = real(carrsig .* signal);
                        qBaseBand = imag(carrsig .* signal);
            
                        %Get the early, late and prompt values
                        I_E = sum(earlyCode .* iBaseBand);
                        Q_E = sum(earlyCode .* qBaseBand);
                        I_P = sum(promptCode .* iBaseBand);
                        Q_P = sum(promptCode .* qBaseBand);
                        I_L = sum(lateCode .* iBaseBand);
                        Q_L = sum(lateCode .* qBaseBand);

                         %% Find PLL error and update carrier NCO

                        %implement carrier loop discriminator (phase detection)
                        carrError = atan(Q_P / I_P) / (2 * pi);
                        %Implement carrier loop filter and generate NCO command
                        carrNco = oldCarrNco + (t2_carr / t1_carr) * (carrError - oldCarrError) + carrError * (PDIcarr / t1_carr);
                        oldCarrNco = carrNco;
                        oldCarrError = carrError;

                        %save carrier frequency for current correlation
                        tracking_results(channelNr).carrFreq(loopCnt) = carrFreq;

                        %modify carrier freq based on NCO command
                        carrFreq = carrFreqBasis + carrNco;
                    
                        %% Find DLL error and update code NCO
                        codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
                        % Implement code loop filter and generate NCO command
                        codeNco = oldCodeNco + (t2_code/t1_code) * (codeError - oldCodeError) + codeError * (PDIcode/t1_code);
                        oldCodeNco   = codeNco;
                        oldCodeError = codeError;

                        % Save code frequency for current correlation
                        tracking_results(channelNr).codeFreq(loopCnt) = codeFreq;
            
                        % Modify code freq based on NCO command
                        codeFreq = obj.fchip - codeNco;

                        %% Record various measures to show in postprocessing ----------------------
                        tracking_results(channelNr).dllDiscr(loopCnt)       = codeError;
                        tracking_results(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
                        tracking_results(channelNr).pllDiscr(loopCnt)       = carrError;
                        tracking_results(channelNr).pllDiscrFilt(loopCnt)   = carrNco;
            
                        tracking_results(channelNr).I_E(loopCnt) = I_E;
                        tracking_results(channelNr).I_P(loopCnt) = I_P;
                        tracking_results(channelNr).I_L(loopCnt) = I_L;
                        tracking_results(channelNr).Q_E(loopCnt) = Q_E;
                        tracking_results(channelNr).Q_P(loopCnt) = Q_P;
                        tracking_results(channelNr).Q_L(loopCnt) = Q_L;

                        %% CNo calculation --------------------------------------
                        if (rem(loopCnt,obj.CN0_VSMinterval)==0)
                            vsmCnt=vsmCnt+1;
                            CNoValue=obj.CNoVSM(tracking_results(channelNr).I_P(loopCnt-obj.CN0_VSMinterval+1:loopCnt),...
                                tracking_results(channelNr).Q_P(loopCnt-obj.CN0_VSMinterval+1:loopCnt),obj.CN0_accTime);
                            tracking_results(channelNr).CNo.VSMValue(vsmCnt)=CNoValue;
                            tracking_results(channelNr).CNo.VSMIndex(vsmCnt)=loopCnt;
                        end

                        clc;
                        fprintf("\nTracking.....%2d\n",tracking_results(channelNr).PRN);
                    end
                    tracking_results(channelNr).status  = channel(channelNr).status;
                end
            end
        end

        function [nav_sol, ephem] = compute_position(obj, tracking_results)
            %Warning a minimum of 30s is needed to calculate position,
            %considering that at least 3 subframes are needed to find the
            %satellites coordinates, each subframe is 6s long 

            if (obj.msToProcess < 36000)
                error("Record is too short to process, Need at least 36s to compute navigation solution");
            end

            %pre allocate space for nagiation results
            %this process needs to happen for each satellite being tracked
            %Preable of each subframe
            subFrameStart = inf(1, obj.numChannels);

            %Time of week (GPS) Time being decoded from the satellites
            TOW = inf(1, obj.numChannels);

            %create a list of satellites that actually tracked to the end
            trackedChan = find([tracking_results.status] ~= '-');

            % Decode the ephemeris from the subframes
            for channelNr = trackedChan
                %get the current PRN number
                PRN = tracking_results(channelNr).PRN;
                fprintf("Now decoding PRN %02d navigation message.......\n", PRN);

                %Decode subframes to get ephemeris data
                [ephem(PRN), subFrameStart(channelNr), TOW(channelNr)] = obj.decode(tracking_results(channelNr).I_P);
            end

        end 

        function [eph, subframeStart, TOW] = decode(obj, I_P) 
            %create ephermeris data structure 
            eph.idValid(1:5) = zeros(1, 5); %valid subframe
            eph.PRN = [];
            eph.week = [];
            eph.TOW = [];

            %subframe 1 elements. contains WN, clock corrections, health
            %and accuracy
            eph.weekNumber = [];
            eph.accuracy = [];
            eph.health = [];
            eph.T_GD = [];
            eph.IODC = [];
            eph.t_oc = [];
            eph.a_f2 = [];
            eph.a_f1 = [];
            eph.a_f0 = [];

            %subframe 2. contains first part of ephemeris parameters
            eph.IODE_2 = []; %Issue of data (ephemeris)
            eph.C_rs = []; % Amplitude of sine harmonic correction term to the orbit radius
            eph.deltan = []; %mean motion difference from computed value
            eph.M_0 = []; %Mean anomoly
            eph.C_uc = []; %Amplitude of cosine harmonic correction term to the argument of latitiude
            eph.e = []; %ecentricity
            eph.C_us = []; %Amplitude of sine harmonic correction term to the argument of latitude
            eph.t_oe = []; %reference time ephemeris

            %subframe 3. contains second part of ephemeris parameters
            eph.C_ic = []; %Amplitude of cosine harmonic correction term to the angle of inclination
            eph.omega_0 = []; %Longitude of Ascending node of orbit plane at weekly epoch
            eph.C_is = []; %Amplitude of sine harmonic correction term to the angle of inclination
            eph.i_0 = []; %Inclination angle at reference time
            eph.C_rc = []; %Amplitude of cosine harmonic correction term to the orbit radius
            eph.omega = []; %argument of perigee
            eph.omega_dot = []; %rate of rigt ascension
            eph.IODE_3 = []; %Issue of data (ephemeris)
            eph.IDOT = []; %Rate of inclination angle

            subframeStart = inf;
            TOW = inf;

            preamble = [1 -1 -1 -1 1 -1 1 1];
            preamble_ms = kron(preamble, ones(1, 20));

            %just in case we want to start searching later in the tracking
            %loops
            searchStartOffset = 0;
            %grab valid bits in I_P
            bits = I_P(1 + searchStartOffset : end);
            
            %create vector 1 and -1 to match the preamble
            bits(bits > 0)  =  1;
            bits(bits <= 0) = -1;

            %correlate the tracking results with preamble to find instance
            %of it in the bits 
            clear index
            clear index2
            tlmXcorrResult = xcorr(bits, preamble_ms);
            xcorrLength = (length(tlmXcorrResult) +  1) /2;

            %get the index where the preamble starts
            index = find(...
                abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 153)' + ...
                searchStartOffset;
            %throw out indecies that might cause boundary conditions
            index = index((index>40 & index<obj.msToProcess - (20 * 60 -1)) == 1);

            for i = 1:size(index, 1)
                %find the differences in time for the preable detection if
                %it is approx 6s worth then we can move on to further
                %checks
                index2 = index - index(i);
                if(~isempty(find(index2 == 6000, 1)))
                    %re read bit values for preamble verification 
                    %preamble is verified by checking the parity of the
                    %first two words in the subframe. In total 62 bits
                    %need to be read:
                    % 2 bits from the previous subframe are needed for
                    % parity check
                    %60 bits for the first two 60bit words (TLM and HOW words)
                    %index is to the start of the TLM word
                    bits = I_P(index(i)-40 : index(i) + 20 * 60 -1)';

                    %combine 20 values of each bit
                    bits = reshape(bits, 20, (size(bits, 1) / 20));
                    bits = sum(bits);

                    %create vector 1 and -1 to match the preamble
                    bits(bits > 0) = 1;
                    bits(bits <= 0) = -1;

                    if (obj.navPartyChk(bits(1:32)) ~= 0) && ...
                            (obj.navPartyChk(bits(31:62)) ~= 0)
                        %if the parity check is fine then save the pramble
                        %starting positing for this channel
                        subframeStart = index(i);
                        break;
                    end
                end   
            end

            %Exclude the channel if there is no preamble
            if subframeStart == inf
                disp('Could not find valid preambles in channel! ');
                return
            end
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
                
                if any(acq_results.carrFreq)
                    channel = obj.sortAcquisition(acq_results);
                    obj.dispFineCAResults(channel);

                else 
                    warning("No Satellites to track");
                    return;
                end
                %get rid of acquisiton plots to start tracking
                close all;

                %% Tracking
                track_results = obj.tracking(fid, channel);
                obj.plotTracking(1:obj.numChannels, track_results, obj.msToProcess);
                disp("Enter to continue to naviagtion solution");
                pause;
                close all;

                %% Position solution

            end
        end
    end
end