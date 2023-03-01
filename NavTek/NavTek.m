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
        fdmin = -6000; %[Hz]
        fdmax = 6000; %[Hz]
        delta_fd = 500; %[Hz]
        %NonCoherent integration times after 1ms coherent integration
        NonCohTime = 30; 
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

        %% navigation solution settings
        navSolPeriod = 500;
        elMask = 5; %degress [0-90] Must be able to see satellites at low elevations

        %% Constants
        c = 299792458;
        startOffset = 68.802; %[ms] Initial sign. travel time
        mu = 3.986005e14; %(m^3/s^2) WGS 84 value of the earth's gravitational constant for GPS User
        omega_e = 7.29212e-5; %(rad/s) WGS 84 value of the earth's rotation rate
        gpsPi = 3.1415926535898; %pi used in GPS coordinate system
        F = -4.442807633e-10; %Constant, [sec/(meter)^(1/2)]
        half_week = 302400; %seconds
        a = 6378137.0; %mean radius of the earth
        flat = 1/298.257223563; % flattening factor of the earth
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

        %Calculate the noise equivalent bandwidth
        function [tau1, tau2] = calcLoppCoeff(f_n, zeta, k)
            %Calculate the natural frequency
            B = (f_n * 8 * zeta) / (4 * zeta.^2 + 1);
            tau1 = k / (B * B);
            tau2 = 2 * zeta / B;
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
                end 
    
            end 
        end

        function CNo = CNoVSM(I, Q, T) 
            X_p = I.^2 + Q.^2;
            X_p_mean = mean(X_p); %calculation of first moment
            X_p_var = var(X_p); %calculation of second moment

            Pavg = sqrt(X_p_mean^2 - X_p_var); %average for one std away from the mean
            Nv = 0.5*(X_p_mean - Pavg);
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

        function num = twosComp2dec(binNum)
            num = bin2dec(binNum);
            if(num(1) == '1')
                num = num - 2^size(binNum, 2); %takes two's comp if negative
            end
        end

        function result = invert(data)
            dataLength = length(data);
            temp(1:dataLength) = '1';

            invertMask = bin2dec(char(temp));

            result = dec2bin(bitxor(bin2dec(data), invertMask), dataLength);
        end

        function word = checkPhase(word, D30Star)
            if D30Star == '1'
                % Data bits must be inverted
                word(1:24) = NavTek.invert(word(1:24));
            end
        end
    end

    methods
        function obj = NavTek()
            addpath("googleearth_matlab/googleearth/");
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
                [X_p, count] = fread(fid, [1, 2 * 100 * sampsPerCode], obj.dataType);

                fclose(fid);

                if (count < 2 * 100 * sampsPerCode)
                    error("Error not enough data to load in");
                end

                %initialize plots
                timeScale = 0: 1/obj.fsamp : 5e-3; 
                f = figure;
                f.Position(3:4) = [1000, 1000];

                %% Time domain plot
                X_p = X_p(1:2:end) + 1i .* X_p(2:2:end);
                subplot(3, 2, 4);
                plot(1000 * timeScale(1:round(sampsPerCode/2)), ...
                real(X_p(1:round(sampsPerCode/2))));

                axis tight;    grid on;
                title ('Time domain plot (I)');
                xlabel('Time (ms)'); ylabel('Amplitude');

                subplot(3, 2, 3);
                plot(1000 * timeScale(1:round(sampsPerCode/2)), ...
                imag(X_p(1:round(sampsPerCode/2))));

                axis tight;    grid on;
                title ('Time domain plot (Q)');
                xlabel('Time (ms)'); ylabel('Amplitude');

                %% Frequency domain plot
                subplot(3,2,1:2);
                [sigspec,freqv]=pwelch(X_p, 32768, 2048, 32768, obj.fsamp,'twosided');
                plot(([-(freqv(length(freqv)/2:-1:1));freqv(1:length(freqv)/2)])/1e6, ...
                10*log10([sigspec(length(freqv)/2+1:end);
                sigspec(1:length(freqv)/2)]));
                axis tight;
                grid on;
                title ('Frequency domain plot');
                xlabel('Frequency (MHz)'); ylabel('Magnitude');

                %% Histogram
                subplot(3, 2, 6);
                histogram(real(X_p), -128:128)
                dmax = max(abs(X_p)) + 1;
                axis tight;     adata = axis;
                axis([-dmax dmax adata(3) adata(4)]);
                grid on;        title ('Histogram (I)');
                xlabel('Bin');  ylabel('Number in bin');
            
                subplot(3, 2, 5);
                histogram(imag(X_p), -128:128)
                dmax = max(abs(X_p)) + 1;
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
            sampsPerCode = round(obj.fsamp / (obj.fchip / obj.codeLength));
            %find sampling period 
            tau = 1 / obj.fsamp;
            %Find phase of 2ms 
            phase = (0: (sampsPerCode*2-1)) * 2 * pi * tau;
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
            fineStep = 25;
            %Number of frequency bins for fine aquisition
            numFineBins = round(obj.delta_fd / fineStep) + 1;
            %Carrier frequencies of fine frequency bins
            fineFreqBins = zeros(1, numFineBins);
            %Correlation values for all find freqeucny bins
            fineResult = zeros(1, numFineBins);
            %Coherent integration of 40 codes
            RFaccum = zeros(1, 40);
            %phase points of the local carrier wave 
            finePhase = (0 : (40 * sampsPerCode-1)) * 2 * pi * tau;

            %==Input signal power for GRLT statistic calculation
            sigPower = sqrt(var(data(1:sampsPerCode)) * sampsPerCode);

            %Preform search for all listed satellite numbers....
            fprintf('(');
            for PRN = obj.SVn
                %% Course Aquisition
                %generate PRN sqequences
                caCodes = obj.generateGoldSeq(obj.fsamp, obj.fchip, obj.codeLength, PRN, true);
                %add zero padding samples
                caCode2ms = [caCodes, zeros(1, sampsPerCode)];
%                 caCode2ms = resample(caCodes, obj.fsamp*sampsPerCode, obj.fchip);
                %search results of all frequency bins and code shifts
                results = zeros(numFreqBins, sampsPerCode * 2);
                %perform DFT of CA codes
                C_c = conj(fft(caCode2ms));

                %make correlation for all frequency bins
                %frequnecy bin values:
                fdVec = -obj.fdmax:obj.delta_fd:obj.fdmax;
                ii = 1;
                for freqBinIndex = 1:numFreqBins
                    %generate carrier wave frequency grid
                    courseFreqBins(freqBinIndex) = obj.fIF + fdVec(ii);
                    
                    %generate local sine and cosine
                    Xp_c = exp(-1i * courseFreqBins(freqBinIndex) * phase);

                    %Complete Correlation
                    for nonCohIndex = 1: obj.NonCohTime
                        %Take 2ms of input data to do correlation
                        Xp = data((nonCohIndex - 1) * sampsPerCode + 1: (nonCohIndex + 1) * sampsPerCode);
                        
                        %remove the carrier from the signal
                        F_X = fft(Xp_c .* Xp);
                        %Multiplication in frequency domain corresponds to
                        %correlation in time domain
                        X_c = F_X .* C_c;
                        %perform DFT and store correlation results
                        S_k = abs(ifft(X_c));
                        %Non Coherent Integration
                        results(freqBinIndex, :) = results(freqBinIndex, :) + S_k;
                    end
                    ii = ii + 1;
                end 

                S = mesh(results);
                grid on;
                ylabel('f_d bin [Hz]')
                xlabel('t_s bin [sec]')
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
                    codeIdx = floor((tau * (0 : 40*sampsPerCode -1)) / (1/obj.fchip));
                    % C/A code samples
                    caCode40ms = caCode(rem(codeIdx, obj.codeLength) + 1); % add one for matlab indexing
                    % Take 40ms incoming signal for fine acquisition
                    X_p40cm = data(codePhase:codePhase + 40*sampsPerCode - 1);

                    %search fine freqency bins
                    fine_fdVec = -obj.delta_fd/2:fineStep:obj.delta_fd/2;
                    ii = 1;
                    for fineBinIdx = 1 : numFineBins
                        %--- Correlation for each code --------------------------------
                        % Carrier frequencies of the frequency bins
                        fineFreqBins(fineBinIdx) = courseFreqBins(acqCoarseBin) + fine_fdVec(ii);

                        % Local carrier signal
                        X_p_c40cm = exp(-1i * fineFreqBins(fineBinIdx) * finePhase);
                        % Wipe off code and carrier from incoming signals
                        X_p = X_p40cm .* caCode40ms .* X_p_c40cm;

                        % Coherent integration for each code
                        for jj = 1:40
                            RFaccum(jj) = sum(X_p(sampsPerCode * (jj - 1) + 1 : sampsPerCode * jj));
                        end

                        % 20 cases of Nav bit edge
                        maxPower = 0;
                        cases = 20;
                        for comIndex = 1:cases
                            % Power for 20ms coherent integration
                            comPower = abs(sum(RFaccum(comIndex:comIndex+cases-1)));
                            % Max integration power
                            maxPower = max(maxPower,comPower);
                        end % Search different NH code combiniations
                        fineResult(fineBinIdx) = maxPower;
                        ii = ii + 1;
                    end 

                    [~, maxFinBin] = max(fineResult); %fine carrier frequency
                    acq_results.carrFreq(PRN) = fineFreqBins(maxFinBin);
                    % Save code phase acquisition result
                    acq_results.codePhase(PRN) = codePhase;
                    if(acq_results.carrFreq(PRN) == 0) % signal has been found
                        acq_results.carrFreq(PRN) = 1;
                    end
                else
                    fprintf('. '); %No signal found
                end
            end
            fprintf(')\n');
        end

        function sat = rankAcquisition(obj, acqResults)
            sat = [];
            sat.PRN = 0;
            sat.acquiredFreq = 0;
            sat.codePhase = 0;
            sat.status          = '-';

            [~, PRNindexes] = sort(acqResults.peakMetric, 2, 'descend');
            sat = repmat(sat, 1, obj.numChannels);
            for ii = 1:min([obj.numChannels, sum(acqResults.carrFreq ~= 0)])
                sat(ii).PRN          = PRNindexes(ii);
                sat(ii).acquiredFreq = acqResults.carrFreq(PRNindexes(ii));
                sat(ii).codePhase    = acqResults.codePhase(PRNindexes(ii));
                
                % Set tracking into mode (there can be more modes if needed e.g. pull-in)
                sat(ii).status       = 'T';
            end
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
            tracking_results.r_code_phase = inf(1, obj.msToProcess);
            tracking_results.r_carr_phase = inf(1, obj.msToProcess);
            %C/N0
            tracking_results.CN0.VSMValue = zeros(1, floor(obj.msToProcess/obj.CN0_VSMinterval));
            tracking_results.CN0.VSMidx = zeros(1, floor(obj.msToProcess/obj.CN0_VSMinterval));

            %copy struct for all tracked satellites
            tracking_results = repmat(tracking_results, 1, obj.numChannels);
            
            %---------------Intitalize tracking variables
            codePeriods = obj.msToProcess;
            
            %---------------DLL variables
            earlyLateSPC = obj.dllCorrelatorSpacing;
            %summation interval
            int_code = obj.intTime;
            %calculate loop coeffiecient values
            %calculate loop coeffiecient values
            [tau1_code, tau2_code] = obj.calcLoppCoeff(obj.dllNoiseBandwidth, obj.dllDampeningRatio, 1);
            
            %--------------PLL variables
            int_carr = obj.intTime;
            [tau1_carr, tau2_carr] = obj.calcLoppCoeff(obj.pllNoiseBandwidth, obj.pllDampeningRatio, 0.25);

            for ii = 1:obj.numChannels
                %acquistion happened correctly
                if channel(ii).PRN ~= 0
                    tracking_results(ii).PRN = channel(ii).PRN;

                    %move to the appropriate start of the file, this should be after the acquisition results have been grabbed
                    fseek(fid, 2*(obj.skipNumBytes + channel(ii).codePhase - 1), 'bof');
            
                    %generate C/A code sampled 1x/chip
                    caCode_gen = obj.generateGoldSeq(obj.fsamp, obj.fchip, obj.codeLength, channel(ii).PRN, 0);
                    caCode = [caCode_gen(1023), caCode_gen, caCode_gen(1)];

                    %Initialize NCO
                    %define initial code frequency basis for NCO
                    NCO = obj.fchip; %will be updated based on tracking loops
                    %define residual code phase
                    r_cp = 0;
                    %define carrier frequency used over entire tracking period
                    f_c = channel(ii).acquiredFreq;
                    carrFreqBasis = channel(ii).acquiredFreq;
                    %define residual carrier phase
                    r_carr_phase = 0;

                    %code tracking loop parameters
                    prev_CodeNco   = 0;
                    prev_CodeError = 0;
            
                    %carrier/Costas loop parameters
                    prev_CarrNco   = 0;
                    prev_carr_disc = 0;
            
                    %CN0 computation
                    power_cnt = 0; 

                    %-------Process for the specificied code periods 
                    for jj = 1:codePeriods
                        %read in the following block of data
                        tracking_results(ii).absoluteSample(jj) = (ftell(fid)) / 2;
            
                        %update the phase step based on code freq (changing) and
                        %sampling freq (fixed)
                        delta_code_phase = NCO / obj.fsamp;

                        %find the size of code period in whole samples
                        codePeriod = ceil((obj.codeLength - r_cp) / delta_code_phase);

                        %read in the appropriate number of samples based on
                        %the code Period
                        [X_p, samp_num] = fread(fid, 2*codePeriod, obj.dataType);

                        X_p = X_p';
                        X_p = X_p(1:2:end) + 1i.*X_p(2:2:end);
            
                        %if sufficient samples were not read in end the program
                        if(samp_num < 2*codePeriod)
                            disp()
                            disp("Not able to read the specified number of samples for tracking, exiting tracking.")
                            return;
                        end

                        %save r_code_phase
                        tracking_results(ii).r_code_phase(jj) = r_cp;

                        % Define early code vector
                        code = (r_cp - earlyLateSPC) : delta_code_phase : ((codePeriod-1) * delta_code_phase + r_cp -earlyLateSPC);
                        code2      = ceil(code) + 1;
                        early   = caCode(code2);

                        % Define late code vector
                        code = (r_cp + earlyLateSPC) : delta_code_phase : ((codePeriod-1) * delta_code_phase +r_cp + earlyLateSPC);
                        code2 = ceil(code) + 1;
                        late = caCode(code2);

                        %Define prompt code vector
                        code       = r_cp : delta_code_phase : ((codePeriod-1) * delta_code_phase + r_cp);
                        code2      = ceil(code) + 1;
                        prompt  = caCode(code2);

                        %remaining code phase for each tracking update
                        r_cp = (code(codePeriod) + delta_code_phase) - obj.codeLength;

                        %% Generate the carrier frequency 
                        tracking_results(ii).r_carr_phase(jj) = r_carr_phase;

                        %Get the argument to the sin and cos functions
                        t_k = (0:codePeriod) ./ obj.fsamp; %should be at the sampling rate
                        phi = ((f_c * 2 * pi) .* t_k) + r_carr_phase;
                        %Remaining carrier phase for each tracking update
                        r_carr_phase = rem(phi(codePeriod + 1), (2*pi));

                        %compute the signal to mixe the collected data to baseband
                        X_c = exp(-1i .* phi(1:codePeriod));

                        %% Correlate the signal and get the size accumulated values
                        I = real(X_c .* X_p);
                        Q = imag(X_c .* X_p);
            
                        %Get the early, late and prompt values
                        I_E = sum(early .* I);
                        Q_E = sum(early .* Q);
                        I_P = sum(prompt .* I);
                        Q_P = sum(prompt .* Q);
                        I_L = sum(late .* I);
                        Q_L = sum(late .* Q);

                         %% Find PLL error and update carrier NCO

                        %implement carrier loop discriminator (phase detection)
                        carr_disc = atan(Q_P / I_P) / (2 * pi);
                        %Implement carrier loop filter and generate NCO command
                        carrNco = prev_CarrNco + (tau2_carr / tau1_carr) * (carr_disc - prev_carr_disc) + carr_disc * (int_carr / tau1_carr);
                        prev_CarrNco = carrNco;
                        prev_carr_disc = carr_disc;

                        %save carrier frequency for current correlation
                        tracking_results(ii).carrFreq(jj) = f_c;

                        %modify carrier freq based on NCO command
                        f_c = carrFreqBasis + carrNco;
                    
                        %% Find DLL error and update code NCO
                        code_disc = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
                        % Implement code loop filter and generate NCO command
                        codeNco = prev_CodeNco + (tau2_code/tau1_code) * (code_disc - prev_CodeError) + code_disc * (int_code/tau1_code);
                        prev_CodeNco   = codeNco;
                        prev_CodeError = code_disc;

                        % Save code frequency for current correlation
                        tracking_results(ii).codeFreq(jj) = NCO;
            
                        % Modify code freq based on NCO command
                        NCO = obj.fchip - codeNco;

                        %% Record various measures to show in postprocessing ----------------------
                        tracking_results(ii).dllDiscr(jj)       = code_disc;
                        tracking_results(ii).dllDiscrFilt(jj)   = codeNco;
                        tracking_results(ii).pllDiscr(jj)       = carr_disc;
                        tracking_results(ii).pllDiscrFilt(jj)   = carrNco;
            
                        tracking_results(ii).I_E(jj) = I_E;
                        tracking_results(ii).I_P(jj) = I_P;
                        tracking_results(ii).I_L(jj) = I_L;
                        tracking_results(ii).Q_E(jj) = Q_E;
                        tracking_results(ii).Q_P(jj) = Q_P;
                        tracking_results(ii).Q_L(jj) = Q_L;

                        %% CNo calculation --------------------------------------
                        if (rem(jj,obj.CN0_VSMinterval)==0)
                            power_cnt=power_cnt+1;
                            CNoValue=obj.CNoVSM(tracking_results(ii).I_P(jj-obj.CN0_VSMinterval+1:jj), tracking_results(ii).Q_P(jj-obj.CN0_VSMinterval+1:jj),obj.CN0_accTime);
                            tracking_results(ii).CNo.VSMValue(power_cnt)=CNoValue;
                            tracking_results(ii).CNo.VSMIndex(power_cnt)=jj;
                        end

                        clc;
                        fprintf("\nTracking.....%2d\n",tracking_results(ii).PRN);
                    end
                    tracking_results(ii).status  = channel(ii).status;
                end
            end
        end

        function [navSol, ephem] = compute_position(obj, tracking_results)
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
            trackedChan = trackedChan(1:end-1);

            % Decode the ephemeris from the subframes
            for channelNr = trackedChan
                %get the current PRN number
                PRN = tracking_results(channelNr).PRN;
                fprintf("Now decoding PRN %02d navigation message.......\n", PRN);

                %Decode subframes to get ephemeris data
                [ephem(PRN), subFrameStart(channelNr), TOW(channelNr)] = obj.decode(tracking_results(channelNr).I_P);


                %Exclude a satellite if it does not have the necessary
                %navigation data

                if(isempty(ephem(PRN).IODC) || isempty(ephem(PRN).IODE_2) || isempty(ephem(PRN).IODE_3) || (ephem(PRN).health ~=0))
                    %Exclude channel from the active list
                    trackedChan = setdiff(trackedChan, channelNr);
                    fprintf(' decoding fails for PRN %02d !\n', PRN);
                else
                    fprintf(' Three requistite messages for PRN %02d all decoded!\n', PRN);
                end
            end

             %check if the number of satellites is still above 4
             if(isempty(trackedChan) || (size(trackedChan, 2) < 4))
                 disp('Too few satellites with ephemeris data for position calculations Exiting!');
                 nav_sol = [];
                 ephem = [];
             end
    
            %set the measurement-time point and step
            %find start and end measurement point locations in IF signal stream
            %with available measurements
            sampleStart = zeros(1, obj.numChannels);
            sampleEnd = inf(1, obj.numChannels);
    
            for channelNr = trackedChan
                sampleStart(channelNr) = tracking_results(channelNr).absoluteSample(subFrameStart(channelNr));
                sampleEnd(channelNr) = tracking_results(channelNr).absoluteSample(end);
            end
    
            %second terms is to make space to avoid index exceeds matrix
            %dimensions
            sampleStart = max(sampleStart) + 1;
            sampleEnd = min(sampleEnd) + 1;
    
            %measurement step in units of IF samples
            measSampleStep = fix(obj.fsamp * obj.navSolPeriod/1000);
    
            %number of measurement points from measurement start to end
            measNrSum = fix((sampleEnd - sampleStart)/measSampleStep);
    
            %initialization 
            %set the local time to ind for the first satellite calculation of
            %reciever positon. After first fix local time will be updated by
            %measurement sample step
            localTime = inf;
    
            fprintf("Positions are being computed. Please wait....\n");
            for currMeasNr = 1:measNrSum
                %positon index of current measurement time in IF signal stream
                currMeasSample = sampleStart + measSampleStep*(currMeasNr - 1);
                %find psuedo-ranges
                [navSol.rawP(:, currMeasNr), transmitTime, localTime] = obj.calcPsuedorange(tracking_results, subFrameStart, TOW, currMeasSample, localTime, trackedChan);
    
    
                %Find satellites positions
                satPositions = obj.satPos(transmitTime(trackedChan), [tracking_results(trackedChan).PRN], ephem);
    
                %find the reciever position, doing this requires at least 4
                %satellites to be present or else the solution is
                %overdetermined
                if(size(transmitTime, 2) > 3)
                    %calculate reciever position and time bias
                    [recv_state, LLA] = obj.leastSquares(satPositions, navSol.rawP(trackedChan, currMeasNr));
    
    
                    navSol.X(currMeasNr)           = recv_state(1);
                    navSol.Y(currMeasNr)           = recv_state(2);
                    navSol.Z(currMeasNr)           = recv_state(3);
                    navSol.bias_t(currMeasNr)     = recv_state(4);
                    navSol.latitude(currMeasNr)    = LLA(1);
                    navSol.longitude(currMeasNr)   = LLA(2);
                    navSol.height(currMeasNr)      = LLA(3);
                else
                    navSol.X(currMeasNr)           = NaN;
                    navSol.Y(currMeasNr)           = NaN;
                    navSol.Z(currMeasNr)           = NaN;
                    navSol.bias_t(currMeasNr)     = NaN;
                    navSol.latitude(currMeasNr)    = NaN;
                    navSol.longitude(currMeasNr)   = NaN;
                    navSol.height(currMeasNr)      = NaN;
    
                    disp(['   Measurement No. ', num2str(currMeasNr), ...
                           ': Not enough information for position solution.']);
                    return;
                end
            end 
        end

        function [state, LLA] = leastSquares(obj, satPos, psuedorange)

            X = zeros(4, 1);


            P = satPos;
            y = psuedorange;
    
            h = zeros(size(y, 1), 1);
            H = zeros(size(y, 1), 4);
    
            while true
                for kk = 1:size(y, 1)
                    v = (X(1:3) - P(:, kk));
                    f = (v'*v)^(1/2);
    
                    h(kk) = f - X(4);
                    er = v' /f;
                    H(kk, :) = [er, 1];
                end
    
                X_old = X;
    
%                 calculate new estimate of X
                X = X + H \ (y-h);
    
                delta_x = X(1:3) - X_old(1:3);
                if((delta_x' * delta_x)^(1/2) < 0.1)
                    break;
                end
   
            end
           
            state = X;
            [lat, lon, alt] = obj.ECEF_to_LLA(X(1:3));
            LLA = [-lat, lon, alt]';
        end

        function [Az, El] = calcAzEl(obj, Xt, Yt, Zt, Rx0)
            Az = zeros(size(Xt));
            El = Az;
            [lat,lon, ~] = obj.ECEF_to_LLA(Rx0);
            for i = 1:size(Xt,1)
                for j = 1:size(Xt,2)
                    r_Rx0_Svij = [Xt(i,j);Yt(i,j);Zt(i,j)]-Rx0;
                    r_ENU = obj.LLA_to_ENU(lat,lon,r_Rx0_Svij);
                    El(i,j) = asind(r_ENU(3)/norm(r_ENU));  % Calculates the elevation angle
                    Az(i,j) = atan2d(r_ENU(1),r_ENU(2));    % Calculates the azimuth angle
                end
            end
        end

        function r_ENU = LLA_to_ENU(obj, lat, lon, r_Gnd_Sky)
            R=[-sind(lon),cosd(lon),0;-sind(lat)*cosd(lon),-sind(lat)*sind(lon),cosd(lat);cosd(lat)*cosd(lon),cosd(lat)*sind(lon),sind(lat)]; % Build the ECEF to ENU rotation matrix
            r_ENU=R*r_Gnd_Sky;    % Express the vector from the receiver to the SV in ENU coordinates 
        end

        function [lat, lon, alt] = ECEF_to_LLA(obj, r)
            e_sq = (2*obj.flat - obj.flat^2);
            e = sqrt(e_sq); %eccentricity of the earth

            x = r(1, :);
            y = r(2, :);
            z = r(3, :);

            p=sqrt(x.^2+y.^2);
            [n,m] = size(x);
            % Initializes Newton-Raphson
            k_new=1/(1-e^2)*ones(n,m);
            err_threshold=0.0001*ones(n,m);
            err=1*ones(n,m);
            % Iterates Newton-Raphson
            while any(err>err_threshold)
                k_old=k_new;
                ci=(p.^2+(1-e^2)*z.^2.*k_old.^2).^(3/2)/(obj.a*e^2);
                k_new=1+(p.^2+(1-e^2)*z.^2.*k_old.^3)./(ci-p.^2);
                err=abs(k_new-k_old);
            end
            k=k_new;

            lon=atan2(y,x); % Calculate longitude
            lat=atan2(z.*k,p);   % Calculate latitude

            Rn = obj.a./sqrt(1-e^2*sin(lat).^2);
            alt = p./cos(lat) - Rn;

            lat = lat * 180/pi;
            lon = lon * 180/pi;
        end

        function [psuedoranges, transmitTime, localTime] = calcPsuedorange(obj, trackResults, subFrameStart, TOW, currMeasSample, localTime, channelList)
            transmitTime = inf(1, obj.numChannels);

            %for all channels in the list 
            for channelNr = channelList
                %find index I_P stream whose integration contains current
                %measurment point location
                for index = 1:length(trackResults(channelNr).absoluteSample)
                    if(trackResults(channelNr).absoluteSample(index) > currMeasSample)
                        break;
                    end
                end
                index = index - 1;

                %update the phase step based on code freq and sampling
                %frequency
                delta_code_phase = trackResults(channelNr).codeFreq(index) / obj.fsamp;

                %code phase from start of PRN code to current measurement
                %sample location
                codePhase = trackResults(channelNr).r_code_phase(index) + delta_code_phase*(currMeasSample - trackResults(channelNr).absoluteSample(index));

                %transmitting time (in units of s) at current measurement
                %sample locartion codephase/obj.codeLength: fraction part
                %of PRN code index-subframeStart(channelNr): integer
                %number of PRN code

                transmitTime(channelNr) = (codePhase/obj.codeLength + index - subFrameStart(channelNr)) * obj.codeLength/obj.fchip + TOW(channelNr);

            end

            %at first time of fix, local time is initialized by transmit
            %time and obj.startoffset

            if(localTime == inf)
                maxTime = max(transmitTime(channelList));
                localTime = maxTime + obj.startOffset/1000;
            end

            %convert the travel time to distance 
            psuedoranges = (localTime - transmitTime)*obj.c;
        end

        function satPositions = satPos(obj, transmit_time, prnList, ephem)
            numOfSat = size(prnList, 2);
            satPositions = zeros(3, numOfSat);

            for ii = 1:numOfSat
                prn = prnList(ii);
                sqrtA = ephem(prn).sqrtA;        %Square root of the Semi-Major Axis
                t_oe = ephem(prn).t_oe;          %Reference Time Ephemeris
                deltan = ephem(prn).deltan;      %Mean Motion Differene from Computed Value
                M_0 = ephem(prn).M_0;            %Mean Anomaly at Reference Time
                w = ephem(prn).omega;            %Argument of Perigree
                C_us = ephem(prn).C_us;          %Amplitude of Sine Harmonic Correction Term to the Argument of Latitude
                C_uc = ephem(prn).C_uc;          %Amplitude of Cosine Harmonic Correction Term to the Argument of Latitude
                C_rs = ephem(prn).C_rs;          %Amplitude of the Sine Harmonic Correction Term to the Orbit Radius
                C_rc = ephem(prn).C_rc;          %Amplitude of the Cosine Harmonic Correction Term to the Orbit Radius
                C_is = ephem(prn).C_is;          %Amplitude of the Sine Harmonic Correction Term to the Angle of Inclination
                C_ic = ephem(prn).C_ic;          %Amplitude of the Sine Harmonic Correction Term to the Angle of Inclination
                i0 = ephem(prn).i_0;             %Inclination Angle at Reference Time
                IDOT = ephem(prn).IDOT;          %Rate Inclination Angle
                omega_0 = ephem(prn).omega_0;    %Longitude of Ascending Node of Orbit Plane at Weekly Epoch
                omega_dot = ephem(prn).omega_dot;%Rate of Right Ascension
                e = ephem(prn).e; %eccentricity
                T_GD = ephem(prn).T_GD;
                a_f2 = ephem(prn).a_f2;
                a_f1 = ephem(prn).a_f1;
                a_f0 = ephem(prn).a_f0;

                dt = obj.weekroll(transmit_time(ii) - t_oe);
                satClkCorr = (a_f2 * dt + a_f1)*dt + a_f0 - T_GD;
                t = transmit_time(ii) - satClkCorr;

                A = sqrtA^2;                    %Semi-major axis
                t_k = obj.weekroll(t - t_oe);
                n_0 = sqrt(obj.mu/(A^3));           %(rad/s) Corrected mean motion
                n = n_0 + deltan; %mean motion
                M_k = M_0+n*t_k;                %Mean anomaly
                M_k = rem(M_k + 2*obj.gpsPi, 2*obj.gpsPi); % reduce mean anomoly between 0 and 360 degrees

                %Kepler's equation (M_k = E_k - esinE_k) may be solved for E_k by iteration:
                E = zeros(1,10);
                E(1) = M_k;                         %Initial value(radians)
                for j = 2:10
                    E(j) = M_k + e*sin(E(j-1));
                end

                E_k = E(10);                         %Final value (radians)
                E_k = rem(E_k + 2*obj.gpsPi, 2*obj.gpsPi);
               
                v_k = atan2(sqrt(1-e^2)*sin(E_k), cos(E_k)-e);        %True Anomaly
                phi_k = v_k + w;                                   %Argument of Latitude
                phi_k = rem(phi_k + 2*obj.gpsPi, 2*obj.gpsPi);

                %%Second Harmonic Perturbations
                duk = C_us*sin(2*phi_k) + C_uc*cos(2*phi_k);       %Argument of Lattitude Correction
                drk = C_rs*sin(2*phi_k) + C_rc*cos(2*phi_k);       %Radius Correction
                dik = C_is*sin(2*phi_k) + C_ic*cos(2*phi_k);       %Inclination Correction
                uk = phi_k+duk;                                    %Corrected Argument of Latitude
                rk = A*(1-e*cos(E_k))+ drk;                        %Corrected radius
                ik = i0 +dik + (IDOT)*t_k;                         %Corrected Inclination
                %%Positions in orbital
                x_k = rk*cos(uk);
                y_k = rk*sin(uk);
                %Corrected longitude of ascending node
                omega_k = omega_0 + (omega_dot - obj.omega_e)*t_k - obj.omega_e*t_oe;
                omega_k = rem(omega_k + 2*obj.gpsPi, 2*obj.gpsPi);
                %Earth-fixed
                satPositions(1, ii) = x_k*cos(omega_k) - y_k*cos(ik)*sin(omega_k);
                satPositions(2, ii) = x_k*sin(omega_k) + y_k*cos(ik)*cos(omega_k);
                satPositions(3, ii) = y_k*sin(ik);
            end
        end

        function t = weekroll(obj, time)
            t = time;
            if time > obj.half_week
                t = time - 2*obj.half_week;
            elseif time < -obj.half_week
                t = time + 2*obj.halfweek;
            end
        end

        function eph = ephem_init(obj)
            %create ephermeris data structure 
            eph.idValid(1:5) = zeros(1, 5); %valid subframe
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
            eph.sqrtA = []; %
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
        end

        function [eph, subframeStart, TOW] = decode(obj, I_P)     
            eph = obj.ephem_init();
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

            %copy 5-subframes long

            navBitsSamples = I_P(subframeStart - 20 : subframeStart + (1500 * 20) -1)';
            %reshape into 20 different values of bits
            navBitsSamples = reshape(navBitsSamples, 20, size(navBitsSamples, 1) / 20);
            navBits = sum(navBitsSamples);

            navBits = (navBits > 0);

            navBitsBin = dec2bin(navBits);


            [eph, TOW] = obj.ephermeris(navBitsBin(2:1501)', navBitsBin(1));
        end

        function [eph, TOW] = ephermeris(obj, bits, D30star) 

            eph = obj.ephem_init();
            %decode all five subframes
            for i = 1:5
                %grab a single subframe out of all of the elements
                %--- "Cut" one sub-frame's bits ---------------------------------------
                subframe = bits(300*(i-1)+1 : 300*i);

                %correct polarity of data bits in all 10 words
                for j = 1:10
                    [subframe(30*(j-1)+1 : 30*j)] = obj.checkPhase(subframe(30*(j-1)+1 : 30*j), D30star);
                    
                    D30star = subframe(30*j);
                end

                %decode subframe id
                subframeID = bin2dec(subframe(50:52));
                
                switch subframeID
                    case 1
                        eph.weekNumber = bin2dec(subframe(61:70)) + 1024;
                        eph.accuracy = bin2dec(subframe(73:76));
                        eph.health = bin2dec(subframe(77:82));
                        eph.T_GD = obj.twosComp2dec(subframe(197:204)) * 2^(-31);
                        eph.IODC = bin2dec([subframe(83:84) subframe(197:204)]);
                        eph.t_oc = bin2dec(subframe(219:234)) * 2^(4);
                        eph.a_f2 = obj.twosComp2dec(subframe(241:248)) * 2^(-55);
                        eph.a_f1 = obj.twosComp2dec(subframe(249:264)) * 2^(-43);
                        eph.a_f0 = obj.twosComp2dec(subframe(271:292)) * 2^(-31);
                        eph.idValid(1) = 1;
                    case 2
                        eph.IODE_2 = bin2dec(subframe(61:68));
                        eph.C_rs = obj.twosComp2dec(subframe(69:84)) * 2^(-5);
                        eph.deltan = obj.twosComp2dec(subframe(91:106)) * 2^(-43) * obj.gpsPi;
                        eph.M_0 = obj.twosComp2dec([subframe(107:114) subframe(121:144)]) * 2^(-31) * obj.gpsPi;
                        eph.C_uc = obj.twosComp2dec(subframe(151:166)) * 2^(-29);
                        eph.e = bin2dec([subframe(167:174) subframe(181:204)]) * 2^(-33);
                        eph.C_us = obj.twosComp2dec(subframe(211:226)) * 2^(-29);
                        eph.sqrtA = bin2dec([subframe(227:234) subframe(241:264)]) * 2^(-19);
                        eph.t_oe = bin2dec(subframe(271:286)) * 2^(4);
                        eph.idValid(1) = 2;
                    case 3
                        eph.C_ic = obj.twosComp2dec(subframe(61:76)) * 2^(-29);
                        eph.omega_0 = obj.twosComp2dec([subframe(77:84) subframe(91:114)]) * 2^(-31) * obj.gpsPi;
                        eph.C_is = obj.twosComp2dec(subframe(121:136)) * 2^(-29);
                        eph.i_0 = obj.twosComp2dec([subframe(137:144) subframe(151:174)]) * 2^(-31) * obj.gpsPi;
                        eph.C_rc = obj.twosComp2dec(subframe(181:196)) * 2^(-5);
                        eph.omega = obj.twosComp2dec([subframe(197:204) subframe(211:234)]) * 2^(-31) * obj.gpsPi;
                        eph.omega_dot = obj.twosComp2dec(subframe(241:264)) * 2^(-43) * obj.gpsPi;
                        eph.IODE_3 = bin2dec(subframe(271:278));
                        eph.IDOT = obj.twosComp2dec(subframe(279:292)) * 2^(-43) * obj.gpsPi;
                        eph.idValid(1) = 3;
                    case 4
                        %Almanac amd ionospheric model, UTC parameters
                    case 5
                        %SV almanac and health
                end
            end
            TOW = bin2dec(subframe(31:47)) * 6 - 30;
            eph.TOW = TOW;
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
                    channel = obj.rankAcquisition(acq_results);
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

                fclose(fid);

                %% Position solution
                [navResults, ephem] = obj.compute_position(track_results);
                kml_str = ge_point(navResults.longitude + 2, -navResults.latitude - .48, navResults.height);
                ge_output(['./NavTek.kml'], kml_str);

            end
        end
    end
end