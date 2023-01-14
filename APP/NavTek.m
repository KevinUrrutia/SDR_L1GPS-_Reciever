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
        fchip = 1.203e6; %[Hz]
        %Code Length
        codeLength = 1023;

        %% Acquisition Settings
        %satellite PRN numbers
        SVn = 1:32;   
        %minimum, maximum, resolution of doppler search
        fdmin = -7000; %[Hz]
        fdmax = 7000; %[Hz]
        delta_fd = 100; %[Hz]
        %NonCoherent integration times after 1ms coherent integration
        NonCohTime = 20; 
        %Sampling rate for downsampling 
        resamplingThresh = 8e6; %[Hz]
        %acquisition threshold for desision rule
        acqThresh = 3.5;

        %% Constants
        c = 299792458;
    end

    methods
        function obj = NavTek()
%             prompt = "Enter full IQ file path\n";
%             file_loc = input(prompt, 's');
            obj.fnameIQ = gui.iq_file;
        end
        
        function [data, sampsPerCode, sigspec, freqv] = probeData(obj)
            %% Generate and plot raw data
            fid = fopen(obj.fnameIQ, 'rb');

            if (fid > 0)
                %move to the start of the where we want to start post
                %processing
                fseek(fid, obj.skipNumBytes, 'bof');

                %find the number of samples per spreading code
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
%                 subplot(3, 2, 4);
%                 plot(1000 * timeScale(1:round(sampsPerCode/2)), ...
%                 real(data(1:round(sampsPerCode/2))));
                % return data

                axis tight;    grid on;
                title ('Time domain plot (I)');
                xlabel('Time (ms)'); ylabel('Amplitude');

                subplot(3, 2, 3);
%                 plot(1000 * timeScale(1:round(sampsPerCode/2)), ...
%                 imag(data(1:round(sampsPerCode/2))));
                % return sampsPerCode, data

                axis tight;    grid on;
                title ('Time domain plot (Q)');
                xlabel('Time (ms)'); ylabel('Amplitude');

                %% Frequency domain plot
                subplot(3,2,1:2);
                [sigspec,freqv]=pwelch(data, 32768, 2048, 32768, obj.fsamp,'twosided');
                plot(([-(freqv(length(freqv)/2:-1:1));freqv(1:length(freqv)/2)])/1e6, ...
                10*log10([sigspec(length(freqv)/2+1:end);
                sigspec(1:length(freqv)/2)]));
                % return sigspec, freqv
                axis tight;
                grid on;
                title ('Frequency domain plot');
                xlabel('Frequency (MHz)'); ylabel('Magnitude');

                %% Histogram
                subplot(3, 2, 6);
                histogram(real(data), -128:128)
                dmax = max(abs(data)) + 1;
                % return data
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
            end
        end

        function postProcess(obj)
            obj.probeData();
            % Acquisition Function
        end
    end
end