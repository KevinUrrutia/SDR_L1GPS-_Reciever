clear; close all; clc;

% ========= Initial Upload of IQ file ================
gui = logpose_app; % Declared class for logpose application
gui.TabGroup.SelectedTab = gui.HOMETab;
pause(10)
% % uiwait(gui.run_IQ_button)
% uploaded = 0;
% while uploaded == 0
%     disp('HELLO')
%     if gui.uploaded == 1
%         uploaded = 1;
%     end
% end 

nav = NavTek(gui);
% 
% % ========= Updates IQ Characteristics plot ==========
% % nav.postProcess(gui);
% 
% gui.TabGroup.SelectedTab = gui.IQTab; % Automatically opens IQ Tab
% 
% % ========= TEST PLOT =================================
% x = 0:pi/100:2*pi;
% y = sin(x);
% plot(gui.iq_frq,x,y)

% ================ PLOT IQ Characteristics ===================
[data, sampsPerCode, sigspec, freqv] = nav.probeData();

 %% Time domain plot
 plot(gui.time_Q,1000 * timeScale(1:round(sampsPerCode/2)), ...
                real(data(1:round(sampsPerCode/2))));
 plot(gui.time_I, 1000 * timeScale(1:round(sampsPerCode/2)), ...
                imag(data(1:round(sampsPerCode/2))));

 %% Frequency domain plot
plot(gui.iq_frq,([-(freqv(length(freqv)/2:-1:1));freqv(1:length(freqv)/2)])/1e6, ...
10*log10([sigspec(length(freqv)/2+1:end);
sigspec(1:length(freqv)/2)]));

 %% Histogram

histogram(gui.hist_I, real(data), -128:128)
dmax = max(abs(data)) + 1;
% return data
axis tight;     adata = axis;
axis([-dmax dmax adata(3) adata(4)]);

histogram(gui.hist_Q, imag(data), -128:128)
dmax = max(abs(data)) + 1;
axis tight;     adata = axis;
axis([-dmax dmax adata(3) adata(4)]);


disp('End of script')

