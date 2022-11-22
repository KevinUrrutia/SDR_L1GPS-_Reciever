function [ts_hat, fd_hat, C_N0] = plotPRN(SK2, Tc, fsamp, N0, startTimeVec, fdVec, Dopplerm, timem, SVn)
    %inputs: SK2: Power spectra of each satellite
    %        Tc: chipping period
    %        fsamp: sampling rate
    %        N0: noise desity

    %outputs: ts_hat: code estimate
    %         fd_hat: doppler est
    %         C_N0: carrier to noise ratio
 
    ts_hat = zeros(1, 8);
    fd_hat = zeros(1, 8);
    C_N0 = zeros(1, 8);

    for i = SVn
         [maxSk2Vec,maxSk2indVec] = max(SK2{i});
         [maxSk2, ind_fd_hat] = max(maxSk2Vec);
         ind_ts_hat = maxSk2indVec(ind_fd_hat);
         ts_hat(i) = mod(startTimeVec(ind_ts_hat),floor(Tc*fsamp))/fsamp;
         fd_hat(i) = fdVec(ind_fd_hat);
         C_N0(i) = 10*log10((maxSk2-N0)/(N0*Tc));
         fig = surf(Dopplerm,timem,SK2{i});
         grid on;
         xlabel('f_d [Hz]')
         ylabel('t_s [sec]')
         title(strcat(['|S_k|^2 for PRN ' int2str(i)]))
         pause(1);
         clf(fig); 
    end
end