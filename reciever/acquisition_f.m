function acquisition_f(fdVec, Xp, C, Dopplerm, timem, Ns, f_samp, fIF, SVn, Tc, startTimeVec)
    %inputs: 
       %fdVec -> possible doppler values to search
       %Xp -> complex signal
       %C -> C/A codes
       %Dopplerm -> linspace of doppler values
       %timem -> linspace of time values
       %SVn -> satellite numbers to search
       %Ns -> number of samples in code period
       %f_samp -> sampling rate
       %Tc -> chipping periond
       %startTimeVec ->plot start time

       %calculate Noise Density
       N0 = 1.291271654429430e4;

       %---PRN Search
       SK2 = acq_svid(fdVec, Ns, f_samp, C, Xp, fIF, SVn);

       %--Plot PRN
       [ts_hat, fd_hat, C_N0] = plotPRN(SK2, Tc, f_samp, N0, startTimeVec, fdVec, Dopplerm, timem, SVn);
      
end