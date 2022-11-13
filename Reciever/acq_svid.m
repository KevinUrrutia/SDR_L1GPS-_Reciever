function [SK2] = acq_svid(fdmax, fdmin, Ns, deltafd, f_samp, C)
    %input -> fdmax -> maximum doppler
    %         fdmin -> minimum doppler
    %         deltafd -> difference between doppler values
    %         Ns -> code length
    %         f_samp -> sampling frequency
    %         C     -> oversampled PRN
    %
    %output-> SK2-> power spectra for each satellite

   SK2 = cell([1, 605121]);
   
   for i = [1, 2, 5, 10, 12, 15, 29, 30]
        j = 0;
        for fd = fdmin:deltafd:fdmax
            X_hat = cos(2*pi*(1:Ns) * fd * (1 / f_samp)) - 1j * sin(2*pi*(1:Ns)* fd * (1/f_samp)); 
            j = j + 1;
            for ts = 0:length(C{i}) - 1
                C_hat = circshift(C{i}', ts);
                R = X_hat * C_hat;
                j = j + 1
                SK2{i}(ts + 1, j) = norm(R'*X_hat')^2;
            end
        end
    end
end