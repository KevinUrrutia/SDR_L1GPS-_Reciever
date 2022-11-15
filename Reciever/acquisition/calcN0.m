function [N0, SNk2] = calcN0(fdmax, fdmin, deltafd, Ns, f_samp, C, Xp)
    %input -> fdmax -> maximum doppler
    %         fdmin -> minimum doppler
    %         deltafd -> difference between doppler values
    %         Ns -> code length
    %         f_samp -> sampling frequency
    %         C     -> oversampled PRN
    %         Xp actual data
    %
    %output -> N0 -> noise density
    %          SNk2-> power spectra for each satellite

    SNk2 = zeros(1, 605121);
    i = 0;
    for fd = fdmin:deltafd:fdmax
        X_hat = cos(2*pi*(1:Ns)* fd * (1 / f_samp)) - 1j * sin(2*pi*(1:Ns)* fd *(1/f_samp)); 
        i = i + 1;
        for ts = 0:length(C) - 1
            C_hat = circshift(C', ts);
            [R, ~] = ccorr(C_hat, X_hat);
            i = i + 1;
            SNk2(i) = norm(R'*Xp')^2;
        end
    end

    N0 = mean(mean(SNk2));
end
