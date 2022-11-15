function [SK2] = acq_svid(fdmax, fdmin, Ns, deltafd, f_samp, C, Xp)
    %input -> fdmax -> maximum doppler
    %         fdmin -> minimum doppler
    %         deltafd -> difference between doppler values
    %         Ns -> code length
    %         f_samp -> sampling frequency
    %         C     -> oversampled PRN
    %         Xp     actual data
    %
    %output-> SK2-> power spectra for each satellite

   SK2 = cell([1, 30]);
   
   for i = [1 2 5 10 12 15 29 30]
        %code start time
        cnt = 0;
        for fd = fdmin:deltafd:fdmax 
            Xp_c = cos(2*pi*fd*(1:Ns)*(1/f_samp)) - 1j*sin(2*pi*fd*(1:Ns)*(1/f_samp));
            cnt = cnt +1;
            for ts = 0:length(C{i})-1
                C_hat = circshift(C{i}',ts);
                Xp_hat = Xp_c.*C_hat';
                SK2{i}(ts+1,cnt) = norm(Xp_hat*Xp')^2;
            end
        end 
   end
end