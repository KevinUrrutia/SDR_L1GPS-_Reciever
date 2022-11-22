function [X] = iq2if(I, Q, Tl, fIF)
    I_os = zeros(length(I) * 2, 1);
    I_os(1:2:end) = I;

    Q_os = zeros(length(Q) * 2, 1);
    Q_os(1:2:end) = Q;

    T = Tl / 2;

    X = zeros(1, length(I_os));

    for n = 1:length(I_os)
        X(n) = I_os(n) * cos(2*pi*fIF*n*T) - Q_os(n) * sin(2*pi*fIF*n*T);
    end
end