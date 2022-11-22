function [CA] = generateGoldSeq(svid)
    load G2_taps.mat;

    n_stages = 10;

    %init registers
    G1 = ones(1, n_stages);
    G2 = ones(1, n_stages);
    
    %pre-allocate
    CA = zeros(1, 1023);

    %create sequence
    for i = 1:1023
        [g1, G1] = shift(G1, [3, 10], [10]);
        [g2, G2] = shift(G2, [2, 3, 6, 8, 9, 10], G2_taps{svid});

        CA(i) = mod((g1 + g2), 2);
    end
  CA(CA == 0) = -1;
end