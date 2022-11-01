%Generate a 1/0 valued linear feddback shift register (LFSR) sequence
%
%INPUT
% 
%n: Number of stages in the linear feedbackshift register
%
%ciVec: Nc by 1 vector whose elements give the indices of the Nc nonzero
%       connection elements. For example, if the characteristic polynomial of LFSR
%       is f(D)=1 + D^2 + D^3, then ciVec = [2,3]' or [3, 2]
%
%a0Vec: n by 1 1/0 valued initial state of the LFSR, where a0Vec = [a(-1),
%       a(-2),.....,a(-n)]'. In defining the initial LFSR staten a Fibonacci LFSR
%       sequence is assumed
%
%OUTPUT
%
%lsfrSeq: m by 1 vector whose elements are the 1/0 valued LFSR sequence
%         corresponding to n, ciVec, and a0Vec, where m = 2n^2 - 1. If the sequence
%         is a maximal length sequence, then ther is no repition in the sequence
%         elements

function [lsfrSeq] = generateLFSRSeq(n, ciVec, a0Vec)
    %the taps are defined by ciVec, ciVec is the index location for the
    %taps in the m length sequence

    %create a tap vector length n 
    ciVec = unique(ciVec);
    c = zeros(n, 1);
    c(ciVec) = 1;

    m = 2^n - 1; 

    lsfrSeq = zeros(m, 1);

    %the initial state of the register is defined by a0Vec

    for i = 1 : m
        temp = mod(c'*a0Vec, 2);
        lsfrSeq(i) = a0Vec(end);
        a0Vec = circshift(a0Vec, 1);
        a0Vec(1) = temp;
    end

end