%[lfsrSeq] = generateLfsrSequence(n,ciVec,a0Vec)
%
% Generate a 1/0-valued linear feedback shift register (LFSR) sequence.
%
% INPUTS
%
% n ------ Number of stages in the linear feedback shift register.
%
% ciVec -- Nc-by-1 vector whose elements give the indices of the Nc nonzero
% connection elements. For example, if the characteristic polynomial
% of an LFSR is f(D) = 1 + D^2 + D^3, then ciVec = [2,3]’ or [3,2]’.
%
% a0Vec -- n-by-1 1/0-valued initial state of the LFSR, where a0Vec = [a(-1),
% a(-2), ..., a(-n)]’. In defining the initial LFSR state, a
% Fibonacci LFSR implementation is assumed.
%
% OUTPUTS
%
% lfsrSeq -- m-by-1 vector whose elements are the 1/0-valued LFSR sequence
% corresponding to n, ciVec, and a0Vec, where m = 2^n - 1. If the
% sequence is a maximal-length sequence, then there is no
% repetition in the m sequence elements.
%
%+------------------------------------------------------------------------------+
% References:
%
%+==============================================================================+
function [lfsrSeq] = generateLfsrSequence(n,ciVec,a0Vec)

    ciVec = unique(ciVec);
    c = zeros(n,1);
    c(ciVec) = 1;
    a = a0Vec;
    for m = 1:2^n -1
        temp = mod(c'*a,2);
        lfsrSeq(m) = a(end);
        a = circshift(a,1);
        a(1) = temp;
    end
end

