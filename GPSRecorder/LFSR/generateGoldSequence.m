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
%
%
%+==============================================================================+
function [GoldSeq] = generateGoldSequence(n,ciVec1,ciVec2,a0Vec1,a0Vec2)

    ciVec1 = unique(ciVec1);
    ciVec2 = unique(ciVec2);
    c1 = zeros(n,1);
    c2 = zeros(n,1);
    c1(ciVec1) = 1;
    c2(ciVec2) = 1;
    a1 = a0Vec1;
    a2 = a0Vec2;
    for m = 1:2^n -1
        temp1 = mod(c1'*a1,2);
        temp2 = mod(c2'*a2,2);
        
        GoldSeq(m) = mod(a1(end),a2(end));
        a1 = circshift(a1,1);
        a2 = circshift(a2,1);
        a1(1) = temp1;
        a2(1) = temp2;
    end
end

