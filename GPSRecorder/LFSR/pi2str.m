function [s,v] = pi2str(m)

%PI2STR Calculate decimals of pi with Machin's formula.
%   [S, V] = PI2STR(M) gives pi truncated to M decimals for M < 2^(53/2).
%   Machin's formula, pi/4 = 4*acot(5) - acot(239), and a number
%   system with base 10^14 is used.
%
%   S is a string output and V is a vector output.
%
%   WARNING! The execution time is some minutes for PI2STR(100000)
%   and some hours for PI2STR(1000000).

%   Author: jonas.lundgren@saabgroup.com, 2008.


nb = 14;                                            % 1 < nb < 15
base = 10^nb;
n = ceil(m/nb);

% Machin's formula
x = 16*xacot(5,n+1,base) - 4*xacot(239,n+1,base);

% Truncate
carry = floor(x(n+2)/base);
x(n+2) = [];

% Canonical form (take care of overflow/underflow)
for k = n+1:-1:2
    xk = x(k) + carry;
    carry = floor(xk/base);
    x(k) = xk - carry*base;                         % x(k) = mod(xk,base)
end
x(1) = x(1) + carry;

v = (int2digits(x(2:end)))';
v = [x(1);v(:)];

% Output string
s = int2str(x(2:end))';
s(isspace(s)) = '0';
s = s(:)';
s = [int2str(x(1)), '.', s(1:m)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = xacot(a,n,base)

% Calculate ACOT(A) to N digits in base BASE for an integer A

b = a^2;
m = floor(n*log(base)/log(b));

x = [1; zeros(n,1)];                                % x = 1
x = xdiv(x,a,base);                                 % x = x/a

y = 0;
for k = 2*m+1:-2:1
	y = xdiv(x,k,base) - xdiv(y,b,base);            % y = x/k - y/b
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = xdiv(x,a,base)

% Calculate X/A in base BASE for an integer A

if a*base < 8e15
    r = 0;
    for k = 1:numel(x)
        d = x(k) + r*base;
        x(k) = floor(d/a);
        r = d - a*x(k);                             % r = mod(d,a)
    end
else
    r = 0;
    b = floor(base/a);
    c = base - a*b;
    for k = 1:numel(x)
        d = x(k) + r*c;
        e = floor(d/a);
        x(k) = e + r*b;
        r = d - a*e;
    end
end



function d = int2digits(i)
%INT2DIGITS Split an integer into digits.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-10-13 14:30:59 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   i = i(:);

   if any(i < 0)
      error('All input values must be non-negative.');
   end

   p = prevpowof10(i);
   p = max(p);
   p = max(p, 0);
   d = rem(floor(i*10.^(-p:0)), 10);


function p = prevpowof10(x)
%PREVPOWOF10 Previous power of 10.
%
%   P = PREVPOWOF10(X) returns the largest integer P such that 10^P <= abs(X).
%
%   Essentially, PREVPOWOF10(X) is the same as FLOOR(LOG(ABS(X)) / LOG(10)),
%   but special care is taken to catch round-off errors.
%
%   See also NEXTPOWOF10, PREVPOW.

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-11-17 11:38:03 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 1, nargin));

   if ~isreal(x)
      error('Input must be real.');
   end

   x = abs(x);
   p = floor(log(x) / log(10));         % estimate
   k = x >= 10.^(p + 1);
   p(k) = p(k) + 1;                     % correction

