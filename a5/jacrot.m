function [ c1, s1, c2, s2, sig1, sig2 ] = jacrot( A )
% Given  2 x 2 matrix A, returns the SVD U^TAV = S
% U = [ c1 s1 ;
%      -s1 c1 ]
% V = [ c2 s2 ;
%      -s2 c2]
% S = [ sig1    0 ;
%          0 sig2 ]

w = A(1,1); x = A(1,2);
y = A(2,1); z = A(2,2);

% Symmetrize
if(w + z == 0)
    c = 0;
    s = 1;
else
    t = ( x - y ) / ( w + z );
    c = 1 / sqrt(1 + t * t);
    s = t * c;
end

wt = c * w - s * y;
xt = c * x - s * z;
yt = s * w + c * y;
zt = s * x + c * z;

if(yt == 0)
    c2 = 1;
    s2 = 0;
else
    rho = (zt - wt) / (2 * yt);
    t2 = sign(rho) / ( abs(rho) + sqrt(rho * rho + 1) );
    c2 = 1 / sqrt(1 + t2 * t2);
    s2 = t2 * c2;
end

c1 = c2 * c - s2 * s;
s1 = s2 * c + c2 * s;

sig1 = wt * c2 * c2 - 2 * yt * c2 * s2 + zt * s2 * s2;
sig2 = wt * s2 * s2 + 2 * yt * c2 * s2 + zt * c2 * c2;
    


