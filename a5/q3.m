clear;
clc;
% Initialize B Matrix
B = [ -1.1,  -6, -11;
         2,   7,  12;
        -3,  -8, -13;
         4,   9,  14;
        -5, -10, -15 ];
    
m = size(B, 1);
n = size(B, 2);

% SVD Decomposition Using Kogbetliantz Algorithm
[ Ub Sb Vb ] = SVDKog(B);

% Test for Algorithm Stability
C = B + 0.00001;

[ Uc Sc Vc ] = SVDKog(C);

sprintf('Test for Stability')
Ub - Uc
Sb - Sc
Vb - Vc

% SVD Decomposition Using Built-In MATLAB SVD
A = B;
A(1, 1) = -1;

[ Um Sm Vm ] = svd(A);

% Invert Diagonal Matrix
Sinv = zeros(n,m);
invd = 1./diag(Sm);
Sinv(1:n,1:n) = diag(invd);

% Form Moore-Penrose Inverse Matrix
sprintf('Moore-Penrose Inverse')
G = Vm*Sinv*Um'

A*G*A - A
G*A*G - G
(A*G)' - A*G
(G*A)' - G*A
% BAD because A does not have full column rank

% Need to exclude a column of A
sprintf('Moore-Penrose Inverse with less column')
% SVD Decomposition Using Built-In MATLAB SVD
D = B(1:5,1:2);
D(1, 1) = -1;

[ Um Sm Vm ] = svd(D);

% Invert Diagonal Matrix
Sinv = zeros(n-1,m);
invd = 1./diag(Sm);
Sinv(1:n-1,1:n-1) = diag(invd)

% Form Moore-Penrose Inverse Matrix
G = Vm*Sinv*Um'

D*G*D - D
G*D*G - G
(D*G)' - D*G
(G*D)' - G*D
% BAD because A does not have full column rank
% Need to exclude a column of A


% Solve Least Square Problem
