%%
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
%%
% SVD Decomposition Using Kogbetliantz Algorithm
[ Ub Sb Vb ] = SVDKog(B);
%%
% Test for Algorithm Stability
C = B + 0.00001;

[ Uc Sc Vc ] = SVDKog(C);

sprintf('Test for Stability')
norm(Ub - Uc)
norm(Sb - Sc)
norm(Vb - Vc)
%%
% SVD Decomposition Using Built-In MATLAB SVD
A = B;
A(1, 1) = -1;

[ Um Sm Vm ] = svd(A)

% Invert Diagonal Matrix
Sinv2 = zeros(n,m);
Sinv3 = zeros(n,m);
invd = 1./diag(Sm);
Sinv2(1:n-1,1:n-1) = diag(invd(1:n-1));
Sinv3(1:n,1:n) = diag(invd(1:n));

% Form Moore-Penrose Inverse Matrix
sprintf('Moore-Penrose Inverse')
G2 = Vm*Sinv2*Um'
G3 = Vm*Sinv3*Um'

norm(A*G2*A - A)
norm(G2*A*G2 - G2)
norm((A*G2)' - A*G2)
norm((G2*A)' - G2*A)


%%
% Solve Least Square Problem
b = [5:-1:1]';
xls2 = G2*b
xls3 = G3*b