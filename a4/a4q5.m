clear;clc;
% QR Factorization Method
% 1  -  MATLAB Built-in
% 2  -  Classical Gram-Schmidt
% 3  -  Modified Gram-Schmidt
method = 0;
% Matrix Sizes
m = 30;
n = 20;

%% Vandermonde Initialization
A = vanderM(m,n);
if method == 0
    %% MATLAB QR
    [Q,R] = qr(A,0);
elseif method == 1
    %% CGS
    [Q,R] = cgs(A);
elseif method == 2
    %% MGS
    [Q,R] = mgs(A);
end

%% Output
norm(A - Q*R,2)/norm(A,2)
norm(Q'*Q - eye(n), 2)


