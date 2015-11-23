%% a1(t) = 1 and a2(t) = t
clc;clear;
n = 9;

% Initialization
t = 1.11 * ones(n,1);
b = 2.11 * ones(n,1);
for i = 1 : n
    t(i) = t(i) + (i-1)*1e-3;
    b(i) = b(i) + (i-1)*1e-3;
end

b(1:3) = b(1:3) + 0.0001;

A = [ones(n,1), t]

% Single Precision QR
A = single(A);          
b = single(b);

ka = cond(A)

xqr = A\b

res_qr = norm(b - A*xqr)

% Single Precision Normal Equations
kaat = cond(A'*A)

xne = (A'*A)\(A'*b)

res_ne = norm(b - A*xne)

% Double Precision QR
A = double(A);
b = double(b);
kab = cond(A)
xdb = A\b
res_db = norm(b - A*xdb)

%% Repeat with a1(t) = 1, a2(t) = (t-1.11)*100
clear;clc;
n = 9;

t = 1.11 * ones(n,1);
b = 2.11 * ones(n,1);
for i = 1 : n
    t(i) = (t(i) + (i-1)*1e-3 - 1.11)*100;
    b(i) = b(i) + (i-1)*1e-3;
end

b(1:3) = b(1:3) + 0.0001;

A = [ones(n,1), t];

% Single Precision QR
A = single(A)
b = single(b)

ka = cond(A)

xqr = A\b

res_qr = norm(b - A*xqr)

% Single Precision Normal Equations
kaat = cond(A'*A)

xne = (A'*A)\(A'*b)

res_ne = norm(b - A*xne)

% Double Precision QR
A = double(A);
b = double(b);
kab = cond(A)
xdb = A\b
res_db = norm(b - A*xdb)