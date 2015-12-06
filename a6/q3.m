%% Initialize
clear;clc;
icase = 5;
switch icase
    case 1
        A = [5 4 1 1;
             4 5 1 1;
             1 1 4 2;
             1 1 2 4];
    case 2
        A = [6 4 4 1;
             4 6 1 4;
             4 1 6 4;
             1 4 4 6];
    case 3
        A = [6 -3 4 1;
             4  2 4 0;
             4 -2 3 1;
             4  2 3 1];
    case 4
        A = [4 -5  0  3;
             0  4 -3 -5;
             5 -3  4  0;
             3  0  5  4];
    case 5
        A = [10 -19 17 -12 4 1;
              9 -18 17 -12 4 1;
              8 -16 15 -11 4 1;
              6 -12 12 -10 4 1;
              4  -8  8  -6 1 2;
              2  -4  4  -3 1 0];
end
n = size(A,1);

%% Iterate
maxIt = 100000;
tol = 1e-6;
Ak = A
for it = 1 : maxIt
    [Qk, Rk] = qr(Ak);
    Ak = Rk*Qk;
    lowTsum = sum(sum(abs(tril(Ak,-1))))
    if lowTsum < tol;
        break
    end
end

Ak

latexMat = sprintf(['begin{bmatrix}[r] \n']);
for i = 1 : n
    for j = 1 : n
        latexEnt = sprintf('%f & ', Ak(i,j));
        latexMat = sprintf([latexMat, latexEnt]);
    end
    latexMat = sprintf([latexMat,   '\n']);
end
latexMat = sprintf([latexMat,'\\end{bmatrix} \n']);

latexMat