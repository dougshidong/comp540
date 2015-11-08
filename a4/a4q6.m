clear;clc;
m = 5;
n = 4;
A = [     7     4     2    10;
     3     5     3     7;
     5     4     8     9;
     9     8     4     5;
     2    10    10     7]

R = A;

for k = 1 : n
    x = R(k:m,k);
    e = zeros(size(x,1),1);
    e(1) = 1;
    rho = sign(x(1))*norm(x);
    v = x + rho*e;
    v = v./norm(v);
    
    tau = -1/(rho*v(1));

    R(k:m,k:n) = R(k:m,k:n) - 2*v*(v'*R(k:m,k:n));
    
    if(k < n-1)
        x = R(k,k+1:n);
        x = x'
        e = zeros(size(x,1),1);
        e(1) = 1;
        rho = sign(x(1))*norm(x);
        v = x + rho*e;
        v = v./norm(v);

        tau = -1/(rho*v(1));

        R(k:m,k+1:n) = R(k:m,k+1:n) - 2*(R(k:m,k+1:n)*v)*v';
    end
    
end

R

[Q,R] = qr(A,0)