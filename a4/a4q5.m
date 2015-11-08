clear;
method = 0;
% Matrix Sizes
m = 10;
n = 5;

%% Matrix Initialization
A = zeros(m,n);
for i = 1 : m
    for j = 1 : n
        A(i,j) = (j/n)^(i-1);
    end
end

if method == 0
    %% MATLAB QR
    [Q,R] = qr(A,0);
elseif method == 1
    %% CGS
    for k = 1 : n
        for i = 1 : k - 1
            R(i,k) = Q(:,i)'*A(:,k);
        end
        for i = 1 : k - 1
            A(:,k) = A(:,k) - R(i,k)*Q(:,i);
        end
        R(k,k) = norm(A(:,k),2);
        Q(:,k) = A(:,k)/R(k,k);
    end
elseif method == 2
    %% MGS
    for k = 1 : n
        for i = 1 : k - 1
            R(i,k) = Q(:,i)'*A(:,k);
            A(:,k) = A(:,k) - R(i,k)*Q(:,i);
        end
        R(k,k) = norm(A(:,k),2);
        Q(:,k) = A(:,k)/R(k,k);
    end
end
%% Output
Q
R
norm(A - Q*R,2)/norm(A,2)
norm(Q'*Q - eye(n), 2)


