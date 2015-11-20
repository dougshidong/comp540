n = 9;

t = 1.11 * ones(n,1);
b = 2.11 * ones(n,1);
for i = 1 : n
    t(i) = t(i) + (i-1)*1e-3;
    b(i) = b(i) + (i-1)*1e-3;
end

b(1:3) = b(1:3) + 0.0001;

A = [ones(n,1), t];

A = single(A);
b = single(b);

ka = cond(A)

xqr = A\b

res_qr = norm(b - A*xqr)

kaat = cond(A'*A)

xne = (A'*A)\(A'*b)

res_ne = norm(b - A*xne)

A = double(A);
b = double(b);

xdb = A\b

n = 9;

t = 1.11 * ones(n,1);
b = 2.11 * ones(n,1);
for i = 1 : n
    t(i) = (t(i) + (i-1)*1e-3 - 1.11)*100;
    b(i) = b(i) + (i-1)*1e-3;
end

b(1:3) = b(1:3) + 0.0001;

A = [ones(n,1), t];

A = single(A);
b = single(b);

ka = cond(A)

xqr = A\b

res_qr = norm(b - A*xqr)

kaat = cond(A'*A)

xne = (A'*A)\(A'*b)

res_ne = norm(b - A*xne)

A = double(A);
b = double(b);

xdb = A\b