function [ U, S, V ] = SVDKog( A )
% Return the SVD of matrix A [ m x n ]
% U^T S V = A
% U is an orthogonal matrix [ m x m ]
% V is an orthogonal matrix [ n x n ]
% S is a diagonal matrix [ m x n ]

tol = 1e-5;

m = size(A, 1);
n = size(A, 2);

U = eye(m);
V = eye(n);
S = A;

eps = tol * norm(A, 'fro');

while( sum(sum(abs(S))) - trace(abs( S(1 : n, 1 : n) )) > tol )
    for i = 1 : n
        for j = i + 1 : n
            [ c1, s1, c2, s2, sig1, sig2 ] = ...
                jacrot( [ S(i, i), S(i, j); S(j, i), S(j, j) ]);
            
            sri = c1 * S(i, :) - s1 * S(j, :);
            srj = s1 * S(i, :) + c1 * S(j, :);
            
            uci = c1 * U(:, i) - s1 * U(:, j);
            ucj = s1 * U(:, i) + c1 * U(:, j);
            
            S(i, :) = sri;
               S(j, :) = srj;
            
            U(:, i) = uci;
            U(:, j) = ucj;
            
            sci = c2 * S(:, i) - s2 * S(:, j);
            scj = s2 * S(:, i) + c2 * S(:, j);
                        
            vci = c2 * V(:, i) - s2 * V(:, j);
            vcj = s2 * V(:, i) + c2 * V(:, j);          
            
            S(:, i) = sci;
            S(:, j) = scj;
            
            V(:, i) = vci;
            V(:, j) = vcj;
        end
        
        for j = n + 1 : m
            [ c1, s1, c2, s2, sig1, sig2 ] = ...
                jacrot( [ S(i, i), 0 ; S(j, i), 0 ]);
            sri = c1 * S(i, :) - s1 * S(j, :);
            srj = s1 * S(i, :) + c1 * S(j, :);
            
            uci = c1 * U(:, i) - s1 * U(:, j);
            ucj = s1 * U(:, i) + c1 * U(:, j);
            
            S(i, :) = sri;
            S(j, :) = srj;
            
            U(:, i) = uci;
            U(:, j) = ucj;
        end
    end            
end



end
