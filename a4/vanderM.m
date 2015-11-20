function [A] = vanderM(m,n)

A = zeros(m,n);
for i = 1 : m
    for j = 1 : n
        A(i,j) = (j/n)^(i-1);
    end
end

end

