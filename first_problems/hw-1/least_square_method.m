function [x1, y1] = least_square_method(n, x, y)
    for i=1:1:(n+1)
        B(i) = sum(y.* (x.^ (n - i + 1)));
    end
    for i=1:1:(n + 1)
        for j=1:1:(n + 1)
            A(i,j)=sum(x.^(2 * n - i - j + 2));
        end
    end
    X = A^(-1) * B';
    x1 = (x(1):0.01:x(end))';
    y1 = zeros(length(x1), 1);
    for i=1:1:(n+1)
        y1 = y1 + X(i) * x1.^(n - i + 1);
    end
end