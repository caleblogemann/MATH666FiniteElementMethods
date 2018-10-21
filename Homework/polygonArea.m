function [A] = polygonArea(x);
    [n, ~] = size(x);
    A = 0;
    for i = 1:(n-1)
        A = A + x(i,1)*x(i+1,2) - x(i+1,1)*x(i,2);
    end
    A = A + x(n, 1)*x(1,2) - x(1,1)*x(n,2);
    A = abs(A)/2;
end

