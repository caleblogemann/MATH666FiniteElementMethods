function [u] = h02_01(a, b, N)
    % this includes points at both ends
    h = (b-a)/(N-1);

    % the y position of row i is y(i)
    % this is indexed starting at 1
    x = @(j) a + (j-1)*h;
    y = @(i) a + (i-1)*h;

    % here k goes from 1 to N^2, and i and j go from 1 to N
    % point k = point (i, j)
    % this is row-wise odering, first row is 1 - N, second row N+1 - 2N, ...
    % if you are at the kth overall point and want to find the row use iFun(k), to find the colum use jFun(k)
    % if you have the row and column and want to know what point you are at use kFun(i, j)
    kFun = @(i, j) i + (j-1)*N;
    iFun = @(k) floor((k-1)/N) + 1;
    jFun = @(k) mod(k-1,N)+1;

    e = ones(N, 1);
    % D matrix gradients on omega
    % D includes i, i-1, and i+1 terms
    % matrix D for first or last row
    D1 = spdiags([-1/6*e, 4/3*e, -1/6*e], -1:1, N, N);
    % corner points
    D1(1,1) = 2/3;
    D1(end,end) = 2/3;
    % matrix D for middle rows
    D2 = spdiags([-1/3*e, 8/3*e, -1/3*e], -1:1, N, N);
    % left and right boundary points
    D2(1, 1) = 4/3;
    D2(end, end) = 4/3;

    % C matrix gradients on omega
    % C matrix includes i+m, i-m, i+m+1, i+m-1, i-m-1, and i-m+1
    C = spdiags([-1/3*e, -1/3*e, -1/3*e], -1:1, N, N);
    % left and right boundary points
    C(1, 1) = -1/6;
    C(end, end) = -1/6;

    % function on boundary
    E = spdiags();

    I = speye(N, N);
    UpperDiagonal = sparse(1:N-1, 2:N, 1, N, N);
    LowerDiagonal = sparse(2:N, 1:N-1, 1, N, N);
    diags = LowerDiagonal + UpperDiagonal;

    A = kron(I, D) + kron(diags, C);
end
