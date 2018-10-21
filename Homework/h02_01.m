function [u] = h02_01(a, b, N, f, g)
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
    kNodeFun = @(i, j) j + (i-1)*N;
    iNodeFun = @(k) floor((k-1)/N) + 1;
    jNodeFun = @(k) mod(k-1,N)+1;

    e = ones(N, 1);
    % D matrix gradients on omega
    % D includes i, i-1, and i+1 terms
    % matrix D for first or last row
    D1 = spdiags([-1/6*e, 4/3*e, -1/6*e], -1:1, N, N);
    % corner points
    D1(1,1) = 2/3;
    D1(end,end) = 2/3;
    % matrix D for middle rows
    D = spdiags([-1/3*e, 8/3*e, -1/3*e], -1:1, N, N);
    % left and right boundary points
    D(1, 1) = 4/3;
    D(end, end) = 4/3;

    % E matrix gradients on omega
    % E matrix includes i+m, i-m, i+m+1, i+m-1, i-m-1, and i-m+1
    E = spdiags([-1/3*e, -1/3*e, -1/3*e], -1:1, N, N);
    % left and right boundary points
    E(1, 1) = -1/6;
    E(end, end) = -1/6;

    % F matrix functions on boundaries
    % F includes i-1, i, i+1
    F1 = spdiags([h/6*e, 4*h/6*e, h/6*e], -1:1, N, N);
    F2 = sparse(N, N);
    F2(1,1) = 4*h/6;
    F2(end,end) = 4*h/6;

    G = sparse(N, N);
    G(1, 1) = h/6;
    G(end, end) = h/6;

    I = speye(N, N);
    I2 = I;
    I2(1,1) = 1/2;
    I2(end,end) = 1/2;
    I3 = I;
    I3(1, 1) = 0;
    I3(end, end) = 0;
    I4 = sparse(N, N);
    I4(1, 1) = 1;
    I4(end, end) = 1;
    UpperDiagonal = sparse(1:N-1, 2:N, 1, N, N);
    LowerDiagonal = sparse(2:N, 1:N-1, 1, N, N);
    diags = LowerDiagonal + UpperDiagonal;

    B = kron(I2, D) + kron(diags, E);
    C = kron(I4, F1) + kron(I3, F2) + kron(diags, G);
    A = B + C;

    b = zeros(N^2, 1);
    % loop over nodes
    for k = 1:N^2
        i = iNodeFun(k);
        j = jNodeFun(k);
        quadf = h^2*f(x(j), y(i));
        if (k == 1 || k == N || k==N^2 || k==N^2-N) % corner point
            b(k) = quadf/4 + h*g(x(j), y(i));
        elseif (i == 1 || j == 1 || i == N || j == N) % boundary
            b(k) = quadf/2 + h*g(x(j), y(i));
        else % interior
            b(k) = quadf;
        end
    end
    % loop over elements
    %num_elements = (N-1)^2;
    %iElemFun = @(k) floor((k-1)/(N-1)) + 1;
    %jElemFun = @(k) mod(k-1,N-1)+1;
    %for k = 1:num_elements
        %elemRow = iElemFun(k);
        %elemCol = jElemFun(k);
        %nodeRows = [elemRow, elemRow, elemRow+1, elemRow+1];
        %nodeCols = [elemCol, elemCol+1, elemCol, elemCol+1];
        %nodes = kNodeFun(nodeRows, nodeCols);
        %for j = 1:4
            %node = nodes(j);
            %nodeRow = nodeRows(j);
            %nodeCol = nodeCols(j);
            %quadf = 1/4*h^2*f(x(nodeCol), y(nodeRow));
            %b(j) = b(j) + quadf;
        %end
    %end
    u = A\b;
    %u = pcg(A, b, 1e-10, 1000);
end
