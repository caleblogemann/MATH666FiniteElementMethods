function [u] = uniformSquareCrouziexRaviart(f, M, a, b)
    h = (b - a)/M;

    e = ones(2*M-1,1);
    alternating = [repmat([2; 1], M-1,1); 2];
    A = spdiags([-2*e, 4*alternating, -2*e], [-1,0,1], 2*M-1, 2*M-1);
    B = 4*speye(M,M);
    C = sparse(1:2:(2*M-1), 1:M, -2*ones(M,1), 2*M-1, M);

    % location of A matrices in D
    upperDiagonal = sparse(1:M-2,2:M-1,1,M-1,M-1);
    lowerDiagonal = sparse(2:M-1,1:M-2,1,M-1,M-1);

    D = kron(speye(M-1,M-1), [A, C; C', B]) ...
        + kron(upperDiagonal, [sparse(2*M-1,3*M-1); C', sparse(M,M)]) ...
        + kron(lowerDiagonal, [sparse(2*M-1,2*M-1), C; sparse(M,3*M-1)]);

    [nRows, nCols] = size(D);
    D = [D, [sparse(nRows-M, 2*M-1); C'] ;sparse(2*M-1,nCols-M), C, A];

    x1 = a+h/2:h/2:b-h/2;
    x2 = a+h/2:h:b-h/2;
    x = [repmat([x1,x2], 1,M-1), x1];
    y1 = ones(2*M-1,1)*x1(1:2:end);
    y2 = ones(M,1)*x1(2:2:end);
    y = [cell2mat(arrayfun(@(i) [y1(:,i)',y2(:,i)'],1:M-1,'UniformOutput',false)), y1(:,end)'];
    rhs = arrayfun(@(i) f(x(i), y(i))*h^2/3, 1:length(x))';

    u = D\rhs;
end
