function [ xi ] = cg1(Bij, Bi0, B00, Li, L0, a, b, nCells )
    %h = (b - a)/nCells;
    %xj = @(j) a + h*j;
%     phi = @(x, j) (x - xj(j-1))/(xj(j) - xj(j-1)).*(x >= xj(j-1) & x <= xj(j)) + (x - xj(j+1))/(xj(j) - xj(j+1)).*(x >= xj(j) & x <= xj(j+1));
%     phiD = @(x, j) 1/(xj(j) - xj(j-1)).*(x >= xj(j-1) & x <= xj(j)) + 1/(xj(j) - xj(j+1)).*(x >= xj(j) & x <= xj(j+1));
    %phi = @(x, j) (x - xj(j-1))/h.*(x >= xj(j-1) & x < xj(j)) - (x - xj(j+1))/h.*(x >= xj(j) & x <= xj(j+1));
    %phiD = @(x, j) 1/h.*(x >= xj(j-1) & x < xj(j)) - 1/h.*(x >= xj(j) & x <= xj(j+1));
    %phi0 = @(x) -(x - xj(1))/h.*(x >= xj(0) & x < xj(1)) + (x - xj(nCells-1))/h.*(x >= xj(nCells-1) & x <= xj(nCells));
    %phi0D = @(x) -1/h.*(x >= xj(0) & x < xj(1)) + 1/h.*(x >= xj(nCells-1) & x <= xj(nCells));

    A = zeros(nCells);
    for i = 1:nCells
        for j = 1:i
            if( i < nCells && j < nCells)
                %Bij = B(@(x) phi(x, i), @(x) phiD(x, i), @(x) phi(x, j), @(x) phiD(x, j));
                B = Bij(i, j);
            elseif ( i == nCells && j < nCells)
                %Bij = B(@(x) phi0(x), @(x) phi0D(x), @(x) phi(x, j), @(x) phiD(x, j));
                B = Bi0(j);
            elseif ( i == nCells && j == nCells)
                %Bij = B(@(x) phi0(x), @(x) phi0D(x), @(x) phi0(x), @(x) phi0D(x));
                B = B00;
            end
            A(i, j) = B;
            A(j, i) = B;
        end
    end

    rhs = zeros(nCells, 1);
    for i = 1:nCells
        if( i < nCells)
            %rhs(i) = L(@(x) phi(x, i));
            rhs(i) = Li(i);
        else
            %rhs(i) = L(@(x) phi0(x));
            rhs(i) = L0;
        end
    end

    xi = A\rhs;

    % periodic boundary conditions
    %xi = [xi(end); xi];
end

