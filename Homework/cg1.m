function [ xi ] = cg1(Bij, Bi0, B00, Li, L0, a, b, nCells )

    A = zeros(nCells);
    for i = 1:nCells
        for j = 1:i
            if( i < nCells && j < nCells)
                B = Bij(i, j);
            elseif ( i == nCells && j < nCells)
                B = Bi0(j);
            elseif ( i == nCells && j == nCells)
                B = B00;
            end
            A(i, j) = B;
            A(j, i) = B;
        end
    end

    rhs = zeros(nCells, 1);
    for i = 1:nCells
        if( i < nCells)
            rhs(i) = Li(i);
        else
            rhs(i) = L0;
        end
    end

    xi = A\rhs;
end

