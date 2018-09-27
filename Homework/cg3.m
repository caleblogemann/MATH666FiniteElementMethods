function [ xi ] = cg3(B11ij, B12ij, B22ij, L1i, L2j, a, b, nCells )
    A11 = zeros(nCells-1);
    A12 = zeros(nCells-1);
    A22 = zeros(nCells-1);
    for i = 1:nCells-1
        for j = 1:i
            B11 = B11ij(i, j);
            A11(i, j) = B11;
            A11(j, i) = B11;
            B12 = B12ij(i, j);
            A12(i, j) = B12;
            A12(j, i) = B12;
            B22 = B22ij(i, j);
            A22(i, j) = B22;
            A22(j, i) = B22;
        end
    end
    A = [A11, A12; A12', A22];

    rhs1 = zeros(nCells-1,1);
    rhs2 = zeros(nCells-1,1);
    for i = 1:nCells-1
            rhs1(i) = L1i(i);
            rhs2(i) = L2i(i);
    end
    rhs = [rhs1; rhs2]

    xi = A\rhs;
end

