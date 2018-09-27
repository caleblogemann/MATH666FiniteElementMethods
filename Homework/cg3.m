function [ xi ] = cg3(B11ij, B11Mp1j, B11Mp1Mp1, B12ij, B12Mp1j, B12iMp1, B12Mp1Mp1, B22ij, B22Mp1j, B22Mp1Mp1, L1i, L1Mp1, L2i, L2Mp1, a, b, M )
    A11 = zeros(M+1);
    A22 = zeros(M+1);
    for i = 1:M+1
        for j = 1:i
            if(i == M+1)
                B11 = B11Mp1j(j);
                B22 = B22Mp1j(j);
            else
                B11 = B11ij(i, j);
                B22 = B22ij(i, j);
            end
            A11(i, j) = B11;
            A11(j, i) = B11;
            A22(i, j) = B22;
            A22(j, i) = B22;
        end
    end
    A11(M+1,M+1) = B11Mp1Mp1;
    A22(M+1,M+1) = B22Mp1Mp1;
    
    A12 = zeros(M+1);
    for i = 1:M+1
        for j = 1:M+1
            if (i == M+1)
                B12 = B12Mp1j(j);
            elseif (j == M+1)
                B12 = B12iMp1(i);
            else
                B12 = B12ij(i, j);
            end
            A12(i, j) = B12;
        end
    end
    A12(M+1, M+1) = B12Mp1Mp1;
    
    A = [A11, A12; A12', A22];
    % zero out small values
    A(abs(A)<1e-10)=0;

    rhs1 = zeros(M+1,1);
    rhs2 = zeros(M+1,1);
    for i = 1:M
            rhs1(i) = L1i(i);
            rhs2(i) = L2i(i);
    end
    rhs1(M+1) = L1Mp1;
    rhs2(M+1) = L2Mp1;
    rhs = [rhs1; rhs2];

    xi = A\rhs;
end

