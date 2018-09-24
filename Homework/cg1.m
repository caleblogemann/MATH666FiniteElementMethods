function [ x ] = cg1(q, f, a, b, nCells )
    h = (b - a)/nCells;
    xj = @(j) a + h*j;
    phi = @(x, j) (x - xj(j-1))/(xj(j) - xj(j-1))*(x >= xj(j-1) && x <= xj(j)) + (x - xj(j+1))/(xj(j) - xj(j+1))*(x >= xj(j) && x <= xj(j+1));
    
    
    


end

