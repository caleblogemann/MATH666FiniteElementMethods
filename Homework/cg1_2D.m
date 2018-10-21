function [u] = cg1_2D(nodes, elements, NIN, vandermondeMatrix, gradPhiFun, bilinearForm, f)
    % Uses homogenous dirichlet boundary conditions
    % TODO: Allow for other boundary conditions
    A = sparse(NIN, NIN);
    b = zeros(NIN, 1);
    [numElements, nodesPerElement] = size(elements);
    [numNodes, ~] = size(nodes);
    % loop over elements
    for i = 1:numElements
        elementNodes = elements(i, :);
        x = nodes(elementNodes,:);
        a = polygonArea(x);
        M = vandermondeMatrix(x);
        v = eye(nodesPerElement);
        phi = (M\v)';
        grad_phi = gradPhiFun(phi);
        for j = 1:nodesPerElement
            for k = 1:nodesPerElement
                % only update for interior
                % would be different for other bc
                if (elementNodes(j) <= NIN && elementNodes(k) <= NIN)
                    A(elementNodes(j), elementNodes(k)) = ...
                        A(elementNodes(j), elementNodes(k)) + ...
                        bilinearForm(x, phi(j, :), phi(k,:), ...
                        @(x, y) grad_phi(j, x, y), @(x, y) grad_phi(k, x, y), a);
                end
            end
            % form b
            % only interior - force zero on boundary
            if (elementNodes(j) <= NIN)
                b(elementNodes(j)) = b(elementNodes(j)) + 1/nodesPerElement*a*f(x(j,1), x(j,2));
            end
        end
    end

    u = A\b;
    u = [u; zeros(numNodes - NIN, 1)];
end
