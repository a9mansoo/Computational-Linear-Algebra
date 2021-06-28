function [x, iter] = SOR(omega, A, b, x_initial, maxiter, tol)
    
    % Inital condition
    x = x_initial;
    iter = 0;
    r = norm(b - A*x);
    M = ((1/omega)*diag(diag(A))) + tril(A,-1);
    tolerance = tol*norm(b);
    
    while ( iter < maxiter && r > tolerance  )
        % Solve for x
        rhs = (M*x + (b-A*x));
        x = M \ rhs;
        % Update conditions
        iter = iter + 1;
        r = norm(b-A*x);
    end
    
end

