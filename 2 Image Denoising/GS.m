function [x, iter] = GS(A, b, x_initial, maxiter, tol)
    
    % Set initial conditions
    iter = 0;
    x = x_initial; 
    r = norm(b-A*x);
    M = diag(diag(A)) + tril(A, -1) ;
    U = -1 * triu(A,1);
    tolerance = tol*norm(b);
    
    while ( iter < maxiter && r > tolerance )
       % Forward solve to get x
       rhs = (b + U*x);
       x = M \ rhs;
       iter = iter + 1;
       % Calculate new norm
       r = norm(b-A*x);
    end
end

